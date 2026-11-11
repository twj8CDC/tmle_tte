#!/usr/bin/env python3
"""
Python bridge for bart-survival (SurvBART) discrete-time hazard estimation.
Fits treatment initiation hazard g(t | W) and/or outcome hazard Q(t | A, W)
using bart-survival v0.2.5 and returns posterior hazard predictions to R.
"""
import os
import sys
import json

# Configure sandbox-safe cache and local Manager prior to importing pymc_bart
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
os.environ.setdefault("PYTENSOR_FLAGS", "base_compiledir=/tmp/pytensor")

import multiprocessing


class _LocalManager:
    def list(self, iterable=()):
        return list(iterable)

    def dict(self, mapping=()):
        return dict(mapping)


multiprocessing.Manager = lambda: _LocalManager()

import numpy as np
from bart_survival import SurvBART


def fit_and_predict(x_train, y_train, x_test, trees, draws, tune, thin, seed):
    np.random.seed(int(seed))
    model = SurvBART(
        trees=int(trees),
        draws=int(draws),
        tune=int(tune),
        chains=1,
        cores=1,
        thin=int(thin),
        split_method="empirical",
    )
    model.fit(y_train, x_train, progressbar=False, random_seed=int(seed))
    preds = np.asarray(model.predict(x_test), dtype=np.float64)
    return preds


def main():
    if len(sys.argv) < 2:
        raise ValueError("Usage: run_bart_survival.py <work_dir>")
    work_dir = sys.argv[1]
    with open(os.path.join(work_dir, "meta.json"), "r") as f:
        meta = json.load(f)

    mode = meta.get("mode", "pair")
    trees = meta.get("trees", 30)
    draws = meta.get("draws", 200)
    tune = meta.get("tune", 100)
    thin = meta.get("thin", 2)
    seed = meta.get("seed", 99)

    out_meta = {}

    if mode in ("pair", "g_only"):
        n_g_train = int(meta["n_g_train"])
        p_g = int(meta["p_g"])
        n_g_test = int(meta["n_g_test"])

        x_g_train = np.fromfile(
            os.path.join(work_dir, "x_g_train.bin"), dtype=np.float64
        ).reshape((n_g_train, p_g))
        y_g_train = np.fromfile(
            os.path.join(work_dir, "y_g_train.bin"), dtype=np.float64
        )
        x_g_test = np.fromfile(
            os.path.join(work_dir, "x_g_test.bin"), dtype=np.float64
        ).reshape((n_g_test, p_g))

        g_draws = fit_and_predict(
            x_g_train,
            y_g_train,
            x_g_test,
            trees=trees,
            draws=draws,
            tune=tune,
            thin=thin,
            seed=seed,
        )
        g_mean = np.mean(g_draws, axis=0).astype(np.float64)
        g_mean.tofile(os.path.join(work_dir, "g_pred_mean.bin"))
        out_meta["M_g"] = int(g_draws.shape[0])

    if mode in ("pair", "q_only"):
        n_q_train = int(meta["n_q_train"])
        p_q = int(meta["p_q"])
        n_q_test = int(meta["n_q_test"])
        n_train_in_test = int(meta.get("n_train_in_test", 0))

        x_q_train = np.fromfile(
            os.path.join(work_dir, "x_q_train.bin"), dtype=np.float64
        ).reshape((n_q_train, p_q))
        y_q_train = np.fromfile(
            os.path.join(work_dir, "y_q_train.bin"), dtype=np.float64
        )
        x_q_test = np.fromfile(
            os.path.join(work_dir, "x_q_test.bin"), dtype=np.float64
        ).reshape((n_q_test, p_q))

        q_draws = fit_and_predict(
            x_q_train,
            y_q_train,
            x_q_test,
            trees=trees,
            draws=draws,
            tune=tune,
            thin=thin,
            seed=seed + 17,
        )
        out_meta["M_q"] = int(q_draws.shape[0])

        if n_train_in_test > 0:
            q_train_mean = np.mean(q_draws[:, :n_train_in_test], axis=0).astype(
                np.float64
            )
            q_train_mean.tofile(os.path.join(work_dir, "q_train_mean.bin"))
            q_cf_draws = np.ascontiguousarray(
                q_draws[:, n_train_in_test:], dtype=np.float64
            )
            q_cf_draws.tofile(os.path.join(work_dir, "q_cf_draws.bin"))
        else:
            q_cf_draws = np.ascontiguousarray(q_draws, dtype=np.float64)
            q_cf_draws.tofile(os.path.join(work_dir, "q_cf_draws.bin"))

    with open(os.path.join(work_dir, "out_meta.json"), "w") as f:
        json.dump(out_meta, f)


if __name__ == "__main__":
    main()
