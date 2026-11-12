#' Pure G-Computation and Targeted Minimum Loss-Based Estimation via Survival-BART
#'
#' @description Fits nonparametric Bayesian Additive Regression Trees (BART) survival models
#'   using the Python \code{bart-survival} (\code{SurvBART}) backend (with automatic fallback
#'   to \code{BART::mc.surv.bart}). Evaluates pure G-computation and doubly-nonparametric
#'   TMLE + BART (where both the treatment initiation hazard \eqn{g(t \mid W)} and the outcome
#'   hazard \eqn{Q(t \mid A, W)} are estimated via Survival-BART, while IPCW weight assembly,
#'   logistic fluctuation targeting, and two-component EIC inference are performed in R).
#'
#' @name bart_estimators
NULL

#' Locate Python Binary and Helper Script for bart-survival
#' @keywords internal
get_bart_survival_paths <- function() {
  py_candidates <- c(
    file.path(getwd(), ".venv", "bin", "python3"),
    file.path(dirname(getwd()), ".venv", "bin", "python3"),
    Sys.which("python3")
  )
  py_bin <- py_candidates[file.exists(py_candidates) & nzchar(py_candidates)][1]

  script_candidates <- c(
    file.path(getwd(), "inst", "python", "run_bart_survival.py"),
    file.path(dirname(getwd()), "inst", "python", "run_bart_survival.py"),
    system.file("python", "run_bart_survival.py", package = "tmle.tte")
  )
  script_path <- script_candidates[file.exists(script_candidates) & nzchar(script_candidates)][1]

  list(py_bin = py_bin, script_path = script_path)
}

#' Execute Python bart-survival Backend
#' @keywords internal
run_python_surv_bart <- function(mode = c("pair", "q_only", "g_only"),
                                 X_g_train = NULL,
                                 Y_g_train = NULL,
                                 X_g_test = NULL,
                                 X_q_train = NULL,
                                 Y_q_train = NULL,
                                 X_q_test = NULL,
                                 n_train_in_test = 0L,
                                 ntree = 30L,
                                 ndpost = 100L,
                                 nskip = 100L,
                                 thin = 2L,
                                 chains = 1L,
                                 cores = 1L,
                                 seed = 99L) {
  mode <- match.arg(mode)
  paths <- get_bart_survival_paths()
  if (is.na(paths$py_bin) || is.na(paths$script_path)) {
    stop("Python bart-survival environment or inst/python/run_bart_survival.py not found.")
  }

  work_dir <- tempfile(pattern = "bart_surv_")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  meta <- list(
    mode = mode,
    trees = as.integer(ntree),
    draws = as.integer(ndpost * thin),
    tune = as.integer(nskip),
    thin = as.integer(thin),
    chains = as.integer(chains),
    cores = as.integer(cores),
    seed = as.integer(if (is.null(seed)) 99L else seed)
  )

  if (mode %in% c("pair", "g_only")) {
    X_g_train <- as.matrix(X_g_train)
    X_g_test <- as.matrix(X_g_test)
    meta$n_g_train <- nrow(X_g_train)
    meta$p_g <- ncol(X_g_train)
    meta$n_g_test <- nrow(X_g_test)
    writeBin(as.double(t(X_g_train)), file.path(work_dir, "x_g_train.bin"), size = 8)
    writeBin(as.double(Y_g_train), file.path(work_dir, "y_g_train.bin"), size = 8)
    writeBin(as.double(t(X_g_test)), file.path(work_dir, "x_g_test.bin"), size = 8)
  }

  if (mode %in% c("pair", "q_only")) {
    X_q_train <- as.matrix(X_q_train)
    X_q_test <- as.matrix(X_q_test)
    meta$n_q_train <- nrow(X_q_train)
    meta$p_q <- ncol(X_q_train)
    meta$n_q_test <- nrow(X_q_test)
    meta$n_train_in_test <- as.integer(n_train_in_test)
    writeBin(as.double(t(X_q_train)), file.path(work_dir, "x_q_train.bin"), size = 8)
    writeBin(as.double(Y_q_train), file.path(work_dir, "y_q_train.bin"), size = 8)
    writeBin(as.double(t(X_q_test)), file.path(work_dir, "x_q_test.bin"), size = 8)
  }

  # Write simple JSON metadata without external package dependency
  kv <- vapply(names(meta), function(k) {
    val <- meta[[k]]
    if (is.character(val)) sprintf('"%s": "%s"', k, val) else sprintf('"%s": %d', k, as.integer(val))
  }, character(1))
  writeLines(paste0("{", paste(kv, collapse = ", "), "}"), file.path(work_dir, "meta.json"))

  res_code <- system2(paths$py_bin, args = shQuote(c(paths$script_path, work_dir)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(res_code, "status")) && attr(res_code, "status") != 0) {
    stop("Python bart-survival failed:\n", paste(res_code, collapse = "\n"))
  }

  out <- list()
  if (mode %in% c("pair", "g_only")) {
    out$g_pred_mean <- readBin(
      file.path(work_dir, "g_pred_mean.bin"),
      what = "double",
      n = meta$n_g_test,
      size = 8
    )
  }

  if (mode %in% c("pair", "q_only")) {
    if (n_train_in_test > 0L) {
      out$q_train_mean <- readBin(
        file.path(work_dir, "q_train_mean.bin"),
        what = "double",
        n = n_train_in_test,
        size = 8
      )
    }
    n_cf_cols <- meta$n_q_test - n_train_in_test
    cf_file <- file.path(work_dir, "q_cf_draws.bin")
    n_floats <- file.info(cf_file)$size / 8
    M_q <- as.integer(round(n_floats / n_cf_cols))
    raw_cf <- readBin(cf_file, what = "double", n = M_q * n_cf_cols, size = 8)
    out$q_cf_draws <- matrix(raw_cf, nrow = M_q, ncol = n_cf_cols, byrow = TRUE)
  }

  out
}

#' Summarize Posterior Counterfactual Hazard Draws for G-COMP (BART)
#' @keywords internal
summarize_bart_gcomp_draws <- function(q_cf_draws,
                                       n_patients,
                                       tau,
                                       w_mat = NULL,
                                       k_vals = NULL) {
  n_cf <- n_patients * tau
  M <- nrow(q_cf_draws)

  P0 <- q_cf_draws[, seq_len(n_cf), drop = FALSE]
  P1 <- q_cf_draws[, (n_cf + 1L):(2L * n_cf), drop = FALSE]

  dim(P0) <- c(M, tau, n_patients)
  dim(P1) <- c(M, tau, n_patients)

  log_1m_0 <- log(pmax(1 - P0, 1e-6))
  log_1m_1 <- log(pmax(1 - P1, 1e-6))

  # Control arm: always untreated (P0) across t = 1..tau
  log_ctrl <- log_1m_0
  if (tau > 1L) {
    for (tt in 2:tau) {
      log_ctrl[, tt, ] <- log_ctrl[, tt - 1L, ] + log_ctrl[, tt, ]
    }
  }
  S_ctrl_draws <- exp(log_ctrl)
  marg_ctrl <- rowMeans(S_ctrl_draws, dims = 2L)

  if (is.null(w_mat) || is.null(k_vals)) {
    # Legacy fallback when q_cf_draws[, (n_cf+1):(2*n_cf)] already encodes a single regime
    log_treat <- log_1m_1
    if (tau > 1L) {
      for (tt in 2:tau) {
        log_treat[, tt, ] <- log_treat[, tt - 1L, ] + log_treat[, tt, ]
      }
    }
    S_treat_draws <- exp(log_treat)
    marg_treat <- rowMeans(S_treat_draws, dims = 2L)
  } else {
    # Exact analytical integration over initiation policy weights w_mat (n_patients x n_k)
    marg_treat <- matrix(0, nrow = M, ncol = tau)
    for (idx in seq_along(k_vals)) {
      wk <- w_mat[, idx]
      if (!any(wk > 0)) next
      k <- as.integer(k_vals[idx])
      log_k <- log_1m_0
      if (k <= 1L) {
        log_k <- log_1m_1
      } else if (k <= tau) {
        log_k[, k:tau, ] <- log_1m_1[, k:tau, ]
      }
      if (tau > 1L) {
        for (tt in 2:tau) {
          log_k[, tt, ] <- log_k[, tt - 1L, ] + log_k[, tt, ]
        }
      }
      S_k <- exp(log_k)
      dim(S_k) <- c(M * tau, n_patients)
      marg_treat <- marg_treat + matrix(as.vector(S_k %*% wk), nrow = M, ncol = tau) / n_patients
    }
  }

  rd_draws <- marg_treat - marg_ctrl
  rmstd_draws <- rowSums(rd_draws)

  mean_s0 <- colMeans(marg_ctrl)
  mean_s1 <- colMeans(marg_treat)
  mean_rd <- colMeans(rd_draws)
  sd_rd <- apply(rd_draws, 2, stats::sd)

  cri_lower <- apply(rd_draws, 2, stats::quantile, probs = 0.025, names = FALSE)
  cri_upper <- apply(rd_draws, 2, stats::quantile, probs = 0.975, names = FALSE)

  # Simultaneous credible band over active event horizon
  valid_cols <- which(sd_rd > 0.002)
  if (length(valid_cols) < 5L) valid_cols <- seq(from = min(5L, tau), to = tau)
  std_dev <- sweep(abs(sweep(rd_draws, 2, mean_rd, "-")), 2, pmax(sd_rd, 1e-6), "/")
  sup_dev <- apply(std_dev[, valid_cols, drop = FALSE], 1, max)
  crit_simul <- stats::quantile(sup_dev, probs = 0.95, names = FALSE)

  band_lower <- mean_rd - crit_simul * sd_rd
  band_upper <- mean_rd + crit_simul * sd_rd

  # Cumulative HR posterior draws per time t
  log_hr_draws <- log(
    log(pmax(pmin(marg_treat, 1 - 1e-6), 1e-6)) /
      log(pmax(pmin(marg_ctrl, 1 - 1e-6), 1e-6))
  )
  hr_cum <- log(pmax(pmin(mean_s1, 1 - 1e-6), 1e-6)) / log(pmax(pmin(mean_s0, 1 - 1e-6), 1e-6))
  se_log_hr <- apply(log_hr_draws, 2, stats::sd, na.rm = TRUE)
  hr_ci_lower <- exp(log(pmax(hr_cum, 1e-6)) - 1.96 * se_log_hr)
  hr_ci_upper <- exp(log(pmax(hr_cum, 1e-6)) + 1.96 * se_log_hr)

  res_table <- tibble::tibble(
    t = seq_len(tau),
    control_no_tx = mean_s0,
    treat_tx = mean_s1,
    rd = mean_rd,
    se_RD = sd_rd,
    cri_lower = cri_lower,
    cri_upper = cri_upper,
    cri_width = cri_upper - cri_lower,
    band_lower = band_lower,
    band_upper = band_upper,
    band_width = band_upper - band_lower,
    hr_cum = hr_cum,
    se_log_hr = se_log_hr,
    hr_ci_lower = hr_ci_lower,
    hr_ci_upper = hr_ci_upper
  )

  list(
    results = res_table,
    rd_draws = rd_draws,
    rmstd = mean(rmstd_draws),
    se_rmstd = stats::sd(rmstd_draws),
    rmstd_ci_lower = stats::quantile(rmstd_draws, probs = 0.025, names = FALSE),
    rmstd_ci_upper = stats::quantile(rmstd_draws, probs = 0.975, names = FALSE),
    crit_val_simul = crit_simul
  )
}

#' Assemble IPCW Panel and Counterfactual Clever Covariates from BART-g Predictions
#' @keywords internal
compute_bart_ipcw_and_cf <- function(base_ord,
                                     clones,
                                     g_grid_pred,
                                     covars,
                                     grace_period = 10,
                                     tau = 30,
                                     trim_quantile = c(0.001, 0.999),
                                     max_weight = 10,
                                     tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                                     align_grace_weights = TRUE) {
  tx_strategy <- match.arg(tx_strategy)
  if (!align_grace_weights && tx_strategy != "natural_grace") {
    tx_strategy <- "natural_grace"
  }

  n_patients <- nrow(base_ord)
  n_grace_int <- max(as.integer(grace_period), 1L)
  # g_grid_pred is ordered by sorted patient id (base_ord), and within patient by t = 1..n_grace_int
  g_mat <- matrix(bound_precision(g_grid_pred), nrow = n_grace_int, ncol = n_patients)
  colnames(g_mat) <- as.character(base_ord$id)

  # Marginal stabilization numerators (matching compute_ipcw_weights)
  obs_panel <- clones %>%
    dplyr::filter(clone_arm == 1L) %>%
    dplyr::mutate(
      tmax = pmin(dplyr::if_else(tx_ever == 1L, t_switch, time), time, grace_period),
      n_int = pmax(ceiling(tmax), 1L)
    ) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_tx = dplyr::if_else(tx_ever == 1L & t == t_switch & time >= t_switch & t <= grace_period, 1L, 0L)
    )

  treat_day10_data <- clones %>%
    dplyr::filter(arm == "treat_tx", time > grace_period) %>%
    dplyr::mutate(took_tx = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))

  df_t <- min(3L, max(n_grace_int - 1L, 1L))
  num_ctrl_form <- stats::as.formula(paste("event_tx ~ splines::ns(t, df =", df_t, ")"))
  m_num_ctrl <- suppressWarnings(stats::glm(num_ctrl_form, family = stats::binomial(), data = obs_panel))
  m_num_treat <- suppressWarnings(stats::glm(took_tx ~ 1, family = stats::binomial(), data = treat_day10_data))

  grace_surv <- compute_grace_survival_grid(
    data_patients = base_ord,
    covars = covars,
    grace_period = grace_period,
    dt = 1,
    m_denom_ctrl = NULL,
    m_num_ctrl = m_num_ctrl,
    use_sl = FALSE,
    g_mat_override = g_mat
  )
  p_star <- compute_empirical_complier_probs(
    clones = clones,
    grace_period = grace_period,
    n_grace_int = grace_surv$n_grace_int
  )
  p_treat_denom_vec <- bound_precision(1 - grace_surv$S_denom_mat[, grace_surv$n_grace_int])
  p_treat_num_vec <- bound_precision(stats::predict(m_num_treat, newdata = base_ord, type = "response"))

  # Expand clones into longitudinal person-time panel
  panel <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      t0 = t - 1L,
      t1 = t,
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > t0 & time_cc <= t1, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
    ) %>%
    dplyr::arrange(id, clone_arm, t)

  panel <- compute_arm_ipcw_weights(
    df_panel = panel,
    grace_surv = grace_surv,
    covars = covars,
    grace_period = grace_period,
    m_denom_treat = NULL,
    m_num_treat = m_num_treat,
    use_sl = FALSE,
    tx_strategy = tx_strategy,
    p_star = p_star,
    p_treat_denom_vec = p_treat_denom_vec,
    p_treat_num_vec = p_treat_num_vec,
    is_cf = FALSE
  ) %>%
    dplyr::mutate(
      w_ipcw_time = w_stab,
      w_unstab_time = w_unstab
    )

  q_lim <- compute_arm_trim_bounds(panel$w_ipcw_time, panel$arm, probs = trim_quantile, max_weight = max_weight)
  q_lim_unstab <- compute_arm_trim_bounds(
    panel$w_unstab_time,
    panel$arm,
    probs = trim_quantile,
    max_weight = if (is.finite(max_weight)) 2 * max_weight else Inf
  )

  panel_w <- panel %>%
    dplyr::mutate(
      w_trim = apply_arm_weight_trimming(w_ipcw_time, arm, q_lim),
      w_unstab_trim = apply_arm_weight_trimming(w_unstab_time, arm, q_lim_unstab),
      post_grace = dplyr::if_else(t > grace_period, 1L, 0L)
    )

  strat <- compute_tx_strategy_weights(
    grace_surv = grace_surv,
    tx_strategy = tx_strategy,
    p_star = p_star
  )
  cf_w <- compute_cf_weight_matrices(
    grace_surv = grace_surv,
    data_patients = base_ord,
    covars = covars,
    grace_period = grace_period,
    tau = tau,
    tx_strategy = tx_strategy,
    k_vals = strat$k_vals,
    q_lim = q_lim,
    m_denom_treat = NULL,
    m_num_treat = m_num_treat,
    use_sl = FALSE,
    p_star = p_star,
    p_treat_denom_vec = p_treat_denom_vec,
    p_treat_num_vec = p_treat_num_vec
  )

  list(
    panel_w = panel_w,
    grace_surv = grace_surv,
    strat = strat,
    cf_H_ctrl = cf_w$cf_H_ctrl,
    cf_H_treat_list = cf_w$cf_H_treat_list,
    q_lim = q_lim,
    tx_strategy = tx_strategy,
    p_star = p_star
  )
}

#' Fit Both G-COMP (BART) and TMLE + BART in a Single Shared Python bart-survival Pass
#'
#' @description Trains \code{SurvBART} for the treatment initiation hazard \eqn{g(t \mid W)}
#'   and outcome hazard \eqn{Q(t \mid A, W)} via \code{bart-survival}, then evaluates both
#'   \code{G-COMP (BART)} and \code{TMLE + BART} in R.
#'
#' @param data Raw standardized data frame (after \code{prepare_tte_data}).
#' @param clones Cloned data frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param ntree Number of trees in the BART ensemble (default 50, matching \code{SurvBART.preset_standard()}).
#' @param ndpost Number of retained posterior draws per chain (default 250, i.e., 1,000 total draws with thin = 4).
#' @param nskip Number of burn-in tune iterations per chain (default 250).
#' @param thin Thinning factor (default 4).
#' @param chains Number of MCMC chains in Python \code{SurvBART} (default 4).
#' @param cores Number of worker cores inside Python \code{SurvBART} (default 1).
#' @param n_boot Multiplier bootstrap replicates (default 5000).
#' @param tx_strategy Treatment-arm grace-period policy: \code{"conditional_complier"} (default),
#'   \code{"empirical_complier"}, or \code{"natural_grace"}.
#' @param align_grace_weights Legacy alias: if \code{FALSE}, overrides \code{tx_strategy} to
#'   \code{"natural_grace"}.
#' @param max_weight Hard upper cap applied after quantile trimming on stabilized IPCW weights
#'   (default 10; \code{Inf} disables hard capping).
#' @param seed Optional random seed.
#' @return A list with elements \code{gcomp} (output matching \code{fit_bart_gcomp}) and
#'   \code{tmle} (output matching \code{fit_bart_tmle}).
#' @export
fit_bart_pair <- function(data,
                          clones,
                          covars = c("w1", "w2"),
                          grace_period = 10,
                          tau = 30,
                          ntree = 50,
                          ndpost = 250,
                          nskip = 250,
                          thin = 4,
                          chains = 4L,
                          cores = 1L,
                          n_boot = 5000,
                          tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                          align_grace_weights = TRUE,
                          max_weight = 10,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tx_strategy <- match.arg(tx_strategy)
  if (!align_grace_weights && tx_strategy != "natural_grace") {
    tx_strategy <- "natural_grace"
  }

  base_ord <- data %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::arrange(id)
  n_patients <- nrow(base_ord)
  n_grace_int <- max(as.integer(grace_period), 1L)

  # 1. Build grace-period treatment initiation training panel for BART-g
  obs_panel_g <- clones %>%
    dplyr::filter(clone_arm == 1L) %>%
    dplyr::mutate(
      tmax = pmin(dplyr::if_else(tx_ever == 1L, t_switch, time), time, grace_period),
      n_int = pmax(ceiling(tmax), 1L)
    ) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_tx = dplyr::if_else(tx_ever == 1L & t == t_switch & time >= t_switch & t <= grace_period, 1L, 0L)
    )

  X_g_train <- obs_panel_g %>%
    dplyr::select(t, dplyr::all_of(covars)) %>%
    as.matrix()
  Y_g_train <- obs_panel_g$event_tx

  # Prediction grid for BART-g: every sorted patient i across t = 1..n_grace_int
  g_grid <- base_ord %>%
    dplyr::select(id, dplyr::all_of(covars)) %>%
    tidyr::uncount(weights = n_grace_int, .id = "t")
  X_g_test <- g_grid %>%
    dplyr::select(t, dplyr::all_of(covars)) %>%
    as.matrix()

  # 2. Build cloned person-time panel for BART-Q (ordered by id, clone_arm, t)
  panel_q <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > (t - 1L) & time_cc <= t, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
    ) %>%
    dplyr::arrange(id, clone_arm, t)

  X_q_train <- panel_q %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()
  Y_q_train <- panel_q$event_cc_day
  n_train <- nrow(X_q_train)

  # 3. Build pure untreated (is_treated = 0) and pure treated (is_treated = 1) grids on sorted patients
  cf_base <- base_ord %>%
    dplyr::select(id, dplyr::all_of(covars)) %>%
    tidyr::uncount(weights = tau, .id = "t")

  cf_0 <- cf_base %>%
    dplyr::mutate(is_treated = 0L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  cf_1 <- cf_base %>%
    dplyr::mutate(is_treated = 1L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  X_q_test <- rbind(X_q_train, cf_0, cf_1)

  # 4. Run Python bart-survival for both g and Q
  t0_py <- proc.time()["elapsed"]
  py_out <- run_python_surv_bart(
    mode = "pair",
    X_g_train = X_g_train,
    Y_g_train = Y_g_train,
    X_g_test = X_g_test,
    X_q_train = X_q_train,
    Y_q_train = Y_q_train,
    X_q_test = X_q_test,
    n_train_in_test = n_train,
    ntree = ntree,
    ndpost = ndpost,
    nskip = nskip,
    thin = thin,
    chains = chains,
    cores = cores,
    seed = if (is.null(seed)) 99L else as.integer(seed)
  )
  t_py_elapsed <- as.numeric(proc.time()["elapsed"] - t0_py)

  # 5. Assemble IPCW weights, policy mixture weights w_k(W_i), and cf_H matrices from BART-g
  t0_tmle <- proc.time()["elapsed"]
  ipcw_bart <- compute_bart_ipcw_and_cf(
    base_ord = base_ord,
    clones = clones,
    g_grid_pred = py_out$g_pred_mean,
    covars = covars,
    grace_period = grace_period,
    tau = tau,
    max_weight = max_weight,
    tx_strategy = tx_strategy,
    align_grace_weights = align_grace_weights
  )

  # 6. Evaluate G-COMP (BART) in R via analytical integration across k_vals
  t0_gcomp <- proc.time()["elapsed"]
  gcomp_res <- summarize_bart_gcomp_draws(
    q_cf_draws = py_out$q_cf_draws,
    n_patients = n_patients,
    tau = tau,
    w_mat = ipcw_bart$strat$w_mat,
    k_vals = ipcw_bart$strat$k_vals
  )
  gcomp_res$tx_strategy <- ipcw_bart$tx_strategy
  gcomp_res$elapsed_sec <- 0.5 * t_py_elapsed + as.numeric(proc.time()["elapsed"] - t0_gcomp)

  # 7. Evaluate TMLE + BART in R (BART-Q initial hazards + logistic targeting + analytical integration + EIC)
  panel_w <- ipcw_bart$panel_w
  n_cf <- nrow(cf_base)
  Q_init_train <- py_out$q_train_mean
  Q_init_0 <- colMeans(py_out$q_cf_draws[, seq_len(n_cf), drop = FALSE])
  Q_init_1 <- colMeans(py_out$q_cf_draws[, (n_cf + 1L):(2L * n_cf), drop = FALSE])

  targ_res <- target_hazard(
    panel_w = panel_w,
    pred_Q_init = Q_init_train,
    grace_period = grace_period
  )
  panel_fluc <- targ_res$panel_w
  eps_ctrl <- targ_res$eps_ctrl
  eps_treat_pre <- targ_res$eps_treat_pre
  eps_treat_post <- targ_res$eps_treat_post

  Q0_logit_mat <- stats::qlogis(bound_precision(matrix(Q_init_0, nrow = n_patients, ncol = tau, byrow = TRUE)))
  Q1_logit_mat <- stats::qlogis(bound_precision(matrix(Q_init_1, nrow = n_patients, ncol = tau, byrow = TRUE)))

  # Control arm fluctuated hazards & cumulative survival (n_patients x tau)
  Q_star_ctrl <- stats::plogis(Q0_logit_mat + eps_ctrl * ipcw_bart$cf_H_ctrl)
  S_mat_ctrl <- 1 - Q_star_ctrl
  if (tau > 1L) {
    for (tt in 2:tau) {
      S_mat_ctrl[, tt] <- S_mat_ctrl[, tt - 1L] * S_mat_ctrl[, tt]
    }
  }

  # Treated arm: exact analytical expectation sum_k w_k(W_i) S_i^*(t | T_A = k, W_i)
  eps_treat_vec <- ifelse(seq_len(tau) <= grace_period, eps_treat_pre, eps_treat_post)
  eps_treat_mat <- matrix(eps_treat_vec, nrow = n_patients, ncol = tau, byrow = TRUE)
  S_mat_treat <- matrix(0, nrow = n_patients, ncol = tau)
  strat <- ipcw_bart$strat
  for (idx in seq_along(strat$k_vals)) {
    wk <- strat$w_mat[, idx]
    if (!any(wk > 0)) next
    k <- as.integer(strat$k_vals[idx])
    Q_init_k <- Q0_logit_mat
    if (k <= 1L) {
      Q_init_k <- Q1_logit_mat
    } else if (k <= tau) {
      Q_init_k[, k:tau] <- Q1_logit_mat[, k:tau, drop = FALSE]
    }
    Q_star_k <- stats::plogis(Q_init_k + eps_treat_mat * ipcw_bart$cf_H_treat_list[[idx]])
    S_k <- 1 - Q_star_k
    if (tau > 1L) {
      for (tt in 2:tau) {
        S_k[, tt] <- S_k[, tt - 1L] * S_k[, tt]
      }
    }
    S_mat_treat <- S_mat_treat + sweep(S_k, 1L, wk, "*")
  }

  mean_s0 <- colMeans(S_mat_ctrl)
  mean_s1 <- colMeans(S_mat_treat)
  mean_rd <- mean_s1 - mean_s0

  cf_ctrl_df <- cf_base %>%
    dplyr::mutate(arm = "control_no_tx", indiv_survival = as.vector(t(S_mat_ctrl)))
  cf_treat_df <- cf_base %>%
    dplyr::mutate(arm = "treat_tx", indiv_survival = as.vector(t(S_mat_treat)))
  cf_data <- dplyr::bind_rows(cf_ctrl_df, cf_treat_df)

  tmle_surv <- tibble::tibble(
    t = rep(seq_len(tau), 2),
    arm = c(rep("control_no_tx", tau), rep("treat_tx", tau)),
    survival = c(mean_s0, mean_s1)
  )
  res_tmle <- tibble::tibble(
    t = seq_len(tau),
    control_no_tx = mean_s0,
    treat_tx = mean_s1,
    rd = mean_rd
  )

  inf_res <- compute_eic_inference(
    panel_w = panel_fluc,
    cf_data = cf_data,
    tmle_surv = tmle_surv,
    res_tmle = res_tmle,
    eps_ctrl = eps_ctrl,
    eps_treat_pre = eps_treat_pre,
    eps_treat_post = eps_treat_post,
    n_boot = n_boot,
    alpha = 0.05,
    seed = seed
  )

  tmle_out <- list(
    results = inf_res$res_tmle,
    eps = c(eps_ctrl = eps_ctrl, eps_treat_pre = eps_treat_pre, eps_treat_post = eps_treat_post),
    crit_val_simul = inf_res$crit_val_simul,
    se_rmstd = inf_res$se_rmstd,
    m_fluc = targ_res$m_fluctuate,
    tx_strategy = ipcw_bart$tx_strategy,
    elapsed_sec = 0.5 * t_py_elapsed + as.numeric(proc.time()["elapsed"] - t0_tmle)
  )

  list(
    gcomp = gcomp_res,
    tmle = tmle_out
  )
}

#' Pure G-Computation via Nonparametric Survival-BART
#'
#' @param data Raw standardized data frame.
#' @param clones Cloned data frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param ntree Number of trees in the BART ensemble (default 50, matching \code{SurvBART.preset_standard()}).
#' @param ndpost Number of retained posterior MCMC draws per chain (default 250, i.e., 1,000 total draws with thin = 4).
#' @param nskip Number of MCMC burn-in iterations per chain (default 250).
#' @param thin Thinning factor (default 4).
#' @param chains Number of MCMC chains in Python \code{SurvBART} (default 4).
#' @param cores Number of worker cores inside Python \code{SurvBART} (default 1).
#' @param mc.cores Unused when Python \code{bart-survival} backend is active.
#' @param p_form Optional treatment-initiation propensity formula used to compute policy weights
#'   when \code{tx_strategy != "empirical_complier"}. Defaults to natural cubic spline in \code{t} plus \code{covars}.
#' @param tx_strategy Treatment-arm grace-period policy: \code{"conditional_complier"} (default),
#'   \code{"empirical_complier"}, or \code{"natural_grace"}.
#' @param align_grace_weights Legacy alias: if \code{FALSE}, overrides \code{tx_strategy} to
#'   \code{"natural_grace"}.
#' @param seed Optional random seed.
#' @export
fit_bart_gcomp <- function(data,
                           clones,
                           covars = c("w1", "w2"),
                           grace_period = 10,
                           tau = 30,
                           ntree = 50,
                           ndpost = 250,
                           nskip = 250,
                           thin = 4,
                           chains = 4L,
                           cores = 1L,
                           mc.cores = 1,
                           p_form = NULL,
                           tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                           align_grace_weights = TRUE,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tx_strategy <- match.arg(tx_strategy)
  if (!align_grace_weights && tx_strategy != "natural_grace") {
    tx_strategy <- "natural_grace"
  }
  t0_all <- proc.time()["elapsed"]

  base_ord <- data %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::arrange(id)
  n_patients <- nrow(base_ord)
  n_grace_int <- max(as.integer(grace_period), 1L)

  panel_bart <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > (t - 1L) & time_cc <= t, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
    ) %>%
    dplyr::arrange(id, clone_arm, t)

  X_train <- panel_bart %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()
  Y_train <- panel_bart$event_cc_day

  cf_base <- base_ord %>%
    dplyr::select(id, dplyr::all_of(covars)) %>%
    tidyr::uncount(weights = tau, .id = "t")

  cf_0 <- cf_base %>%
    dplyr::mutate(is_treated = 0L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  cf_1 <- cf_base %>%
    dplyr::mutate(is_treated = 1L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  X_test <- rbind(cf_0, cf_1)

  py_out <- run_python_surv_bart(
    mode = "q_only",
    X_q_train = X_train,
    Y_q_train = Y_train,
    X_q_test = X_test,
    n_train_in_test = 0L,
    ntree = ntree,
    ndpost = ndpost,
    nskip = nskip,
    thin = thin,
    chains = chains,
    cores = cores,
    seed = if (is.null(seed)) 99L else as.integer(seed)
  )

  # Policy mixture weights w_k(W_i) on sorted patients
  obs_panel <- clones %>%
    dplyr::filter(clone_arm == 1L) %>%
    dplyr::mutate(
      tmax = pmin(dplyr::if_else(tx_ever == 1L, t_switch, time), time, grace_period),
      n_int = pmax(ceiling(tmax), 1L)
    ) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_tx = dplyr::if_else(tx_ever == 1L & t == t_switch & time >= t_switch & t <= grace_period, 1L, 0L)
    )

  df_t <- min(3L, max(n_grace_int - 1L, 1L))
  p_form_use <- if (is.null(p_form)) {
    stats::as.formula(paste("event_tx ~ splines::ns(t, df =", df_t, ") +", paste(covars, collapse = " + ")))
  } else {
    p_form
  }
  num_form <- stats::as.formula(paste("event_tx ~ splines::ns(t, df =", df_t, ")"))
  m_denom_ctrl <- suppressWarnings(stats::glm(p_form_use, family = stats::binomial(), data = obs_panel))
  m_num_ctrl <- suppressWarnings(stats::glm(num_form, family = stats::binomial(), data = obs_panel))

  grace_surv <- compute_grace_survival_grid(
    data_patients = base_ord,
    covars = covars,
    grace_period = grace_period,
    dt = 1,
    m_denom_ctrl = m_denom_ctrl,
    m_num_ctrl = m_num_ctrl,
    use_sl = FALSE
  )
  p_star <- compute_empirical_complier_probs(
    clones = clones,
    grace_period = grace_period,
    n_grace_int = grace_surv$n_grace_int
  )
  strat <- compute_tx_strategy_weights(
    grace_surv = grace_surv,
    tx_strategy = tx_strategy,
    p_star = p_star
  )

  res <- summarize_bart_gcomp_draws(
    q_cf_draws = py_out$q_cf_draws,
    n_patients = n_patients,
    tau = tau,
    w_mat = strat$w_mat,
    k_vals = strat$k_vals
  )
  res$tx_strategy <- tx_strategy
  res$elapsed_sec <- as.numeric(proc.time()["elapsed"] - t0_all)
  res
}

#' Targeted Minimum Loss-Based Estimation via Survival-BART for Both g and Q (TMLE + BART)
#'
#' @param data Raw standardized data frame.
#' @param clones Cloned data frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param ntree Number of trees in the BART ensemble (default 50, matching \code{SurvBART.preset_standard()}).
#' @param ndpost Number of retained posterior MCMC draws per chain (default 250, i.e., 1,000 total draws with thin = 4).
#' @param nskip Number of MCMC burn-in iterations per chain (default 250).
#' @param thin Thinning factor (default 4).
#' @param chains Number of MCMC chains in Python \code{SurvBART} (default 4).
#' @param cores Number of worker cores inside Python \code{SurvBART} (default 1).
#' @param mc.cores Unused when Python \code{bart-survival} backend is active.
#' @param n_boot Multiplier bootstrap replicates (default 5000).
#' @param tx_strategy Treatment-arm grace-period policy: \code{"conditional_complier"} (default),
#'   \code{"empirical_complier"}, or \code{"natural_grace"}.
#' @param align_grace_weights Legacy alias: if \code{FALSE}, overrides \code{tx_strategy} to
#'   \code{"natural_grace"}.
#' @param max_weight Hard upper cap applied after quantile trimming on stabilized IPCW weights
#'   (default 10; \code{Inf} disables hard capping).
#' @param seed Optional random seed.
#' @export
fit_bart_tmle <- function(data,
                          clones,
                          covars = c("w1", "w2"),
                          grace_period = 10,
                          tau = 30,
                          ntree = 50,
                          ndpost = 250,
                          nskip = 250,
                          thin = 4,
                          chains = 4L,
                          cores = 1L,
                          mc.cores = 1,
                          n_boot = 5000,
                          tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                          align_grace_weights = TRUE,
                          max_weight = 10,
                          seed = NULL) {
  pair_res <- fit_bart_pair(
    data = data,
    clones = clones,
    covars = covars,
    grace_period = grace_period,
    tau = tau,
    ntree = ntree,
    ndpost = ndpost,
    nskip = nskip,
    thin = thin,
    chains = chains,
    cores = cores,
    n_boot = n_boot,
    tx_strategy = tx_strategy,
    align_grace_weights = align_grace_weights,
    max_weight = max_weight,
    seed = seed
  )
  pair_res$tmle
}
