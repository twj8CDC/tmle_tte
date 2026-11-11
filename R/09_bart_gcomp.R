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
summarize_bart_gcomp_draws <- function(q_cf_draws, n_patients, tau) {
  n_cf <- n_patients * tau
  M <- nrow(q_cf_draws)

  P_ctrl <- q_cf_draws[, seq_len(n_cf), drop = FALSE]
  P_treat <- q_cf_draws[, (n_cf + 1L):(2L * n_cf), drop = FALSE]

  dim(P_ctrl) <- c(M, tau, n_patients)
  dim(P_treat) <- c(M, tau, n_patients)

  log_1m_ctrl <- log(pmax(1 - P_ctrl, 1e-6))
  log_1m_treat <- log(pmax(1 - P_treat, 1e-6))

  cum_log_ctrl <- aperm(apply(log_1m_ctrl, c(1, 3), cumsum), c(2, 1, 3))
  cum_log_treat <- aperm(apply(log_1m_treat, c(1, 3), cumsum), c(2, 1, 3))

  S_ctrl_draws <- exp(cum_log_ctrl)
  S_treat_draws <- exp(cum_log_treat)

  marg_ctrl <- apply(S_ctrl_draws, c(1, 2), mean)
  marg_treat <- apply(S_treat_draws, c(1, 2), mean)

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
compute_bart_ipcw_and_cf <- function(data,
                                     clones,
                                     cf_base,
                                     g_grid_pred,
                                     covars,
                                     grace_period = 10,
                                     tau = 30,
                                     trim_quantile = c(0.001, 0.999)) {
  n_patients <- nrow(data)
  # g_grid_pred is ordered by patient i (1..n_patients), and within patient i by t = 1..grace_period
  g_mat <- matrix(bound_precision(g_grid_pred), nrow = grace_period, ncol = n_patients)
  colnames(g_mat) <- as.character(data$id)

  # Cumulative probability of initiating treatment within grace_period: 1 - prod_{s=1}^{grace_period} (1 - h_g(s | W_i))
  surv_no_tx_grace <- apply(1 - g_mat, 2, prod)
  p_treat_denom_vec <- bound_precision(1 - surv_no_tx_grace)
  names(p_treat_denom_vec) <- as.character(data$id)

  # Lookup table for (id, t) -> p_ctrl_denom
  g_lookup <- tibble::tibble(
    id = rep(data$id, each = grace_period),
    t = rep(seq_len(grace_period), times = n_patients),
    p_ctrl_denom_bart = as.vector(g_mat)
  )
  treat_lookup <- tibble::tibble(
    id = data$id,
    p_treat_denom_bart = as.numeric(p_treat_denom_vec)
  )

  # Marginal stabilization numerators (matching compute_ipcw_weights)
  obs_panel <- clones %>%
    dplyr::filter(clone_arm == 1L) %>%
    dplyr::mutate(
      tmax = pmin(dplyr::if_else(tx_ever == 1L, t_switch, time), time, grace_period),
      n_int = pmax(ceiling(tmax), 1L)
    ) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      t_sq = t^2,
      event_tx = dplyr::if_else(tx_ever == 1L & t == t_switch & time >= t_switch & t <= grace_period, 1L, 0L)
    )

  treat_day10_data <- clones %>%
    dplyr::filter(arm == "treat_tx", time > grace_period) %>%
    dplyr::mutate(took_tx = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))

  m_num_ctrl <- suppressWarnings(stats::glm(event_tx ~ t + t_sq, family = stats::binomial(), data = obs_panel))
  m_num_treat <- suppressWarnings(stats::glm(took_tx ~ 1, family = stats::binomial(), data = treat_day10_data))

  # Expand clones into longitudinal person-time panel
  panel <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      t0 = t - 1L,
      t1 = t,
      t_sq = t^2,
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > t0 & time_cc <= t1, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
    ) %>%
    dplyr::left_join(g_lookup, by = c("id", "t")) %>%
    dplyr::left_join(treat_lookup, by = "id") %>%
    dplyr::mutate(
      p_ctrl_denom = bound_precision(dplyr::coalesce(p_ctrl_denom_bart, 0.01)),
      p_ctrl_num = bound_precision(stats::predict(m_num_ctrl, newdata = ., type = "response")),
      p_treat_denom = bound_precision(p_treat_denom_bart),
      p_treat_num = bound_precision(stats::predict(m_num_treat, newdata = ., type = "response"))
    ) %>%
    dplyr::mutate(
      w_inc = dplyr::case_when(
        arm == "control_no_tx" & t <= grace_period ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom),
        arm == "control_no_tx" ~ 1,
        arm == "treat_tx" & t == grace_period + 1L ~ p_treat_num / p_treat_denom,
        TRUE ~ 1
      ),
      w_inc_unstab = dplyr::case_when(
        arm == "control_no_tx" & t <= grace_period ~ 1 / (1 - p_ctrl_denom),
        arm == "control_no_tx" ~ 1,
        arm == "treat_tx" & t == grace_period + 1L ~ 1 / p_treat_denom,
        TRUE ~ 1
      )
    ) %>%
    dplyr::arrange(id, clone_arm, t) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(
      w_ipcw_time = cumprod(w_inc),
      w_unstab_time = cumprod(w_inc_unstab)
    ) %>%
    dplyr::ungroup()

  q_lim <- stats::quantile(panel$w_ipcw_time, probs = trim_quantile, na.rm = TRUE)
  q_lim_unstab <- stats::quantile(panel$w_unstab_time, probs = trim_quantile, na.rm = TRUE)

  panel_w <- panel %>%
    dplyr::mutate(
      w_trim = pmin(pmax(w_ipcw_time, q_lim[1]), q_lim[2]),
      w_unstab_trim = pmin(pmax(w_unstab_time, q_lim_unstab[1]), q_lim_unstab[2]),
      post_grace = dplyr::if_else(t > grace_period, 1L, 0L)
    )

  # Counterfactual clever covariates cf_H on cf_base (t = 1..tau for each patient)
  cf_grid <- cf_base %>%
    dplyr::mutate(t_sq = t^2) %>%
    dplyr::left_join(g_lookup, by = c("id", "t")) %>%
    dplyr::left_join(treat_lookup, by = "id") %>%
    dplyr::mutate(
      cf_p_ctrl_denom = bound_precision(dplyr::coalesce(p_ctrl_denom_bart, 0.01)),
      cf_p_ctrl_num = bound_precision(stats::predict(m_num_ctrl, newdata = ., type = "response")),
      cf_p_treat_denom = bound_precision(p_treat_denom_bart),
      cf_p_treat_num = bound_precision(stats::predict(m_num_treat, newdata = ., type = "response")),
      w_inc_ctrl = dplyr::if_else(t <= grace_period, (1 - cf_p_ctrl_num) / (1 - cf_p_ctrl_denom), 1),
      w_inc_treat = dplyr::if_else(t == grace_period + 1L, cf_p_treat_num / cf_p_treat_denom, 1)
    ) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      cf_H_ctrl = pmin(pmax(cumprod(w_inc_ctrl), q_lim[1]), q_lim[2]),
      cf_H_treat = pmin(pmax(cumprod(w_inc_treat), q_lim[1]), q_lim[2])
    ) %>%
    dplyr::ungroup()

  list(
    panel_w = panel_w,
    cf_H_ctrl = cf_grid$cf_H_ctrl,
    cf_H_treat = cf_grid$cf_H_treat,
    q_lim = q_lim
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
#' @param ntree Number of trees in the BART ensemble (default 30 for fast bart-survival).
#' @param ndpost Number of retained posterior draws (default 100, i.e., 200 draws with thin = 2).
#' @param nskip Number of burn-in tune iterations (default 100).
#' @param thin Thinning factor (default 2).
#' @param n_boot Multiplier bootstrap replicates (default 5000).
#' @param seed Optional random seed.
#' @return A list with elements \code{gcomp} (output matching \code{fit_bart_gcomp}) and
#'   \code{tmle} (output matching \code{fit_bart_tmle}).
#' @export
fit_bart_pair <- function(data,
                          clones,
                          covars = c("w1", "w2"),
                          grace_period = 10,
                          tau = 30,
                          ntree = 30,
                          ndpost = 100,
                          nskip = 100,
                          thin = 2,
                          n_boot = 5000,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

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

  # Prediction grid for BART-g: every patient i across t = 1..grace_period
  g_grid <- data %>%
    dplyr::select(id, dplyr::all_of(covars)) %>%
    tidyr::uncount(weights = grace_period, .id = "t")
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

  # 3. Build counterfactual test grids (1..tau per patient)
  complier_times <- clones %>%
    dplyr::filter(arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period) %>%
    dplyr::pull(t_switch)
  if (length(complier_times) == 0) complier_times <- grace_period

  cf_base <- data %>%
    dplyr::select(id, dplyr::all_of(covars), tx_ever, t_switch) %>%
    dplyr::mutate(sampled_t = sample(complier_times, dplyr::n(), replace = TRUE)) %>%
    tidyr::uncount(weights = tau, .id = "t")

  cf_ctrl <- cf_base %>%
    dplyr::mutate(is_treated = 0L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  cf_treat <- cf_base %>%
    dplyr::mutate(
      t_sw = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, t_switch, sampled_t),
      is_treated = dplyr::if_else(t >= t_sw, 1L, 0L)
    ) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  X_q_test <- rbind(X_q_train, cf_ctrl, cf_treat)

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
    seed = if (is.null(seed)) 99L else as.integer(seed)
  )
  t_py_elapsed <- as.numeric(proc.time()["elapsed"] - t0_py)

  # 5. Evaluate G-COMP (BART) in R from posterior counterfactual draws
  t0_gcomp <- proc.time()["elapsed"]
  n_patients <- nrow(data)
  gcomp_res <- summarize_bart_gcomp_draws(
    q_cf_draws = py_out$q_cf_draws,
    n_patients = n_patients,
    tau = tau
  )
  gcomp_res$elapsed_sec <- 0.5 * t_py_elapsed + as.numeric(proc.time()["elapsed"] - t0_gcomp)

  # 6. Evaluate TMLE + BART in R (BART-g IPCW weights + BART-Q initial hazards + logistic targeting + EIC)
  t0_tmle <- proc.time()["elapsed"]
  ipcw_bart <- compute_bart_ipcw_and_cf(
    data = data,
    clones = clones,
    cf_base = cf_base,
    g_grid_pred = py_out$g_pred_mean,
    covars = covars,
    grace_period = grace_period,
    tau = tau
  )
  panel_w <- ipcw_bart$panel_w

  n_cf <- nrow(cf_base)
  Q_init_train <- py_out$q_train_mean
  Q_init_ctrl_cf <- colMeans(py_out$q_cf_draws[, seq_len(n_cf), drop = FALSE])
  Q_init_treat_cf <- colMeans(py_out$q_cf_draws[, (n_cf + 1L):(2L * n_cf), drop = FALSE])

  # TMLE logistic fluctuation step in R
  targ_res <- target_hazard(
    panel_w = panel_w,
    pred_Q_init = Q_init_train,
    grace_period = grace_period
  )
  panel_fluc <- targ_res$panel_w
  eps_ctrl <- targ_res$eps_ctrl
  eps_treat_pre <- targ_res$eps_treat_pre
  eps_treat_post <- targ_res$eps_treat_post

  # Update counterfactual hazards using BART-g clever covariates cf_H
  Q_star_ctrl <- stats::plogis(
    stats::qlogis(bound_precision(Q_init_ctrl_cf)) + eps_ctrl * ipcw_bart$cf_H_ctrl
  )
  Q_star_treat <- stats::plogis(
    stats::qlogis(bound_precision(Q_init_treat_cf)) +
      ifelse(cf_base$t <= grace_period, eps_treat_pre, eps_treat_post) * ipcw_bart$cf_H_treat
  )

  dim(Q_star_ctrl) <- c(tau, n_patients)
  dim(Q_star_treat) <- c(tau, n_patients)

  S_mat_ctrl <- apply(1 - Q_star_ctrl, 2, cumprod)
  S_mat_treat <- apply(1 - Q_star_treat, 2, cumprod)

  mean_s0 <- rowMeans(S_mat_ctrl)
  mean_s1 <- rowMeans(S_mat_treat)
  mean_rd <- mean_s1 - mean_s0

  cf_ctrl_df <- cf_base %>%
    dplyr::mutate(arm = "control_no_tx", indiv_survival = as.vector(S_mat_ctrl))
  cf_treat_df <- cf_base %>%
    dplyr::mutate(arm = "treat_tx", indiv_survival = as.vector(S_mat_treat))
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
#' @param ntree Number of trees in the BART ensemble (default 30).
#' @param ndpost Number of retained posterior MCMC draws (default 100).
#' @param nskip Number of MCMC burn-in iterations (default 100).
#' @param thin Thinning factor (default 2).
#' @param mc.cores Unused when Python \code{bart-survival} backend is active.
#' @param seed Optional random seed.
#' @export
fit_bart_gcomp <- function(data,
                           clones,
                           covars = c("w1", "w2"),
                           grace_period = 10,
                           tau = 30,
                           ntree = 30,
                           ndpost = 100,
                           nskip = 100,
                           thin = 2,
                           mc.cores = 1,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  panel_bart <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > (t - 1L) & time_cc <= t, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
    )

  X_train <- panel_bart %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()
  Y_train <- panel_bart$event_cc_day

  complier_times <- clones %>%
    dplyr::filter(arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period) %>%
    dplyr::pull(t_switch)
  if (length(complier_times) == 0) complier_times <- grace_period

  cf_base <- data %>%
    dplyr::select(id, dplyr::all_of(covars), tx_ever, t_switch) %>%
    dplyr::mutate(sampled_t = sample(complier_times, dplyr::n(), replace = TRUE)) %>%
    tidyr::uncount(weights = tau, .id = "t")

  cf_ctrl <- cf_base %>%
    dplyr::mutate(is_treated = 0L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  cf_treat <- cf_base %>%
    dplyr::mutate(
      t_sw = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, t_switch, sampled_t),
      is_treated = dplyr::if_else(t >= t_sw, 1L, 0L)
    ) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()

  X_test <- rbind(cf_ctrl, cf_treat)

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
    seed = if (is.null(seed)) 99L else as.integer(seed)
  )

  summarize_bart_gcomp_draws(
    q_cf_draws = py_out$q_cf_draws,
    n_patients = nrow(data),
    tau = tau
  )
}

#' Targeted Minimum Loss-Based Estimation via Survival-BART for Both g and Q (TMLE + BART)
#'
#' @param data Raw standardized data frame.
#' @param clones Cloned data frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param ntree Number of trees in the BART ensemble (default 30).
#' @param ndpost Number of retained posterior MCMC draws (default 100).
#' @param nskip Number of MCMC burn-in iterations (default 100).
#' @param thin Thinning factor (default 2).
#' @param mc.cores Unused when Python \code{bart-survival} backend is active.
#' @param n_boot Multiplier bootstrap replicates (default 5000).
#' @param seed Optional random seed.
#' @export
fit_bart_tmle <- function(data,
                          clones,
                          covars = c("w1", "w2"),
                          grace_period = 10,
                          tau = 30,
                          ntree = 30,
                          ndpost = 100,
                          nskip = 100,
                          thin = 2,
                          mc.cores = 1,
                          n_boot = 5000,
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
    n_boot = n_boot,
    seed = seed
  )
  pair_res$tmle
}
