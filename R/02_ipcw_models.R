#' Compute Grace-Period Initiation Survival Grid
#'
#' @description Predicts daily treatment initiation hazards across the grace window
#'   for each patient and computes cumulative untreated survival probabilities
#'   \eqn{S_{0,A}(t \mid W)}, daily unconditional initiation probabilities
#'   \eqn{p_0(k \mid W)}, and end-of-grace probabilities \eqn{S_{0,A}(m \mid W)}.
#'
#' @keywords internal
#' @noRd
compute_grace_survival_grid <- function(data_patients,
                                        covars,
                                        grace_period = 10,
                                        dt = 1,
                                        m_denom_ctrl,
                                        m_num_ctrl,
                                        use_sl = FALSE,
                                        g_mat_override = NULL) {
  n_grace_int <- max(1L, as.integer(ceiling(grace_period / dt)))
  base_pat <- data_patients %>%
    dplyr::select(id, dplyr::all_of(covars)) %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::arrange(id)

  ids <- base_pat$id
  n_pat <- length(ids)

  g_grid <- tidyr::expand_grid(
    id = ids,
    t = seq_len(n_grace_int)
  ) %>%
    dplyr::left_join(base_pat, by = "id") %>%
    dplyr::mutate(t_sq = t^2)

  if (!is.null(g_mat_override)) {
    # g_mat_override is (n_grace_int x n_pat) ordered by columns matching ids
    if (!is.null(colnames(g_mat_override))) {
      g_mat_override <- g_mat_override[, as.character(ids), drop = FALSE]
    }
    lam_denom <- as.vector(g_mat_override)
  } else if (use_sl) {
    X_pred <- g_grid %>% dplyr::select(dplyr::all_of(covars), t, t_sq) %>% as.data.frame()
    lam_denom <- stats::predict(m_denom_ctrl, newdata = X_pred, onlySL = TRUE)$pred[, 1]
  } else {
    lam_denom <- stats::predict(m_denom_ctrl, newdata = g_grid, type = "response")
  }
  lam_num <- stats::predict(m_num_ctrl, newdata = g_grid, type = "response")

  g_grid <- g_grid %>%
    dplyr::mutate(
      lam_denom = bound_precision(lam_denom),
      lam_num   = bound_precision(lam_num)
    ) %>%
    dplyr::arrange(id, t) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      S_A_denom = cumprod(1 - lam_denom),
      S_A_num   = cumprod(1 - lam_num),
      p0_day    = dplyr::lag(S_A_denom, default = 1) * lam_denom
    ) %>%
    dplyr::ungroup()

  s_m_lookup <- g_grid %>%
    dplyr::filter(t == n_grace_int) %>%
    dplyr::select(id, S_m_denom = S_A_denom, S_m_num = S_A_num)

  # Build (n_pat x n_grace_int) matrices sorted by (id, t) for fast vectorized integration
  S_denom_mat <- matrix(g_grid$S_A_denom, nrow = n_pat, ncol = n_grace_int, byrow = TRUE)
  S_num_mat   <- matrix(g_grid$S_A_num,   nrow = n_pat, ncol = n_grace_int, byrow = TRUE)
  p0_mat      <- matrix(g_grid$p0_day,    nrow = n_pat, ncol = n_grace_int, byrow = TRUE)

  list(
    g_grid = g_grid %>% dplyr::select(id, t, S_A_denom, S_A_num, p0_day),
    s_m_lookup = s_m_lookup,
    ids = ids,
    S_denom_mat = S_denom_mat,
    S_num_mat = S_num_mat,
    p0_mat = p0_mat,
    n_grace_int = n_grace_int
  )
}

#' Compute Empirical Grace-Period Complier Initiation Probabilities
#'
#' @keywords internal
#' @noRd
compute_empirical_complier_probs <- function(clones, grace_period = 10, n_grace_int = 10L) {
  comp_days <- clones %>%
    dplyr::filter(
      arm == "treat_tx",
      tx_ever == 1L,
      !is.na(t_switch),
      t_switch >= 1,
      t_switch <= grace_period,
      time >= t_switch
    ) %>%
    dplyr::pull(t_switch)

  k_seq <- seq_len(n_grace_int)
  if (length(comp_days) == 0L) {
    return(rep(1 / n_grace_int, n_grace_int))
  }
  comp_int <- pmax(1L, pmin(n_grace_int, as.integer(round(comp_days))))
  tab <- tabulate(comp_int, nbins = n_grace_int)
  tab / sum(tab)
}

#' Compute Treatment-Arm Counterfactual Switch-Day Strategy Weights
#'
#' @keywords internal
#' @noRd
compute_tx_strategy_weights <- function(grace_surv,
                                        tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                                        p_star = NULL) {
  tx_strategy <- match.arg(tx_strategy)
  m <- grace_surv$n_grace_int
  n_pat <- length(grace_surv$ids)

  if (tx_strategy == "conditional_complier") {
    p0_mat <- grace_surv$p0_mat
    rs <- pmax(rowSums(p0_mat), 1e-8)
    w_mat <- p0_mat / rs
    k_vals <- seq_len(m)
  } else if (tx_strategy == "empirical_complier") {
    if (is.null(p_star) || length(p_star) != m) {
      stop("p_star of length n_grace_int is required for empirical_complier mode.")
    }
    p_star <- pmax(p_star, 0)
    p_star <- p_star / sum(p_star)
    w_mat <- matrix(p_star, nrow = n_pat, ncol = m, byrow = TRUE)
    k_vals <- seq_len(m)
  } else {
    # Hernan natural_grace: natural initiation probability p0(k | W) for k = 1..m,
    # plus non-initiator mass S_{0,A}(m | W) initiating immediately at grace_period + 1
    p0_mat <- grace_surv$p0_mat
    s_m <- pmax(1 - rowSums(p0_mat), 0)
    w_mat <- cbind(p0_mat, s_m)
    k_vals <- seq_len(m + 1L)
  }

  list(w_mat = w_mat, k_vals = k_vals)
}

#' Compute Counterfactual Clever Covariate Matrices for Analytical Splicing
#'
#' @keywords internal
#' @noRd
compute_cf_weight_matrices <- function(grace_surv,
                                       data_patients,
                                       covars,
                                       grace_period = 10,
                                       tau = 30,
                                       tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                                       k_vals,
                                       q_lim,
                                       m_denom_treat = NULL,
                                       m_num_treat = NULL,
                                       use_sl = FALSE,
                                       p_star = NULL,
                                       p_treat_denom_vec = NULL,
                                       p_treat_num_vec = NULL) {
  tx_strategy <- match.arg(tx_strategy)
  m <- grace_surv$n_grace_int
  n_pat <- length(grace_surv$ids)
  t_m <- min(m, tau)

  S_denom_full <- matrix(grace_surv$S_denom_mat[, m], nrow = n_pat, ncol = tau)
  S_num_full   <- matrix(grace_surv$S_num_mat[, m],   nrow = n_pat, ncol = tau)
  S_denom_full[, seq_len(t_m)] <- grace_surv$S_denom_mat[, seq_len(t_m), drop = FALSE]
  S_num_full[, seq_len(t_m)]   <- grace_surv$S_num_mat[, seq_len(t_m), drop = FALSE]

  # Control counterfactual clever covariate matrix (n_pat x tau)
  cf_H_ctrl_raw <- S_num_full / bound_precision(S_denom_full)
  cf_H_ctrl <- matrix(
    apply_arm_weight_trimming(as.vector(cf_H_ctrl_raw), rep("control_no_tx", n_pat * tau), q_lim),
    nrow = n_pat, ncol = tau
  )

  # Treated counterfactual clever covariate matrix per switch day k (list of n_pat x tau)
  cf_H_treat_list <- vector("list", length(k_vals))
  s_m_vec <- grace_surv$S_denom_mat[, m]

  if (tx_strategy == "conditional_complier") {
    post_col <- 1 / bound_precision(1 - s_m_vec)
    pre_mat <- sweep(pmax(S_denom_full - s_m_vec, 0), 1, bound_precision(1 - s_m_vec), "/") /
      bound_precision(S_denom_full)

    for (idx in seq_along(k_vals)) {
      k <- k_vals[idx]
      raw_mat <- matrix(post_col, nrow = n_pat, ncol = tau)
      pre_t <- seq_len(min(k - 1L, t_m))
      if (length(pre_t) > 0L) {
        raw_mat[, pre_t] <- pre_mat[, pre_t, drop = FALSE]
      }
      cf_H_treat_list[[idx]] <- matrix(
        apply_arm_weight_trimming(as.vector(raw_mat), rep("treat_tx", n_pat * tau), q_lim),
        nrow = n_pat, ncol = tau
      )
    }
  } else if (tx_strategy == "empirical_complier") {
    # Marginal probability of remaining untreated past day t: S*(t) = sum_{j > t} p*(j)
    S_star_vec <- c(rev(cumsum(rev(p_star[-1L]))), 0)
    if (length(p_star) == 1L) S_star_vec <- 0
    pre_mat <- sweep(1 / bound_precision(S_denom_full[, seq_len(t_m), drop = FALSE]), 2, S_star_vec[seq_len(t_m)], "*")

    for (idx in seq_along(k_vals)) {
      k <- k_vals[idx]
      post_col <- if (p_star[k] > 0) p_star[k] / bound_precision(grace_surv$p0_mat[, k]) else rep(0, n_pat)
      raw_mat <- matrix(post_col, nrow = n_pat, ncol = tau)
      pre_t <- seq_len(min(k - 1L, t_m))
      if (length(pre_t) > 0L) {
        raw_mat[, pre_t] <- pre_mat[, pre_t, drop = FALSE]
      }
      cf_H_treat_list[[idx]] <- matrix(
        apply_arm_weight_trimming(as.vector(raw_mat), rep("treat_tx", n_pat * tau), q_lim),
        nrow = n_pat, ncol = tau
      )
    }
  } else {
    # natural_grace: w = 1 for t <= m, p_treat_num / p_treat_denom for t > m
    base_ord <- data_patients %>%
      dplyr::select(id, dplyr::all_of(covars)) %>%
      dplyr::distinct(id, .keep_all = TRUE) %>%
      dplyr::arrange(id)

    if (is.null(p_treat_denom_vec)) {
      if (use_sl) {
        X_t <- base_ord %>% dplyr::select(dplyr::all_of(covars)) %>% as.data.frame()
        p_treat_denom_vec <- stats::predict(m_denom_treat, newdata = X_t, onlySL = TRUE)$pred[, 1]
      } else {
        p_treat_denom_vec <- stats::predict(m_denom_treat, newdata = base_ord, type = "response")
      }
    }
    if (is.null(p_treat_num_vec)) {
      p_treat_num_vec <- stats::predict(m_num_treat, newdata = base_ord, type = "response")
    }
    p_treat_denom_vec <- bound_precision(p_treat_denom_vec)
    p_treat_num_vec   <- bound_precision(p_treat_num_vec)
    post_col <- p_treat_num_vec / p_treat_denom_vec

    raw_mat <- matrix(1, nrow = n_pat, ncol = tau)
    if (tau > m) {
      raw_mat[, (m + 1L):tau] <- matrix(post_col, nrow = n_pat, ncol = tau - m)
    }
    trimmed_ng <- matrix(
      apply_arm_weight_trimming(as.vector(raw_mat), rep("treat_tx", n_pat * tau), q_lim),
      nrow = n_pat, ncol = tau
    )
    for (idx in seq_along(k_vals)) {
      cf_H_treat_list[[idx]] <- trimmed_ng
    }
  }

  list(cf_H_ctrl = cf_H_ctrl, cf_H_treat_list = cf_H_treat_list)
}

#' Compute Arm-Stratified IPCW / Stochastic Grace-Period Weights
#'
#' @keywords internal
#' @noRd
compute_arm_ipcw_weights <- function(df_panel,
                                     grace_surv,
                                     covars,
                                     grace_period = 10,
                                     m_denom_treat = NULL,
                                     m_num_treat = NULL,
                                     use_sl = FALSE,
                                     tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                                     p_star = NULL,
                                     p_treat_denom_vec = NULL,
                                     p_treat_num_vec = NULL,
                                     is_cf = FALSE) {
  tx_strategy <- match.arg(tx_strategy)
  n_grace_int <- grace_surv$n_grace_int

  df_joined <- df_panel %>%
    dplyr::left_join(grace_surv$g_grid %>% dplyr::select(id, t, S_A_denom, S_A_num), by = c("id", "t")) %>%
    dplyr::left_join(grace_surv$s_m_lookup, by = "id") %>%
    dplyr::mutate(
      S_A_denom = dplyr::if_else(t > n_grace_int, S_m_denom, S_A_denom),
      S_A_num   = dplyr::if_else(t > n_grace_int, S_m_num, S_A_num)
    )

  sw_col <- if (is_cf) df_joined$cf_t_switch else df_joined$t_switch
  in_grace_switch <- if (is_cf) {
    sw_col <= grace_period
  } else {
    df_joined$tx_ever == 1L & !is.na(sw_col) & sw_col <= grace_period
  }

  if (tx_strategy == "conditional_complier") {
    # Pre-switch uses (S_A_denom(t) - S_m_denom) / (S_A_denom(t) * (1 - S_m_denom));
    # post-switch uses 1 / (1 - S_m_denom). E[H | W] == 1 identically (self-stabilized).
    df_joined <- df_joined %>%
      dplyr::mutate(
        in_grace_sw = in_grace_switch,
        w_stab = dplyr::case_when(
          arm == "control_no_tx" ~ S_A_num / bound_precision(S_A_denom),
          arm == "treat_tx" & is_treated == 0 & t <= n_grace_int ~
            pmax(S_A_denom - S_m_denom, 0) / (bound_precision(S_A_denom) * bound_precision(1 - S_m_denom)),
          arm == "treat_tx" & is_treated == 1 & in_grace_sw ~
            1 / bound_precision(1 - S_m_denom),
          TRUE ~ 0
        ),
        w_unstab = dplyr::case_when(
          arm == "control_no_tx" ~ 1 / bound_precision(S_A_denom),
          arm == "treat_tx" ~ w_stab,
          TRUE ~ 0
        )
      ) %>%
      dplyr::select(-in_grace_sw)
  } else if (tx_strategy == "empirical_complier") {
    # Marginal empirical complier policy:
    # Pre-switch: sum_{j > t} p*(j) / S_{0,A}(t | W); Post-switch: p*(T_A) / p_0(T_A | W).
    if (is.null(p_star) || length(p_star) != n_grace_int) {
      stop("p_star is required for tx_strategy = 'empirical_complier'.")
    }
    S_star_vec <- c(rev(cumsum(rev(p_star[-1L]))), 0)
    if (length(p_star) == 1L) S_star_vec <- 0
    s_star_tbl <- tibble::tibble(t = seq_len(n_grace_int), S_star_t = S_star_vec)
    p_star_tbl <- tibble::tibble(sw_day = seq_len(n_grace_int), p_star_sw = p_star)

    sw_day_int <- dplyr::if_else(
      in_grace_switch,
      pmax(1L, pmin(n_grace_int, as.integer(round(dplyr::if_else(in_grace_switch, sw_col, 1))))),
      NA_integer_
    )

    p0_sw_tbl <- grace_surv$g_grid %>%
      dplyr::select(id, sw_day = t, p0_sw = p0_day)

    df_joined <- df_joined %>%
      dplyr::mutate(in_grace_sw = in_grace_switch, sw_day = sw_day_int) %>%
      dplyr::left_join(s_star_tbl, by = "t") %>%
      dplyr::left_join(p_star_tbl, by = "sw_day") %>%
      dplyr::left_join(p0_sw_tbl, by = c("id", "sw_day")) %>%
      dplyr::mutate(
        S_star_t = tidyr::replace_na(S_star_t, 0),
        w_stab = dplyr::case_when(
          arm == "control_no_tx" ~ S_A_num / bound_precision(S_A_denom),
          arm == "treat_tx" & is_treated == 0 & t <= n_grace_int ~
            S_star_t / bound_precision(S_A_denom),
          arm == "treat_tx" & is_treated == 1 & in_grace_sw ~
            tidyr::replace_na(p_star_sw, 0) / bound_precision(tidyr::replace_na(p0_sw, 1e-4)),
          TRUE ~ 0
        ),
        w_unstab = dplyr::case_when(
          arm == "control_no_tx" ~ 1 / bound_precision(S_A_denom),
          arm == "treat_tx" ~ w_stab,
          TRUE ~ 0
        )
      ) %>%
      dplyr::select(-in_grace_sw, -sw_day, -S_star_t, -p_star_sw, -p0_sw)
  } else {
    # Textbook Hernan natural_grace: w = 1 for t <= m, IPCW P(T_A <= m)/P(T_A <= m | W) for t > m
    if (!is.null(p_treat_denom_vec) && !is.null(p_treat_num_vec)) {
      t_lkp <- tibble::tibble(
        id = grace_surv$ids,
        p_treat_denom = bound_precision(p_treat_denom_vec),
        p_treat_num   = bound_precision(p_treat_num_vec)
      )
      df_joined <- df_joined %>% dplyr::left_join(t_lkp, by = "id")
    } else {
      if (use_sl) {
        X_treat <- df_joined %>% dplyr::select(dplyr::all_of(covars)) %>% as.data.frame()
        pred_treat_denom <- stats::predict(m_denom_treat, newdata = X_treat, onlySL = TRUE)$pred[, 1]
      } else {
        pred_treat_denom <- stats::predict(m_denom_treat, newdata = df_joined, type = "response")
      }
      df_joined$p_treat_denom <- bound_precision(pred_treat_denom)
      df_joined$p_treat_num <- bound_precision(stats::predict(m_num_treat, newdata = df_joined, type = "response"))
    }

    df_joined <- df_joined %>%
      dplyr::mutate(
        w_stab = dplyr::case_when(
          arm == "control_no_tx" ~ S_A_num / bound_precision(S_A_denom),
          arm == "treat_tx" & t > n_grace_int ~ p_treat_num / p_treat_denom,
          TRUE ~ 1
        ),
        w_unstab = dplyr::case_when(
          arm == "control_no_tx" ~ 1 / bound_precision(S_A_denom),
          arm == "treat_tx" & t > n_grace_int ~ 1 / p_treat_denom,
          TRUE ~ 1
        )
      )
  }

  df_joined
}

#' Compute Arm-Stratified Trimming Quantile Bounds
#'
#' @keywords internal
#' @noRd
compute_arm_trim_bounds <- function(w_vec, arm_vec, probs = c(0.001, 0.999), max_weight = Inf) {
  get_bounds <- function(vals) {
    pos <- vals[vals > 0 & !is.na(vals)]
    if (length(pos) == 0L) return(c(0, 1))
    b <- stats::quantile(pos, probs = probs, na.rm = TRUE)
    if (is.finite(max_weight)) b[2] <- min(b[2], max_weight)
    b
  }
  list(
    control_no_tx = get_bounds(w_vec[arm_vec == "control_no_tx"]),
    treat_tx      = get_bounds(w_vec[arm_vec == "treat_tx"]),
    pooled        = get_bounds(w_vec)
  )
}

#' Apply Arm-Stratified Weight Trimming (Preserving Structural Zeros)
#'
#' @keywords internal
#' @noRd
apply_arm_weight_trimming <- function(w_vec, arm_vec, q_lim) {
  if (is.list(q_lim) && all(c("control_no_tx", "treat_tx") %in% names(q_lim))) {
    b_ctrl <- q_lim$control_no_tx
    b_trt  <- q_lim$treat_tx
    lo <- dplyr::if_else(arm_vec == "control_no_tx", b_ctrl[1], b_trt[1])
    hi <- dplyr::if_else(arm_vec == "control_no_tx", b_ctrl[2], b_trt[2])
  } else {
    lo <- q_lim[1]
    hi <- q_lim[2]
  }
  dplyr::if_else(w_vec > 0, pmin(pmax(w_vec, lo), hi), 0)
}

#' Fit Compliance Models and Compute Inverse Probability of Censoring Weights (IPCW)
#'
#' @description Fits pooled logistic or SuperLearner compliance models on observational person-time,
#'   expands clones into a discrete person-time panel, and calculates stabilized, arm-stratified
#'   trimmed cumulative inverse probability of censoring weights (and aligned stochastic
#'   grace-period weights for the treated arm when \code{align_grace_weights = TRUE}).
#'
#' @param clones Cloned data.frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param dt Discrete time interval step size (default 1).
#' @param p_form Optional formula string for the compliance/propensity models. If \code{NULL},
#'   defaults to additive linear terms of \code{covars}.
#' @param use_sl Logical indicating whether to use \code{SuperLearner} for denominator compliance models.
#' @param sl_lib Character vector of SuperLearner algorithm wrappers.
#' @param trim_quantile Numeric vector of length 2 specifying lower and upper percentiles for weight trimming.
#' @param align_grace_weights Deprecated alias for \code{tx_strategy}. If explicitly set to
#'   \code{FALSE} while \code{tx_strategy} is omitted, selects \code{"natural_grace"}.
#' @param tx_strategy Treatment-arm initiation policy during the grace period:
#'   \code{"conditional_complier"} (default), \code{"empirical_complier"}, or \code{"natural_grace"}.
#' @param max_weight Numeric upper cap on stabilized weights to prevent spline tail separation
#'   (default 10; unstabilized weights are capped at \code{2 * max_weight}).
#'
#' @return A list containing:
#'   \describe{
#'     \item{panel_w}{Longitudinal person-time panel of cloned observations with trimmed weights \code{w_trim}.}
#'     \item{fit_ccw}{Survival fit object from \code{survival::survfit} using IPCW weights.}
#'     \item{m_denom_ctrl}{Fitted compliance model for the control arm.}
#'     \item{m_denom_treat}{Fitted compliance model for the treated arm.}
#'     \item{m_num_ctrl}{Fitted numerator model for control stabilization.}
#'     \item{m_num_treat}{Fitted numerator model for treated stabilization.}
#'     \item{grace_surv}{Precomputed grace-period survival grid from \code{m_denom_ctrl} and \code{m_num_ctrl}.}
#'     \item{tx_strategy}{Resolved treatment-arm policy.}
#'     \item{p_star}{Empirical complier initiation probability vector across days \code{1..m}.}
#'     \item{align_grace_weights}{Logical flag (\code{TRUE} for complier modes, \code{FALSE} for \code{"natural_grace"}).}
#'     \item{q_lim}{Arm-stratified weight trimming limits applied.}
#'   }
#' @export
#' @importFrom dplyr select mutate filter arrange group_by ungroup if_else case_when all_of distinct left_join
#' @importFrom tidyr uncount expand_grid replace_na
#' @importFrom stats as.formula binomial glm predict quantile
#' @importFrom survival survfit Surv
#' @importFrom splines ns
compute_ipcw_weights <- function(clones,
                                 covars = c("x1"),
                                 grace_period = 10,
                                 dt = 1,
                                 p_form = NULL,
                                 use_sl = FALSE,
                                 sl_lib = c("SL.mean", "SL.glm", "SL.hal.time", "SL.earth.int", "SL.ranger"),
                                 trim_quantile = c(0.001, 0.999),
                                 align_grace_weights = TRUE,
                                 tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                                 max_weight = 10) {
  if (missing(tx_strategy) && !isTRUE(align_grace_weights)) {
    tx_strategy <- "natural_grace"
  } else {
    tx_strategy <- match.arg(tx_strategy)
  }
  align_grace_weights <- (tx_strategy != "natural_grace")

  if (is.null(p_form)) {
    p_form <- paste(covars, collapse = " + ")
  }

  if (use_sl) {
    if (!require("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required when use_sl = TRUE.")
    }
    register_sl_earth_int()
  }

  # 1. Observational panel for control compliance hazard (restricted to grace period while alive)
  obs_panel <- clones %>%
    dplyr::filter(clone_arm == 1L) %>%
    dplyr::mutate(
      tmax = pmin(dplyr::if_else(tx_ever == 1L, t_switch, time), time, grace_period),
      n_int = pmax(ceiling(tmax / dt), 1L)
    ) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      t_sq = t^2,
      event_tx = dplyr::if_else(tx_ever == 1L & t == t_switch & time >= t_switch & t <= grace_period, 1L, 0L)
    )

  # 2. Treated cohort at grace period cutoff
  treat_day10_data <- clones %>%
    dplyr::filter(arm == "treat_tx", time > grace_period) %>%
    dplyr::mutate(took_tx = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))

  # Basis for discrete time inside grace period (includes t x W interaction for dynamic depletion)
  df_t <- min(3L, max(1L, length(unique(obs_panel$t)) - 1L))
  t_basis <- if (df_t > 1L) paste0("splines::ns(t, df = ", df_t, ")") else "t"

  # 3. Fit Denominator Compliance Models
  if (use_sl) {
    X_ctrl <- obs_panel %>% dplyr::select(dplyr::all_of(covars), t, t_sq) %>% as.data.frame()
    m_denom_ctrl <- SuperLearner::SuperLearner(
      Y = obs_panel$event_tx,
      X = X_ctrl,
      family = stats::binomial(),
      SL.library = sl_lib,
      cvControl = list(V = 2)
    )

    if (tx_strategy == "natural_grace") {
      X_treat <- treat_day10_data %>% dplyr::select(dplyr::all_of(covars)) %>% as.data.frame()
      m_denom_treat <- SuperLearner::SuperLearner(
        Y = treat_day10_data$took_tx,
        X = X_treat,
        family = stats::binomial(),
        SL.library = sl_lib,
        cvControl = list(V = 2)
      )
    } else {
      m_denom_treat <- NULL
    }
  } else {
    ctrl_formula <- stats::as.formula(paste0("event_tx ~ ", t_basis, " * (", p_form, ")"))
    treat_formula <- stats::as.formula(paste("took_tx ~", p_form))
    m_denom_ctrl <- suppressWarnings(stats::glm(ctrl_formula, family = stats::binomial(), data = obs_panel))
    m_denom_treat <- suppressWarnings(stats::glm(treat_formula, family = stats::binomial(), data = treat_day10_data))
  }

  # 4. Numerator Models for Stabilization
  m_num_ctrl  <- suppressWarnings(stats::glm(stats::as.formula(paste("event_tx ~", t_basis)), family = stats::binomial(), data = obs_panel))
  m_num_treat <- suppressWarnings(stats::glm(took_tx ~ 1, family = stats::binomial(), data = treat_day10_data))

  # 5. Compute patient-level grace-period survival probabilities S_{0,A}(t | W)
  base_patients <- clones %>% dplyr::filter(clone_arm == 1L)
  grace_surv <- compute_grace_survival_grid(
    data_patients = base_patients,
    covars = covars,
    grace_period = grace_period,
    dt = dt,
    m_denom_ctrl = m_denom_ctrl,
    m_num_ctrl = m_num_ctrl,
    use_sl = use_sl
  )

  p_star <- compute_empirical_complier_probs(
    clones = clones,
    grace_period = grace_period,
    n_grace_int = grace_surv$n_grace_int
  )

  # 6. Expand Clones into Longitudinal Person-Time Panel
  panel <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc / dt), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      t0 = (t - 1L) * dt,
      t1 = t * dt,
      t_sq = t^2,
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > t0 & time_cc <= t1, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1, 0)
    )

  # 7. Calculate Arm-Specific Weights via Shared Helper
  panel <- compute_arm_ipcw_weights(
    df_panel = panel,
    grace_surv = grace_surv,
    covars = covars,
    grace_period = grace_period,
    m_denom_treat = m_denom_treat,
    m_num_treat = m_num_treat,
    use_sl = use_sl,
    tx_strategy = tx_strategy,
    p_star = p_star
  ) %>%
    dplyr::rename(
      w_ipcw_time = w_stab,
      w_unstab_time = w_unstab
    ) %>%
    dplyr::arrange(clone_id, t)

  # 8. Arm-Stratified Weight Trimming (preserving structural zeros and enforcing max_weight ceiling)
  q_lim <- compute_arm_trim_bounds(panel$w_ipcw_time, panel$arm, probs = trim_quantile, max_weight = max_weight)
  q_lim_unstab <- compute_arm_trim_bounds(panel$w_unstab_time, panel$arm, probs = trim_quantile, max_weight = max_weight * 2)

  panel_w <- panel %>%
    dplyr::mutate(
      w_trim = apply_arm_weight_trimming(w_ipcw_time, arm, q_lim),
      w_unstab_trim = apply_arm_weight_trimming(w_unstab_time, arm, q_lim_unstab),
      post_grace = dplyr::if_else(t > grace_period, 1, 0)
    )

  fit_ccw <- survival::survfit(
    survival::Surv(t0, t1, event_cc_day) ~ arm,
    data = panel_w[panel_w$w_trim > 0, , drop = FALSE],
    weights = w_trim,
    id = clone_id,
    robust = TRUE
  )

  return(list(
    panel_w = panel_w,
    fit_ccw = fit_ccw,
    m_denom_ctrl = m_denom_ctrl,
    m_denom_treat = m_denom_treat,
    m_num_ctrl = m_num_ctrl,
    m_num_treat = m_num_treat,
    grace_surv = grace_surv,
    tx_strategy = tx_strategy,
    p_star = p_star,
    align_grace_weights = align_grace_weights,
    q_lim = q_lim
  ))
}

