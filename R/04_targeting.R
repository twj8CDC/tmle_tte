#' Decoupled Clever Covariates and TMLE Fluctuation Step
#'
#' @description Implements the TMLE hazard targeting step via either:
#'   \itemize{
#'     \item \code{"pooled_glm"} (default): Decoupled 3-parameter pooled logistic fluctuation
#'       separating control, pre-grace, and post-grace windows.
#'     \item \code{"one_step"}: Universal One-Step Survival TMLE (Cai & van der Laan 2020 / MOSS)
#'       iteratively updating hazard surfaces along the joint pathwise gradient across the entire
#'       follow-up horizon \eqn{\tau}, eliminating mid-follow-up residual cancellation across time.
#'   }
#'
#' @param panel_w Longitudinal person-time panel from \code{\link{compute_ipcw_weights}}.
#' @param pred_Q_init Initial hazard predictions on \code{panel_w} from \code{\link{fit_q_models}}.
#' @param grace_period Numeric grace period duration (default 10).
#' @param targeting_method Targeting algorithm: \code{"pooled_glm"} (default) or \code{"one_step"}.
#' @param base_ord Baseline patient cohort ordered by id (required for \code{"one_step"}).
#' @param strat Strategy weights list from \code{\link{compute_tx_strategy_weights}} (required for \code{"one_step"}).
#' @param cf_w Counterfactual weight matrices from \code{\link{compute_cf_weight_matrices}} (required for \code{"one_step"}).
#' @param logit_ctrl_init Initial logit hazard matrix for control arm (N by tau; required for \code{"one_step"}).
#' @param logit_tr0_init Initial logit hazard matrix for treated arm untreated state (N by tau; required for \code{"one_step"}).
#' @param logit_tr1_init Initial logit hazard matrix for treated arm treated state (N by tau; required for \code{"one_step"}).
#' @param tau Follow-up time horizon (default 30).
#' @param max_iter Maximum iterations for one-step targeting (default 15L).
#' @param tol Convergence tolerance on maximum absolute EIC residual \eqn{\max_{t_0} |c_{\mathrm{RD}}(t_0)|} (default 0.005).
#' @param damping Damping factor for one-step GLM updates to prevent oscillation (default 0.5).
#'
#' @return A list containing:
#'   \describe{
#'     \item{panel_w}{Panel updated with initial/targeted hazard logits and clever covariates.}
#'     \item{m_fluctuate}{Fitted fluctuation regression object (or \code{NULL} for \code{"one_step"}).}
#'     \item{eps_ctrl}{Estimated fluctuation coefficient for control arm (0 for \code{"one_step"}).}
#'     \item{eps_treat_pre}{Estimated fluctuation coefficient for treatment arm pre-grace (0 for \code{"one_step"}).}
#'     \item{eps_treat_post}{Estimated fluctuation coefficient for treatment arm post-grace (0 for \code{"one_step"}).}
#'     \item{eps}{Named vector of fluctuation coefficients.}
#'     \item{Q_star_ctrl_mat}{Pre-targeted control hazard matrix (for \code{"one_step"}).}
#'     \item{Q_star_tr0_mat}{Pre-targeted treated pre-switch hazard matrix (for \code{"one_step"}).}
#'     \item{Q_star_tr1_mat}{Pre-targeted treated post-switch hazard matrix (for \code{"one_step"}).}
#'     \item{c_rd}{Final EIC residual vector across \eqn{t \in 1\dots\tau} (for \code{"one_step"}).}
#'     \item{n_iter}{Number of iterations completed (for \code{"one_step"}).}
#'     \item{converged}{Logical convergence indicator (for \code{"one_step"}).}
#'     \item{targeting_method}{Character string indicating targeting method used.}
#'   }
#' @export
#' @importFrom dplyr mutate if_else
#' @importFrom stats binomial coef glm glm.fit glm.control qlogis plogis
target_hazard <- function(panel_w,
                          pred_Q_init = NULL,
                          grace_period = 10,
                          targeting_method = c("pooled_glm", "one_step"),
                          base_ord = NULL,
                          strat = NULL,
                          cf_w = NULL,
                          logit_ctrl_init = NULL,
                          logit_tr0_init = NULL,
                          logit_tr1_init = NULL,
                          tau = 30,
                          max_iter = 15L,
                          tol = 0.005,
                          damping = 0.5) {

  targeting_method <- match.arg(targeting_method)

  if (targeting_method == "one_step") {
    return(target_hazard_onestep(
      panel_w = panel_w,
      base_ord = base_ord,
      strat = strat,
      cf_w = cf_w,
      logit_ctrl_init = logit_ctrl_init,
      logit_tr0_init = logit_tr0_init,
      logit_tr1_init = logit_tr1_init,
      grace_period = grace_period,
      tau = tau,
      max_iter = max_iter,
      tol = tol,
      damping = damping
    ))
  }

  # Default "pooled_glm" implementation
  panel_w <- panel_w %>%
    dplyr::mutate(
      Q_init_hazard = bound_precision(pred_Q_init),
      Q_init_logit = stats::qlogis(Q_init_hazard),
      H_ctrl = dplyr::if_else(arm == "control_no_tx", w_trim, 0),
      H_treat_pre = dplyr::if_else(arm == "treat_tx" & t <= grace_period, w_trim, 0),
      H_treat_post = dplyr::if_else(arm == "treat_tx" & t > grace_period, w_trim, 0)
    )

  m_fluctuate <- suppressWarnings(stats::glm(
    event_cc_day ~ -1 + H_ctrl + H_treat_pre + H_treat_post + offset(Q_init_logit),
    family = stats::binomial(),
    data = panel_w
  ))

  coefs <- stats::coef(m_fluctuate)

  eps_ctrl <- if ("H_ctrl" %in% names(coefs) && !is.na(coefs["H_ctrl"])) coefs["H_ctrl"] else 0
  eps_treat_pre <- if ("H_treat_pre" %in% names(coefs) && !is.na(coefs["H_treat_pre"])) coefs["H_treat_pre"] else 0
  eps_treat_post <- if ("H_treat_post" %in% names(coefs) && !is.na(coefs["H_treat_post"])) coefs["H_treat_post"] else 0

  names(eps_ctrl) <- NULL
  names(eps_treat_pre) <- NULL
  names(eps_treat_post) <- NULL

  return(list(
    panel_w = panel_w,
    m_fluctuate = m_fluctuate,
    eps_ctrl = eps_ctrl,
    eps_treat_pre = eps_treat_pre,
    eps_treat_post = eps_treat_post,
    eps = c(eps_ctrl = eps_ctrl, eps_treat_pre = eps_treat_pre, eps_treat_post = eps_treat_post),
    targeting_method = "pooled_glm"
  ))
}

#' Universal One-Step Survival TMLE Hazard Targeting
#'
#' @description Updates initial hazard surfaces along the joint pathwise efficient influence curve
#'   gradient across all time horizons \eqn{t_0 \in 1\dots\tau} simultaneously, solving the EIC
#'   estimating equations for both arms without pooling across time periods.
#'
#' @inheritParams target_hazard
#' @export
target_hazard_onestep <- function(panel_w,
                                  base_ord,
                                  strat,
                                  cf_w,
                                  logit_ctrl_init,
                                  logit_tr0_init,
                                  logit_tr1_init,
                                  grace_period = 10,
                                  tau = 30,
                                  max_iter = 15L,
                                  tol = 0.005,
                                  damping = 0.5) {
  n_pat <- nrow(base_ord)
  id_idx <- match(panel_w$id, base_ord$id)
  w_eic_vec <- if ("w_unstab_trim" %in% names(panel_w)) panel_w$w_unstab_trim else panel_w$w_trim

  logit_c <- logit_ctrl_init
  logit_0 <- logit_tr0_init
  logit_1 <- logit_tr1_init

  is_ctrl <- panel_w$arm == "control_no_tx"
  is_tr <- panel_w$arm == "treat_tx"
  idx_c <- which(is_ctrl)
  idx_t <- which(is_tr)
  rc_mat_c <- cbind(id_idx[idx_c], panel_w$t[idx_c])
  rc_mat_t <- cbind(id_idx[idx_t], panel_w$t[idx_t])

  t_c <- panel_w$t[idx_c]
  prev_t_c <- pmax(1L, t_c - 1L)
  t_t <- panel_w$t[idx_t]
  prev_t_t <- pmax(1L, t_t - 1L)
  sw_t <- panel_w$t_switch[idx_t]
  in_post <- t_t >= sw_t
  idx_t_pre <- idx_t[!in_post]
  idx_t_post <- idx_t[in_post]

  cf_w_pre <- matrix(1, nrow = n_pat, ncol = tau)
  cf_w_post <- matrix(0, nrow = n_pat, ncol = tau)
  for (idx in seq_along(strat$k_vals)) {
    cf_w_post <- cf_w_post + sweep(cf_w$cf_H_treat_list[[idx]], 1L, strat$w_mat[, idx], "*")
  }

  fit_1d_step <- function(y, g, offset) {
    if (length(y) == 0L || max(abs(g)) < 1e-10) return(0)
    m <- suppressWarnings(stats::glm.fit(
      x = matrix(g, ncol = 1),
      y = y,
      offset = offset,
      family = stats::binomial(),
      control = stats::glm.control(maxit = 10, epsilon = 1e-4)
    ))
    co <- as.numeric(stats::coef(m)[1])
    if (is.na(co)) 0 else pmin(pmax(co, -5.0), 5.0)
  }

  converged <- FALSE
  c_rd_final <- numeric(tau)
  iter_final <- 0L

  for (iter in seq_len(max_iter)) {
    iter_final <- iter

    # 1. S0
    Q_c <- stats::plogis(logit_c)
    S0 <- 1 - Q_c
    if (tau > 1L) {
      for (tt in 2:tau) S0[, tt] <- S0[, tt - 1L] * S0[, tt]
    }

    # 2. S1 across policy distribution
    S1 <- matrix(0, nrow = n_pat, ncol = tau)
    for (idx in seq_along(strat$k_vals)) {
      wk <- strat$w_mat[, idx]
      if (all(wk <= 0)) next
      k <- strat$k_vals[idx]
      lk <- logit_1
      if (k > 1L) {
        pre_t <- seq_len(min(k - 1L, tau))
        lk[, pre_t] <- logit_0[, pre_t, drop = FALSE]
      }
      qk <- stats::plogis(lk)
      sk <- 1 - qk
      if (tau > 1L) {
        for (tt in 2:tau) sk[, tt] <- sk[, tt - 1L] * sk[, tt]
      }
      S1 <- S1 + sweep(sk, 1L, wk, "*")
    }

    # 3. Predict Q on panel rows
    Q_panel <- numeric(nrow(panel_w))
    S_prev_panel <- numeric(nrow(panel_w))
    offset_panel <- numeric(nrow(panel_w))

    Q_panel[idx_c] <- Q_c[rc_mat_c]
    offset_panel[idx_c] <- logit_c[rc_mat_c]
    S_prev_panel[idx_c] <- ifelse(t_c == 1L, 1.0, S0[cbind(id_idx[idx_c], prev_t_c)])

    off_t <- ifelse(in_post, logit_1[rc_mat_t], logit_0[rc_mat_t])
    offset_panel[idx_t] <- off_t
    Q_panel[idx_t] <- stats::plogis(off_t)
    S_prev_panel[idx_t] <- ifelse(t_t == 1L, 1.0, S1[cbind(id_idx[idx_t], prev_t_t)])

    # 4. Martingale residual steps & cumulative martingales
    denom <- pmax(S_prev_panel * (1 - Q_panel), 1e-4)
    M_step <- w_eic_vec * (panel_w$event_cc_day - Q_panel) / denom

    M_mat_0 <- matrix(0, nrow = n_pat, ncol = tau)
    M_mat_0[rc_mat_c] <- M_step[idx_c]
    cum_M_0 <- t(apply(M_mat_0, 1, cumsum))

    M_mat_1 <- matrix(0, nrow = n_pat, ncol = tau)
    M_mat_1[rc_mat_t] <- M_step[idx_t]
    cum_M_1 <- t(apply(M_mat_1, 1, cumsum))

    c_0 <- colMeans(S0 * cum_M_0)
    c_1 <- colMeans(S1 * cum_M_1)
    c_rd <- c_1 - c_0
    c_rd_final <- c_rd

    if (max(abs(c_rd)) <= tol) {
      converged <- TRUE
      break
    }

    # 5. Gradient & 1D GLM update steps
    # Control
    S0_weighted <- sweep(S0, 2L, c_0, "*")
    V0 <- t(apply(S0_weighted[, tau:1, drop = FALSE], 1, cumsum))[, tau:1, drop = FALSE]
    G_panel_c <- (w_eic_vec[idx_c] / denom[idx_c]) * V0[rc_mat_c]
    eps_c <- fit_1d_step(panel_w$event_cc_day[idx_c], G_panel_c, offset_panel[idx_c]) * damping

    S0_prev_grid <- cbind(1, S0[, -tau, drop = FALSE])
    denom_grid_c <- pmax(S0_prev_grid * (1 - Q_c), 1e-4)
    G_grid_c <- (cf_w$cf_H_ctrl / denom_grid_c) * V0
    logit_c <- logit_c + eps_c * G_grid_c

    # Treated
    S1_weighted <- sweep(S1, 2L, c_1, "*")
    V1 <- t(apply(S1_weighted[, tau:1, drop = FALSE], 1, cumsum))[, tau:1, drop = FALSE]

    if (length(idx_t_pre) > 0L) {
      rc_pre <- cbind(id_idx[idx_t_pre], panel_w$t[idx_t_pre])
      G_pre <- (w_eic_vec[idx_t_pre] / denom[idx_t_pre]) * V1[rc_pre]
      eps_pre <- fit_1d_step(panel_w$event_cc_day[idx_t_pre], G_pre, offset_panel[idx_t_pre]) * damping

      S1_prev_grid <- cbind(1, S1[, -tau, drop = FALSE])
      denom_grid_0 <- pmax(S1_prev_grid * (1 - stats::plogis(logit_0)), 1e-4)
      G_grid_0 <- (cf_w_pre / denom_grid_0) * V1
      logit_0 <- logit_0 + eps_pre * G_grid_0
    }

    if (length(idx_t_post) > 0L) {
      rc_post <- cbind(id_idx[idx_t_post], panel_w$t[idx_t_post])
      G_post <- (w_eic_vec[idx_t_post] / denom[idx_t_post]) * V1[rc_post]
      eps_post <- fit_1d_step(panel_w$event_cc_day[idx_t_post], G_post, offset_panel[idx_t_post]) * damping

      S1_prev_grid <- cbind(1, S1[, -tau, drop = FALSE])
      denom_grid_1 <- pmax(S1_prev_grid * (1 - stats::plogis(logit_1)), 1e-4)
      G_grid_1 <- (cf_w_post / denom_grid_1) * V1
      logit_1 <- logit_1 + eps_post * G_grid_1
    }
  }

  Q_star_ctrl_mat <- stats::plogis(logit_c)
  Q_star_tr0_mat  <- stats::plogis(logit_0)
  Q_star_tr1_mat  <- stats::plogis(logit_1)

  Q_panel[idx_c] <- Q_star_ctrl_mat[rc_mat_c]
  Q_panel[idx_t] <- ifelse(in_post,
                           Q_star_tr1_mat[rc_mat_t],
                           Q_star_tr0_mat[rc_mat_t])

  panel_w <- panel_w %>%
    dplyr::mutate(
      Q_star_hazard = bound_precision(Q_panel),
      H_ctrl = dplyr::if_else(arm == "control_no_tx", w_trim, 0),
      H_treat_pre = dplyr::if_else(arm == "treat_tx" & t <= grace_period, w_trim, 0),
      H_treat_post = dplyr::if_else(arm == "treat_tx" & t > grace_period, w_trim, 0)
    )

  list(
    panel_w = panel_w,
    m_fluctuate = NULL,
    eps_ctrl = 0,
    eps_treat_pre = 0,
    eps_treat_post = 0,
    eps = c(eps_ctrl = 0, eps_treat_pre = 0, eps_treat_post = 0),
    Q_star_ctrl_mat = Q_star_ctrl_mat,
    Q_star_tr0_mat = Q_star_tr0_mat,
    Q_star_tr1_mat = Q_star_tr1_mat,
    c_rd = c_rd_final,
    n_iter = iter_final,
    converged = converged,
    targeting_method = "one_step"
  )
}
