#' Counterfactual G-Computation and Marginal Survival Curves
#'
#' @description Constructs the counterfactual expanded person-time dataset ("ghost population"),
#'   analytically integrates counterfactual survival across the target treatment-arm initiation
#'   policy (\code{"conditional_complier"}, \code{"empirical_complier"}, or \code{"natural_grace"}),
#'   applies the TMLE targeting fluctuation (or evaluates pure G-computation),
#'   and marginalizes across the baseline population to produce counterfactual survival curves and risk differences.
#'
#' @param data Raw standardized data frame.
#' @param clones Cloned data frame from \code{\link{create_clones}}.
#' @param m_Q_ctrl Outcome hazard model for control arm.
#' @param m_Q_treat Outcome hazard model for treated arm.
#' @param m_denom_ctrl IPCW denominator model for control arm.
#' @param m_denom_treat IPCW denominator model for treated arm.
#' @param m_num_ctrl IPCW numerator model for control arm.
#' @param m_num_treat IPCW numerator model for treated arm.
#' @param eps_ctrl Fluctuation coefficient for control arm.
#' @param eps_treat_pre Fluctuation coefficient for treated arm (pre-grace period).
#' @param eps_treat_post Fluctuation coefficient for treated arm (post-grace period).
#' @param grace_period Grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param covars Character vector of baseline covariates.
#' @param q_lim Trimming bounds for weights (list by arm or length-2 numeric vector).
#' @param use_sl Logical indicating if SuperLearner was used.
#' @param do_targeting Logical. If \code{TRUE} (default), applies the TMLE fluctuation update.
#'   If \code{FALSE}, evaluates pure G-computation without targeting.
#' @param grace_surv Optional precomputed grace-period survival grid from \code{\link{compute_ipcw_weights}}.
#' @param tx_strategy Treatment-arm initiation policy (\code{"conditional_complier"},
#'   \code{"empirical_complier"}, or \code{"natural_grace"}).
#' @param p_star Optional precomputed empirical complier initiation distribution across days \code{1..m}.
#' @param align_grace_weights Deprecated compatibility flag (mapped to \code{"natural_grace"} if \code{FALSE}).
#' @param targeting_method Targeting algorithm: \code{"pooled_glm"} (default) or \code{"one_step"}.
#' @param targ_res Optional targeting result list from \code{\link{target_hazard}}.
#' @param logit_ctrl_mat Optional precomputed baseline control logit hazard matrix.
#' @param logit_tr0_mat Optional precomputed baseline treated untreated-state logit hazard matrix.
#' @param logit_tr1_mat Optional precomputed baseline treated treated-state logit hazard matrix.
#' @param strat Optional precomputed strategy weights list.
#' @param cf_w Optional precomputed counterfactual weight matrices list.
#'
#' @return A list containing:
#'   \describe{
#'     \item{cf_data}{Longitudinal counterfactual dataset with individual survival probabilities.}
#'     \item{tmle_surv}{Summary tibble of marginal survival by arm and time.}
#'     \item{res_tmle}{Wide tibble containing marginal survival for both arms and the marginal risk difference \code{rd}.}
#'     \item{strat}{Treatment strategy mixture weights list.}
#'     \item{cf_w}{Counterfactual weight matrices.}
#'   }
#' @export
#' @importFrom dplyr select mutate filter arrange all_of distinct
#' @importFrom tidyr uncount pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats predict plogis qlogis
marginalize_gcomp <- function(data,
                              clones,
                              m_Q_ctrl,
                              m_Q_treat,
                              m_denom_ctrl,
                              m_denom_treat,
                              m_num_ctrl,
                              m_num_treat,
                              eps_ctrl = 0,
                              eps_treat_pre = 0,
                              eps_treat_post = 0,
                              grace_period = 10,
                              tau = 30,
                              covars = c("x1"),
                              q_lim = c(0.1, 10),
                              use_sl = FALSE,
                              do_targeting = TRUE,
                              grace_surv = NULL,
                              tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace"),
                              p_star = NULL,
                              align_grace_weights = TRUE,
                              targeting_method = c("pooled_glm", "one_step"),
                              targ_res = NULL,
                              logit_ctrl_mat = NULL,
                              logit_tr0_mat = NULL,
                              logit_tr1_mat = NULL,
                              strat = NULL,
                              cf_w = NULL) {
  if (missing(tx_strategy) && !isTRUE(align_grace_weights)) {
    tx_strategy <- "natural_grace"
  } else {
    tx_strategy <- match.arg(tx_strategy)
  }
  targeting_method <- match.arg(targeting_method)

  # 1. Baseline patient cohort (ordered by id to align with grace_surv matrices)
  base_ord <- data %>%
    dplyr::select(id, dplyr::all_of(covars)) %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::arrange(id)
  n_pat <- nrow(base_ord)

  if (is.null(grace_surv)) {
    grace_surv <- compute_grace_survival_grid(
      data_patients = base_ord,
      covars = covars,
      grace_period = grace_period,
      dt = 1,
      m_denom_ctrl = m_denom_ctrl,
      m_num_ctrl = m_num_ctrl,
      use_sl = use_sl
    )
  }

  if (tx_strategy == "empirical_complier" && is.null(p_star)) {
    p_star <- compute_empirical_complier_probs(clones, grace_period = grace_period, n_grace_int = grace_surv$n_grace_int)
  }

  # 2. Expanded (id x t) baseline prediction grid (ordered by id, t)
  if (is.null(logit_ctrl_mat) || is.null(logit_tr0_mat) || is.null(logit_tr1_mat)) {
    cf_grid <- base_ord %>%
      tidyr::uncount(weights = tau, .id = "t") %>%
      dplyr::mutate(
        t_sq = t^2,
        post_grace = dplyr::if_else(t > grace_period, 1, 0)
      )

    # Predict base hazard surfaces ONCE: control (A=0), treated untreated-state (A=0), treated treated-state (A=1)
    cf_ctrl_grid <- cf_grid %>% dplyr::mutate(is_treated = 0)
    cf_tr0_grid  <- cf_grid %>% dplyr::mutate(is_treated = 0)
    cf_tr1_grid  <- cf_grid %>% dplyr::mutate(is_treated = 1)

    if (use_sl) {
      X_ctrl <- cf_ctrl_grid %>% dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace) %>% as.data.frame()
      X_tr0  <- cf_tr0_grid %>% dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()
      X_tr1  <- cf_tr1_grid %>% dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()

      q_ctrl_init <- stats::predict(m_Q_ctrl, newdata = X_ctrl, onlySL = TRUE)$pred[, 1]
      q_tr0_init  <- stats::predict(m_Q_treat, newdata = X_tr0, onlySL = TRUE)$pred[, 1]
      q_tr1_init  <- stats::predict(m_Q_treat, newdata = X_tr1, onlySL = TRUE)$pred[, 1]
    } else {
      q_ctrl_init <- stats::predict(m_Q_ctrl, newdata = cf_ctrl_grid, type = "response")
      q_tr0_init  <- stats::predict(m_Q_treat, newdata = cf_tr0_grid, type = "response")
      q_tr1_init  <- stats::predict(m_Q_treat, newdata = cf_tr1_grid, type = "response")
    }

    logit_ctrl_mat <- matrix(stats::qlogis(bound_precision(q_ctrl_init)), nrow = n_pat, ncol = tau, byrow = TRUE)
    logit_tr0_mat  <- matrix(stats::qlogis(bound_precision(q_tr0_init)),  nrow = n_pat, ncol = tau, byrow = TRUE)
    logit_tr1_mat  <- matrix(stats::qlogis(bound_precision(q_tr1_init)),  nrow = n_pat, ncol = tau, byrow = TRUE)
  }

  # 3. Strategy weights & counterfactual clever-covariate matrices
  if (is.null(strat)) {
    strat <- compute_tx_strategy_weights(grace_surv, tx_strategy = tx_strategy, p_star = p_star)
  }
  if (is.null(cf_w)) {
    cf_w <- compute_cf_weight_matrices(
      grace_surv = grace_surv,
      data_patients = base_ord,
      covars = covars,
      grace_period = grace_period,
      tau = tau,
      tx_strategy = tx_strategy,
      k_vals = strat$k_vals,
      q_lim = q_lim,
      m_denom_treat = m_denom_treat,
      m_num_treat = m_num_treat,
      use_sl = use_sl,
      p_star = p_star
    )
  }

  # 4. Control counterfactual survival matrix (n_pat x tau)
  q_star_ctrl <- if (targeting_method == "one_step" && !is.null(targ_res$Q_star_ctrl_mat)) {
    targ_res$Q_star_ctrl_mat
  } else if (do_targeting) {
    stats::plogis(logit_ctrl_mat + eps_ctrl * cf_w$cf_H_ctrl)
  } else {
    stats::plogis(logit_ctrl_mat)
  }
  S0_mat <- 1 - q_star_ctrl
  if (tau > 1L) {
    for (tt in 2:tau) S0_mat[, tt] <- S0_mat[, tt - 1L] * S0_mat[, tt]
  }

  # 5. Treated counterfactual survival matrix via exact analytical integration over k
  S1_mat <- matrix(0, nrow = n_pat, ncol = tau)
  eps_t_vec <- ifelse(seq_len(tau) > grace_period, eps_treat_post, eps_treat_pre)

  for (idx in seq_along(strat$k_vals)) {
    wk <- strat$w_mat[, idx]
    if (all(wk <= 0)) next
    k <- strat$k_vals[idx]

    logit_k <- logit_tr1_mat
    pre_t <- seq_len(min(k - 1L, tau))
    if (length(pre_t) > 0L) {
      logit_k[, pre_t] <- logit_tr0_mat[, pre_t, drop = FALSE]
    }

    q_star_k <- if (targeting_method == "one_step" && !is.null(targ_res$Q_star_tr1_mat)) {
      qk <- targ_res$Q_star_tr1_mat
      if (length(pre_t) > 0L) {
        qk[, pre_t] <- targ_res$Q_star_tr0_mat[, pre_t, drop = FALSE]
      }
      qk
    } else if (do_targeting) {
      stats::plogis(logit_k + sweep(cf_w$cf_H_treat_list[[idx]], 2L, eps_t_vec, "*"))
    } else {
      stats::plogis(logit_k)
    }
    S_k <- 1 - q_star_k
    if (tau > 1L) {
      for (tt in 2:tau) S_k[, tt] <- S_k[, tt - 1L] * S_k[, tt]
    }
    S1_mat <- S1_mat + sweep(S_k, 1L, wk, "*")
  }

  # 6. Assemble patient-level counterfactual tibble matching compute_eic_inference()
  ids_vec <- rep(base_ord$id, each = tau)
  t_vec   <- rep(seq_len(tau), times = n_pat)

  cf_data <- dplyr::bind_rows(
    tibble::tibble(id = ids_vec, arm = "control_no_tx", t = t_vec, indiv_survival = as.vector(t(S0_mat))),
    tibble::tibble(id = ids_vec, arm = "treat_tx",      t = t_vec, indiv_survival = as.vector(t(S1_mat)))
  )

  mean_s0 <- colMeans(S0_mat)
  mean_s1 <- colMeans(S1_mat)

  tmle_surv <- tibble::tibble(
    t = rep(seq_len(tau), 2L),
    arm = c(rep("control_no_tx", tau), rep("treat_tx", tau)),
    survival = c(mean_s0, mean_s1)
  )

  res_tmle <- tibble::tibble(
    t = seq_len(tau),
    control_no_tx = mean_s0,
    treat_tx = mean_s1,
    rd = mean_s1 - mean_s0
  )

  return(list(
    cf_data = cf_data,
    tmle_surv = tmle_surv,
    res_tmle = res_tmle,
    strat = strat,
    cf_w = cf_w,
    logit_ctrl_mat = logit_ctrl_mat,
    logit_tr0_mat = logit_tr0_mat,
    logit_tr1_mat = logit_tr1_mat,
    base_ord = base_ord
  ))
}

