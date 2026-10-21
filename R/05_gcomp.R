#' Counterfactual G-Computation and Marginal Survival Curves
#'
#' @description Constructs the counterfactual expanded person-time dataset ("ghost population"),
#'   imputes dynamic regime treatment initiation times for non-compliers by sampling from empirical compliers,
#'   updates hazard predictions using the TMLE targeting fluctuation (or evaluates pure G-computation),
#'   and marginalizes across the population to produce counterfactual survival curves and risk differences.
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
#' @param q_lim Trimming bounds for weights.
#' @param use_sl Logical indicating if SuperLearner was used.
#' @param do_targeting Logical. If \code{TRUE} (default), applies the TMLE fluctuation update.
#'   If \code{FALSE}, evaluates pure G-computation without targeting.
#'
#' @return A list containing:
#'   \describe{
#'     \item{cf_data}{Longitudinal counterfactual dataset with individual survival probabilities.}
#'     \item{tmle_surv}{Summary tibble of marginal survival by arm and time.}
#'     \item{res_tmle}{Wide tibble containing marginal survival for both arms and the marginal risk difference \code{rd}.}
#'   }
#' @export
#' @importFrom dplyr select mutate filter pull arrange group_by ungroup summarise if_else case_when all_of
#' @importFrom tidyr uncount pivot_wider
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
                              do_targeting = TRUE) {
  
  # 1. Complier switch times empirical distribution
  complier_times <- clones %>%
    dplyr::filter(arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period) %>%
    dplyr::pull(t_switch)
  
  if (length(complier_times) == 0) complier_times <- grace_period
  
  # 2. Counterfactual baseline
  cf_base <- data %>%
    dplyr::select(id, tx_ever, t_switch, dplyr::all_of(covars)) %>%
    dplyr::mutate(sampled_t = sample(complier_times, dplyr::n(), replace = TRUE))
  
  # 3. Longitudinal counterfactual expansion
  cf_data <- cf_base %>%
    tidyr::uncount(weights = tau, .id = "t") %>%
    dplyr::mutate(t_sq = t^2, post_grace = dplyr::if_else(t > grace_period, 1, 0)) %>%
    tidyr::uncount(weights = 2, .id = "clone_arm") %>%
    dplyr::mutate(
      arm = dplyr::if_else(clone_arm == 1L, "control_no_tx", "treat_tx"),
      cf_t_switch = dplyr::case_when(
        arm == "control_no_tx" ~ 999,
        arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period ~ t_switch,
        arm == "treat_tx" ~ sampled_t
      ),
      is_treated = dplyr::if_else(t >= cf_t_switch, 1, 0)
    )
  
  # 4. Predict initial counterfactual hazards & weights
  if (use_sl) {
    X_cf_ctrl <- cf_data %>% dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace) %>% as.data.frame()
    X_cf_treat <- cf_data %>% dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()
    X_cf_denom_ctrl <- cf_data %>% dplyr::select(dplyr::all_of(covars), t, t_sq) %>% as.data.frame()
    X_cf_treat_only <- cf_data %>% dplyr::select(dplyr::all_of(covars)) %>% as.data.frame()
    
    cf_pred_Q_ctrl <- stats::predict(m_Q_ctrl, newdata = X_cf_ctrl, onlySL = TRUE)$pred[, 1]
    cf_pred_Q_treat <- stats::predict(m_Q_treat, newdata = X_cf_treat, onlySL = TRUE)$pred[, 1]
    cf_pred_Q_init <- dplyr::if_else(cf_data$arm == "treat_tx", cf_pred_Q_treat, cf_pred_Q_ctrl)
    
    cf_pred_ctrl_denom <- stats::predict(m_denom_ctrl, newdata = X_cf_denom_ctrl, onlySL = TRUE)$pred[, 1]
    cf_pred_treat_denom <- stats::predict(m_denom_treat, newdata = X_cf_treat_only, onlySL = TRUE)$pred[, 1]
  } else {
    cf_pred_Q_ctrl <- stats::predict(m_Q_ctrl, newdata = cf_data, type = "response")
    cf_pred_Q_treat <- stats::predict(m_Q_treat, newdata = cf_data, type = "response")
    cf_pred_Q_init <- dplyr::if_else(cf_data$arm == "treat_tx", cf_pred_Q_treat, cf_pred_Q_ctrl)
    
    cf_pred_ctrl_denom <- stats::predict(m_denom_ctrl, newdata = cf_data, type = "response")
    cf_pred_treat_denom <- stats::predict(m_denom_treat, newdata = cf_data, type = "response")
  }
  
  # 5. Calculate counterfactual weights and clever covariates
  cf_data <- cf_data %>%
    dplyr::mutate(
      cf_Q_init_hazard = bound_precision(cf_pred_Q_init),
      cf_Q_init_logit  = stats::qlogis(cf_Q_init_hazard),
      cf_p_ctrl_denom  = bound_precision(cf_pred_ctrl_denom),
      cf_p_ctrl_num    = bound_precision(stats::predict(m_num_ctrl, newdata = ., type = "response")),
      cf_p_treat_denom = bound_precision(cf_pred_treat_denom),
      cf_p_treat_num   = bound_precision(stats::predict(m_num_treat, newdata = ., type = "response"))
    ) %>%
    dplyr::mutate(
      w_inc = dplyr::case_when(
        arm == "control_no_tx" ~ (1 - cf_p_ctrl_num) / (1 - cf_p_ctrl_denom),
        arm == "treat_tx" ~ dplyr::case_when(t == grace_period + 1L ~ cf_p_treat_num / cf_p_treat_denom, TRUE ~ 1)
      )
    ) %>%
    dplyr::arrange(id, arm, t) %>%
    dplyr::group_by(id, arm) %>%
    dplyr::mutate(cf_H_raw = cumprod(w_inc)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cf_H = pmin(pmax(cf_H_raw, q_lim[1]), q_lim[2]),
      cf_H_ctrl = dplyr::if_else(arm == "control_no_tx", cf_H, 0),
      cf_H_treat_pre = dplyr::if_else(arm == "treat_tx" & t <= grace_period, cf_H, 0),
      cf_H_treat_post = dplyr::if_else(arm == "treat_tx" & t > grace_period, cf_H, 0)
    )
  
  # 6. Fluctuate hazard (or keep unfluctuated for pure g-computation)
  if (do_targeting) {
    cf_data <- cf_data %>%
      dplyr::mutate(
        Q_star_hazard = stats::plogis(cf_Q_init_logit + 
                                        eps_ctrl * cf_H_ctrl + 
                                        eps_treat_pre * cf_H_treat_pre + 
                                        eps_treat_post * cf_H_treat_post)
      )
  } else {
    cf_data <- cf_data %>%
      dplyr::mutate(Q_star_hazard = cf_Q_init_hazard)
  }
  
  # 7. Cumulative Individual Survival & Population Marginalization
  cf_data <- cf_data %>%
    dplyr::group_by(id, arm) %>%
    dplyr::mutate(indiv_survival = cumprod(1 - Q_star_hazard)) %>%
    dplyr::ungroup()
  
  tmle_surv <- cf_data %>%
    dplyr::group_by(arm, t) %>%
    dplyr::summarise(survival = mean(indiv_survival), .groups = "drop")
  
  res_tmle <- tmle_surv %>%
    dplyr::select(t, arm, survival) %>%
    tidyr::pivot_wider(names_from = arm, values_from = survival) %>%
    dplyr::mutate(rd = treat_tx - control_no_tx)
  
  return(list(
    cf_data = cf_data,
    tmle_surv = tmle_surv,
    res_tmle = res_tmle
  ))
}

