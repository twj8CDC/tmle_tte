#' Fit Compliance Models and Compute Inverse Probability of Censoring Weights (IPCW)
#'
#' @description Fits pooled logistic or SuperLearner compliance models on observational person-time,
#'   expands clones into a discrete person-time panel, and calculates stabilized, trimmed
#'   cumulative inverse probability of censoring weights.
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
#'
#' @return A list containing:
#'   \describe{
#'     \item{panel_w}{Longitudinal person-time panel of cloned observations with trimmed weights \code{w_trim}.}
#'     \item{fit_ccw}{Survival fit object from \code{survival::survfit} using IPCW weights.}
#'     \item{m_denom_ctrl}{Fitted compliance model for the control arm.}
#'     \item{m_denom_treat}{Fitted compliance model for the treated arm.}
#'     \item{m_num_ctrl}{Fitted numerator model for control stabilization.}
#'     \item{m_num_treat}{Fitted numerator model for treated stabilization.}
#'     \item{q_lim}{Weight trimming limits applied.}
#'   }
#' @export
#' @importFrom dplyr select mutate filter arrange group_by ungroup if_else case_when all_of
#' @importFrom tidyr uncount
#' @importFrom stats as.formula binomial glm predict quantile
#' @importFrom survival survfit Surv
compute_ipcw_weights <- function(clones,
                                 covars = c("x1"),
                                 grace_period = 10,
                                 dt = 1,
                                 p_form = NULL,
                                 use_sl = FALSE,
                                 sl_lib = c("SL.mean", "SL.glm", "SL.earth.int", "SL.ranger"),
                                 trim_quantile = c(0.01, 0.99)) {
  
  if (is.null(p_form)) {
    p_form <- paste(covars, collapse = " + ")
  }
  
  if (use_sl) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required when use_sl = TRUE.")
    }
    register_sl_earth_int()
  }
  
  # 1. Observational panel for control compliance hazard
  obs_panel <- clones %>%
    dplyr::filter(clone_arm == 1L) %>%
    dplyr::mutate(
      tmax = dplyr::if_else(tx_ever == 1L, t_switch, time),
      n_int = pmax(ceiling(tmax / dt), 1L)
    ) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      t_sq = t^2,
      event_tx = dplyr::if_else(tx_ever == 1L & t == t_switch, 1L, 0L)
    )
  
  # 2. Treated cohort at grace period cutoff
  treat_day10_data <- clones %>%
    dplyr::filter(arm == "treat_tx", time > grace_period) %>%
    dplyr::mutate(took_tx = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))
  
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
    
    X_treat <- treat_day10_data %>% dplyr::select(dplyr::all_of(covars)) %>% as.data.frame()
    m_denom_treat <- SuperLearner::SuperLearner(
      Y = treat_day10_data$took_tx,
      X = X_treat,
      family = stats::binomial(),
      SL.library = sl_lib,
      cvControl = list(V = 2)
    )
  } else {
    ctrl_formula <- stats::as.formula(paste("event_tx ~", p_form, "+ t + t_sq"))
    treat_formula <- stats::as.formula(paste("took_tx ~", p_form))
    m_denom_ctrl <- suppressWarnings(stats::glm(ctrl_formula, family = stats::binomial(), data = obs_panel))
    m_denom_treat <- suppressWarnings(stats::glm(treat_formula, family = stats::binomial(), data = treat_day10_data))
  }
  
  # 4. Numerator Models for Stabilization
  m_num_ctrl  <- suppressWarnings(stats::glm(event_tx ~ t + t_sq, family = stats::binomial(), data = obs_panel))
  m_num_treat <- suppressWarnings(stats::glm(took_tx ~ 1, family = stats::binomial(), data = treat_day10_data))
  
  # 5. Expand Clones into Longitudinal Person-Time Panel
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
  
  # 6. Predict Denominator Probabilities
  if (use_sl) {
    X_pred_ctrl <- panel %>% dplyr::select(dplyr::all_of(covars), t, t_sq) %>% as.data.frame()
    X_pred_treat <- panel %>% dplyr::select(dplyr::all_of(covars)) %>% as.data.frame()
    
    pred_ctrl_denom <- stats::predict(m_denom_ctrl, newdata = X_pred_ctrl, onlySL = TRUE)$pred[, 1]
    pred_treat_denom <- stats::predict(m_denom_treat, newdata = X_pred_treat, onlySL = TRUE)$pred[, 1]
  } else {
    pred_ctrl_denom <- stats::predict(m_denom_ctrl, newdata = panel, type = "response")
    pred_treat_denom <- stats::predict(m_denom_treat, newdata = panel, type = "response")
  }
  
  # 7. Calculate Incremental & Cumulative Weights
  panel <- panel %>%
    dplyr::mutate(
      p_ctrl_denom = bound_precision(pred_ctrl_denom),
      p_ctrl_num   = bound_precision(stats::predict(m_num_ctrl, newdata = ., type = "response")),
      p_treat_denom = bound_precision(pred_treat_denom),
      p_treat_num   = bound_precision(stats::predict(m_num_treat, newdata = ., type = "response"))
    ) %>%
    dplyr::mutate(
      w_inc = dplyr::case_when(
        arm == "control_no_tx" ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom),
        arm == "treat_tx" ~ dplyr::case_when(
          t == grace_period + 1L ~ p_treat_num / p_treat_denom,
          TRUE ~ 1
        )
      ),
      w_inc_unstab = dplyr::case_when(
        arm == "control_no_tx" ~ 1 / (1 - p_ctrl_denom),
        arm == "treat_tx" ~ dplyr::case_when(
          t == grace_period + 1L ~ 1 / p_treat_denom,
          TRUE ~ 1
        )
      )
    ) %>%
    dplyr::arrange(clone_id, t) %>%
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
      post_grace = dplyr::if_else(t > grace_period, 1, 0)
    )
  
  fit_ccw <- survival::survfit(
    survival::Surv(t0, t1, event_cc_day) ~ arm,
    data = panel_w,
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
    q_lim = q_lim
  ))
}

