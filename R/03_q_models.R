#' Fit Stratified Initial Outcome Hazard Models (Q-Models)
#'
#' @description Fits stratified discrete-time hazard models for each treatment arm
#'   using generalized linear models (GLM) or non-parametric machine learning ensembling (\code{SuperLearner}).
#'
#' @param panel_w Longitudinal person-time panel from \code{\link{compute_ipcw_weights}}.
#' @param covars Character vector of baseline covariate names.
#' @param q_form Optional formula string for the initial hazard model. Defaults to linear additive terms of \code{covars}.
#' @param use_sl Logical indicating whether to use \code{SuperLearner} for hazard modeling.
#' @param sl_lib Character vector of SuperLearner algorithm wrappers.
#'
#' @return A list containing:
#'   \describe{
#'     \item{m_Q_ctrl}{Fitted outcome hazard model for the control arm.}
#'     \item{m_Q_treat}{Fitted outcome hazard model for the treated arm.}
#'     \item{pred_Q_ctrl}{Initial hazard predictions under control across \code{panel_w}.}
#'     \item{pred_Q_treat}{Initial hazard predictions under treatment across \code{panel_w}.}
#'     \item{pred_Q_init}{Initial hazard predictions assigned according to observed arm.}
#'   }
#' @export
#' @importFrom dplyr filter select mutate if_else all_of
#' @importFrom stats as.formula binomial glm predict
fit_q_models <- function(panel_w,
                         covars = c("x1"),
                         q_form = NULL,
                         use_sl = FALSE,
                         sl_lib = c("SL.mean", "SL.glm", "SL.earth.int", "SL.ranger")) {
  
  if (is.null(q_form)) {
    q_form <- paste(covars, collapse = " + ")
  }
  
  panel_ctrl <- panel_w %>% dplyr::filter(arm == "control_no_tx")
  panel_treat <- panel_w %>% dplyr::filter(arm == "treat_tx")
  
  if (use_sl) {
    if (!require("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required when use_sl = TRUE.")
    }
    register_sl_earth_int()
    
    # In control arm, is_treated is identically 0 by design (protocol censoring prevents treated person-time)
    X_Q_ctrl <- panel_ctrl %>% 
      dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace) %>% 
      as.data.frame()
    X_Q_treat <- panel_treat %>% 
      dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace, is_treated) %>% 
      as.data.frame()
    
    m_Q_ctrl <- SuperLearner::SuperLearner(
      Y = panel_ctrl$event_cc_day,
      X = X_Q_ctrl,
      family = stats::binomial(),
      SL.library = sl_lib,
      cvControl = list(V = 2)
    )
    
    m_Q_treat <- SuperLearner::SuperLearner(
      Y = panel_treat$event_cc_day,
      X = X_Q_treat,
      family = stats::binomial(),
      SL.library = sl_lib,
      cvControl = list(V = 2)
    )
    
    X_Q_full_ctrl <- panel_w %>% 
      dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace) %>% 
      as.data.frame()
    X_Q_full_treat <- panel_w %>% 
      dplyr::select(dplyr::all_of(covars), t, t_sq, post_grace, is_treated) %>% 
      as.data.frame()
    
    pred_Q_ctrl <- stats::predict(m_Q_ctrl, newdata = X_Q_full_ctrl, onlySL = TRUE)$pred[, 1]
    pred_Q_treat <- stats::predict(m_Q_treat, newdata = X_Q_full_treat, onlySL = TRUE)$pred[, 1]
  } else {
    Q_formula_ctrl <- stats::as.formula(paste("event_cc_day ~", q_form, "+ t + t_sq + post_grace"))
    Q_formula_treat <- stats::as.formula(paste("event_cc_day ~", q_form, "+ t + t_sq + post_grace + is_treated"))
    
    m_Q_ctrl <- stats::glm(Q_formula_ctrl, family = stats::binomial(), data = panel_ctrl)
    m_Q_treat <- stats::glm(Q_formula_treat, family = stats::binomial(), data = panel_treat)
    
    pred_Q_ctrl <- stats::predict(m_Q_ctrl, newdata = panel_w, type = "response")
    pred_Q_treat <- stats::predict(m_Q_treat, newdata = panel_w, type = "response")
  }
  
  pred_Q_init <- dplyr::if_else(panel_w$arm == "treat_tx", pred_Q_treat, pred_Q_ctrl)
  
  return(list(
    m_Q_ctrl = m_Q_ctrl,
    m_Q_treat = m_Q_treat,
    pred_Q_ctrl = pred_Q_ctrl,
    pred_Q_treat = pred_Q_treat,
    pred_Q_init = pred_Q_init
  ))
}

