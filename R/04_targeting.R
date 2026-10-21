#' Decoupled Clever Covariates and TMLE Fluctuation Step
#'
#' @description Constructs decoupled clever covariates (separating control, pre-grace, and post-grace
#'   time windows) and fits the TMLE fluctuation logistic regression submodel with offset.
#'
#' @param panel_w Longitudinal person-time panel from \code{\link{compute_ipcw_weights}}.
#' @param pred_Q_init Initial hazard predictions from \code{\link{fit_q_models}}.
#' @param grace_period Numeric grace period duration (default 10).
#'
#' @return A list containing:
#'   \describe{
#'     \item{panel_w}{Panel updated with initial hazard logits and clever covariates.}
#'     \item{m_fluctuate}{Fitted fluctuation regression object.}
#'     \item{eps_ctrl}{Estimated fluctuation coefficient for control arm.}
#'     \item{eps_treat_pre}{Estimated fluctuation coefficient for treatment arm (pre-grace period).}
#'     \item{eps_treat_post}{Estimated fluctuation coefficient for treatment arm (post-grace period).}
#'   }
#' @export
#' @importFrom dplyr mutate if_else
#' @importFrom stats binomial coef glm qlogis
target_hazard <- function(panel_w,
                          pred_Q_init,
                          grace_period = 10) {
  
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
    eps_treat_post = eps_treat_post
  ))
}

