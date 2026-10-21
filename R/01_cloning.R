#' Create Clones and Apply Protocol Censoring for Target Trial Emulation
#'
#' @description Expands observational time-to-event data by cloning each individual into
#'   both treatment strategy arms (Treatment within grace period vs. Never Treat), and applies
#'   artificial protocol censoring when an individual deviates from their assigned strategy.
#'
#' @param data A data.frame containing individual-level observational data with standardized columns
#'   \code{id}, \code{time}, \code{event}, \code{tx_ever}, \code{t_switch}, and baseline covariates.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric threshold for treatment initiation in the treatment arm (default 10).
#' @param tau Maximum follow-up horizon (default 30).
#'
#' @return A data.frame of cloned person-level records with added variables:
#'   \describe{
#'     \item{arm}{Treatment strategy arm: \code{"control_no_tx"} or \code{"treat_tx"}.}
#'     \item{clone_id}{Unique identifier for each clone (\code{id_arm}).}
#'     \item{dev_time}{Time of protocol deviation, or \code{NA} if compliant.}
#'     \item{time_cc}{Observed clone follow-up time, minimum of true event time, deviation time, and \code{tau}.}
#'     \item{event_cc}{Binary indicator of event occurring on or before clone censoring time.}
#'   }
#' @export
#' @importFrom dplyr select mutate if_else case_when all_of
#' @importFrom tidyr uncount
create_clones <- function(data, 
                          covars = c("x1"), 
                          grace_period = 10, 
                          tau = 30) {
  
  clones <- data %>%
    dplyr::select(id, tx_ever, t_switch, event, time, dplyr::all_of(covars)) %>%
    tidyr::uncount(weights = 2, .id = "clone_arm") %>%
    dplyr::mutate(
      arm = dplyr::if_else(clone_arm == 1L, "control_no_tx", "treat_tx"),
      clone_id = paste0(id, "_", clone_arm),
      dev_time = dplyr::case_when(
        arm == "control_no_tx" & tx_ever == 1L ~ t_switch - 1,
        arm == "treat_tx" & !(tx_ever == 1L & t_switch <= grace_period) & time > grace_period ~ as.numeric(grace_period),
        TRUE ~ NA_real_
      ),
      time_cc = dplyr::if_else(!is.na(dev_time), pmin(time, dev_time), time),
      time_cc = pmin(time_cc, tau),
      event_cc = dplyr::if_else(event == 1L & time <= time_cc, 1L, 0L)
    )
  
  return(clones)
}

