#' Internal Utilities for ttetmle
#'
#' @description Helper functions for data formatting, bounding probabilities,
#'   and registering default SuperLearner wrappers.
#'
#' @keywords internal
#' @noRd

# Bound probabilities away from 0 and 1 for numerical stability
bound_precision <- function(x, tol = 1e-4) {
  pmin(pmax(x, tol), 1 - tol)
}

# Standardize column names if necessary
prepare_tte_data <- function(data, 
                             id_col = "id", 
                             time_col = "event_time", 
                             event_col = "event", 
                             tx_col = "tx_obs", 
                             tx_time_col = "tx_time_obs") {
  
  # Ensure data is a tibble or data.frame
  data <- as.data.frame(data)
  
  # Check if renamed columns already exist
  col_map <- c(
    id = id_col,
    time = time_col,
    event = event_col,
    tx_ever = tx_col,
    t_switch = tx_time_col
  )
  
  for (orig_name in col_map) {
    if (!orig_name %in% names(data)) {
      # Fallback checks for common alternative names
      if (orig_name == "event_time" && "time" %in% names(data)) col_map["time"] <- "time"
      else if (orig_name == "tx_obs" && "tx_ever" %in% names(data)) col_map["tx_ever"] <- "tx_ever"
      else if (orig_name == "tx_time_obs" && "t_switch" %in% names(data)) col_map["t_switch"] <- "t_switch"
      else stop(paste0("Required column '", orig_name, "' not found in dataset."))
    }
  }
  
  # Standardize
  df <- data
  df$id <- df[[col_map["id"]]]
  df$time <- df[[col_map["time"]]]
  df$event <- as.integer(df[[col_map["event"]]])
  df$tx_ever <- as.integer(df[[col_map["tx_ever"]]])
  df$t_switch <- as.numeric(df[[col_map["t_switch"]]])
  
  return(df)
}

# Ensure SL.earth.int exists for 2-way interaction MARS
register_sl_earth_int <- function() {
  if (!exists("SL.earth.int", mode = "function", envir = .GlobalEnv)) {
    assign("SL.earth.int", function(..., degree = 2) {
      SuperLearner::SL.earth(..., degree = degree)
    }, envir = .GlobalEnv)
  }
}

