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

#' SuperLearner wrapper for Earth (MARS) with 2-way interactions
#'
#' @param ... Arguments passed to SuperLearner::SL.earth
#' @param degree Interaction degree (default 2)
#' @export
SL.earth.int <- function(..., degree = 2) {
  SuperLearner::SL.earth(..., degree = degree)
}

#' SuperLearner wrapper for GLM with Natural Spline Time-by-Treatment Interactions
#'
#' @description Fits a discrete-time hazard model with a flexible natural cubic spline
#'   basis for time interacted with treatment indicator (\code{is_treated}). This allows the
#'   hazard trajectory and treatment effect to vary smoothly and nonparametrically across follow-up
#'   without requiring a priori prespecification of waning periods or change-points.
#'
#' @param Y Outcome vector.
#' @param X Predictor matrix or data frame.
#' @param newX New predictor data for prediction.
#' @param family Error distribution and link function.
#' @param obsWeights Observation weights.
#' @param id Cluster ID.
#' @param df_ns Degrees of freedom for natural splines (default 4).
#' @param ... Additional arguments.
#' @export
SL.glm.time.int <- function(Y, X, newX, family, obsWeights, id, df_ns = 4, ...) {
  if (all(c("t", "is_treated") %in% names(X))) {
    if (length(unique(X$is_treated)) > 1) {
      cov_vars <- setdiff(names(X), c("t", "t_sq", "is_treated", "post_grace"))
      cov_term <- if (length(cov_vars) > 0) paste(cov_vars, collapse = " + ") else "1"
      form_str <- paste0("Y ~ ", cov_term, " + splines::ns(t, df = ", df_ns, ") * is_treated")
    } else {
      cov_vars <- setdiff(names(X), c("t", "t_sq", "is_treated", "post_grace"))
      cov_term <- if (length(cov_vars) > 0) paste(cov_vars, collapse = " + ") else "1"
      form_str <- paste0("Y ~ ", cov_term, " + splines::ns(t, df = ", df_ns, ")")
    }
  } else if ("t" %in% names(X)) {
    cov_vars <- setdiff(names(X), c("t", "t_sq", "post_grace"))
    cov_term <- if (length(cov_vars) > 0) paste(cov_vars, collapse = " + ") else "1"
    form_str <- paste0("Y ~ ", cov_term, " + splines::ns(t, df = ", df_ns, ")")
  } else {
    form_str <- "Y ~ ."
  }
  fit.glm <- stats::glm(stats::as.formula(form_str), data = X, family = family, weights = obsWeights)
  pred <- stats::predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  return(list(pred = pred, fit = fit))
}

#' SuperLearner wrapper for BART (Bayesian Additive Regression Trees)
#'
#' @description Fits a Bayesian probit sum-of-trees model using \code{BART::pbart}
#'   for nonparametric discovery of time-varying hazards and compliance probabilities.
#'
#' @param Y Outcome vector.
#' @param X Predictor matrix or data frame.
#' @param newX New predictor data for prediction.
#' @param family Error distribution and link function.
#' @param obsWeights Observation weights.
#' @param id Cluster ID.
#' @param ntree Number of trees (default 50).
#' @param ndpost Number of posterior draws (default 100).
#' @param nskip Number of burn-in draws (default 50).
#' @param ... Additional arguments.
#' @export
SL.bart <- function(Y, X, newX, family, obsWeights, id, 
                    ntree = 50L, ndpost = 100L, nskip = 50L, ...) {
  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package 'BART' is required for SL.bart.")
  }
  fit_bart <- BART::pbart(
    x.train = X, 
    y.train = Y, 
    x.test = newX, 
    ntree = ntree, 
    ndpost = ndpost, 
    nskip = nskip, 
    printevery = 1000000L
  )
  pred <- colMeans(fit_bart$prob.test)
  fit <- list(object = fit_bart, train_cols = colnames(X), ndpost = ndpost, nskip = nskip, ntree = ntree)
  class(fit) <- "SL.bart"
  return(list(pred = pred, fit = fit))
}

#' Predict method for SL.bart objects
#' @param object An object of class \code{"SL.bart"}.
#' @param newdata New predictor data.
#' @param ... Additional arguments.
#' @export
predict.SL.bart <- function(object, newdata, ...) {
  cols <- object$train_cols
  if (!is.null(cols) && all(cols %in% colnames(newdata))) {
    newdata_sub <- as.data.frame(newdata)[, cols, drop = FALSE]
  } else {
    newdata_sub <- as.data.frame(newdata)
  }
  pred_res <- stats::predict(object$object, newdata = newdata_sub)
  return(colMeans(pred_res$prob.test))
}

# Ensure SuperLearner wrappers exist in global environment for SuperLearner workers
register_sl_earth_int <- function() {
  if (!exists("SL.earth.int", mode = "function", envir = .GlobalEnv)) {
    assign("SL.earth.int", SL.earth.int, envir = .GlobalEnv)
  }
  if (!exists("SL.glm.time.int", mode = "function", envir = .GlobalEnv)) {
    assign("SL.glm.time.int", SL.glm.time.int, envir = .GlobalEnv)
  }
  if (!exists("SL.bart", mode = "function", envir = .GlobalEnv)) {
    assign("SL.bart", SL.bart, envir = .GlobalEnv)
  }
  if (!exists("predict.SL.bart", mode = "function", envir = .GlobalEnv)) {
    assign("predict.SL.bart", predict.SL.bart, envir = .GlobalEnv)
  }
}

