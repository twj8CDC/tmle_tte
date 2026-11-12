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
  cov_vars <- setdiff(names(X), c("t", "t_sq", "is_treated", "post_grace"))
  build_cov_term <- function(df_x, vars) {
    if (length(vars) == 0L) return("1")
    terms <- vapply(vars, function(v) {
      if (length(unique(df_x[[v]])) > 5L) {
        paste0("splines::ns(", v, ", df = 3)")
      } else {
        v
      }
    }, character(1))
    base_str <- paste(terms, collapse = " + ")
    if (length(terms) >= 2L) {
      paste0("(", base_str, ")^2")
    } else {
      base_str
    }
  }
  cov_term <- build_cov_term(X, cov_vars)

  if (all(c("t", "is_treated") %in% names(X))) {
    if (length(unique(X$is_treated)) > 1) {
      form_str <- paste0("Y ~ ", cov_term, " + splines::ns(t, df = ", df_ns, ") * is_treated")
    } else {
      form_str <- paste0("Y ~ ", cov_term, " + splines::ns(t, df = ", df_ns, ")")
    }
  } else if ("t" %in% names(X)) {
    form_str <- paste0("Y ~ ", cov_term, " + splines::ns(t, df = ", df_ns, ")")
  } else {
    form_str <- paste0("Y ~ ", cov_term)
  }
  fit.glm <- suppressWarnings(stats::glm(stats::as.formula(form_str), data = X, family = family, weights = obsWeights))
  pred <- suppressWarnings(stats::predict(fit.glm, newdata = newX, type = "response"))
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

#' SuperLearner wrapper for 0th-Order Highly Adaptive Lasso (HAL) 1-Day Cumulative Step Basis
#'
#' @description Constructs a nonparametric 1-day resolution cumulative step-function basis
#'   \eqn{\{I(t \ge k)\}_{k \ge 2}} interacted with treatment status \code{is_treated}, a 3-knot
#'   natural cubic spline expansion on continuous baseline covariates \eqn{W}, pairwise covariate
#'   interactions, and fits an elastic-net/Lasso penalized binomial regression via \code{glmnet}.
#'   When cluster identifiers \code{id} are provided, internal CV folds are grouped by patient
#'   to prevent leakage of person-time records. Discovers arbitrary hazard change-points or kinks
#'   at 1-day granularity without pre-specified time intervals.
#'
#' @param Y Outcome vector.
#' @param X Predictor matrix or data frame.
#' @param newX New predictor data for prediction.
#' @param family Error distribution and link function.
#' @param obsWeights Observation weights.
#' @param id Cluster ID.
#' @param alpha Elastic-net mixing parameter (default 1.0 for pure L1 Lasso).
#' @param nfolds Number of internal CV folds for glmnet (default 3).
#' @param ... Additional arguments.
#' @export
SL.hal.time <- function(Y, X, newX, family, obsWeights, id, alpha = 1.0, nfolds = 3L, ...) {
  has_t <- ("t" %in% names(X)) && (length(unique(X$t)) > 2L)
  if (has_t && requireNamespace("glmnet", quietly = TRUE)) {
    t_vals <- sort(unique(as.numeric(X$t)))[-1]
    step_train <- vapply(t_vals, function(k) as.numeric(X$t >= k), numeric(nrow(X)))
    step_test  <- vapply(t_vals, function(k) as.numeric(newX$t >= k), numeric(nrow(newX)))
    colnames(step_train) <- colnames(step_test) <- paste0("t_ge_", t_vals)
    
    cov_vars <- setdiff(names(X), c("t", "t_sq", "is_treated", "post_grace"))
    spline_fits <- list()
    X_cov_train <- NULL
    X_cov_test  <- NULL
    
    if (length(cov_vars) > 0L) {
      cov_train_list <- list()
      cov_test_list  <- list()
      for (v in cov_vars) {
        vals <- X[[v]]
        if (is.numeric(vals) && length(unique(vals)) > 4L) {
          ns_basis <- splines::ns(vals, df = 3L)
          spline_fits[[v]] <- ns_basis
          cov_train_list[[v]] <- as.matrix(ns_basis)
          colnames(cov_train_list[[v]]) <- paste0(v, "_ns", 1:3)
          cov_test_list[[v]] <- stats::predict(ns_basis, newX[[v]])
          colnames(cov_test_list[[v]]) <- paste0(v, "_ns", 1:3)
        } else {
          cov_train_list[[v]] <- matrix(as.numeric(vals), ncol = 1, dimnames = list(NULL, v))
          cov_test_list[[v]]  <- matrix(as.numeric(newX[[v]]), ncol = 1, dimnames = list(NULL, v))
        }
      }
      
      if (length(cov_vars) >= 2L) {
        combos <- utils::combn(cov_vars, 2, simplify = FALSE)
        for (cb in combos) {
          col_nm <- paste0(cb[1], "_x_", cb[2])
          cov_train_list[[col_nm]] <- matrix(as.numeric(X[[cb[1]]]) * as.numeric(X[[cb[2]]]), ncol = 1, dimnames = list(NULL, col_nm))
          cov_test_list[[col_nm]]  <- matrix(as.numeric(newX[[cb[1]]]) * as.numeric(newX[[cb[2]]]), ncol = 1, dimnames = list(NULL, col_nm))
        }
      }
      
      X_cov_train <- do.call(cbind, cov_train_list)
      X_cov_test  <- do.call(cbind, cov_test_list)
    }
    
    has_trt <- ("is_treated" %in% names(X)) && (length(unique(X$is_treated)) > 1L)
    if (has_trt) {
      trt_tr <- as.numeric(X$is_treated)
      trt_te <- as.numeric(newX$is_treated)
      step_int_tr <- step_train * trt_tr
      step_int_te <- step_test * trt_te
      colnames(step_int_tr) <- colnames(step_int_te) <- paste0("trt_t_ge_", t_vals)
      if (!is.null(X_cov_train)) {
        x_tr <- cbind(X_cov_train, is_treated = trt_tr, step_train, step_int_tr)
        x_te <- cbind(X_cov_test,  is_treated = trt_te, step_test,  step_int_te)
      } else {
        x_tr <- cbind(is_treated = trt_tr, step_train, step_int_tr)
        x_te <- cbind(is_treated = trt_te, step_test,  step_int_te)
      }
    } else {
      if (!is.null(X_cov_train)) {
        x_tr <- cbind(X_cov_train, step_train)
        x_te <- cbind(X_cov_test,  step_test)
      } else {
        x_tr <- step_train
        x_te <- step_test
      }
    }
    
    foldid <- NULL
    if (!missing(id) && !is.null(id) && length(id) == nrow(X)) {
      u_ids <- unique(id)
      if (length(u_ids) >= nfolds) {
        id_fold_map <- stats::setNames(sample(rep(seq_len(nfolds), length.out = length(u_ids))), u_ids)
        foldid <- id_fold_map[as.character(id)]
      }
    }
    
    fit_cv <- tryCatch({
      if (!is.null(foldid)) {
        suppressWarnings(glmnet::cv.glmnet(
          x = x_tr, y = Y, weights = obsWeights,
          family = "binomial", alpha = alpha, foldid = foldid
        ))
      } else {
        suppressWarnings(glmnet::cv.glmnet(
          x = x_tr, y = Y, weights = obsWeights,
          family = "binomial", alpha = alpha, nfolds = nfolds
        ))
      }
    }, error = function(e) NULL)
    
    if (!is.null(fit_cv)) {
      pred <- as.numeric(stats::predict(fit_cv, newx = x_te, s = "lambda.min", type = "response"))
      fit <- list(object = fit_cv, t_vals = t_vals, cov_vars = cov_vars,
                  spline_fits = spline_fits, has_trt = has_trt, has_t = TRUE)
      class(fit) <- "SL.hal.time"
      return(list(pred = pred, fit = fit))
    }
  }
  
  # Fallback to standard GLM if no time column or glmnet fails
  fit.glm <- suppressWarnings(stats::glm(Y ~ ., data = X, family = family, weights = obsWeights))
  pred <- suppressWarnings(stats::predict(fit.glm, newdata = newX, type = "response"))
  fit <- list(object = fit.glm, has_t = FALSE)
  class(fit) <- "SL.hal.time"
  return(list(pred = pred, fit = fit))
}

#' Predict method for SL.hal.time objects
#' @param object An object of class \code{"SL.hal.time"}.
#' @param newdata New predictor data.
#' @param ... Additional arguments.
#' @export
predict.SL.hal.time <- function(object, newdata, ...) {
  if (isTRUE(object$has_t)) {
    t_vals <- object$t_vals
    step_te <- vapply(t_vals, function(k) as.numeric(newdata$t >= k), numeric(nrow(newdata)))
    colnames(step_te) <- paste0("t_ge_", t_vals)
    
    X_cov_te <- NULL
    if (length(object$cov_vars) > 0L) {
      cov_test_list <- list()
      for (v in object$cov_vars) {
        if (v %in% names(object$spline_fits)) {
          cov_test_list[[v]] <- stats::predict(object$spline_fits[[v]], newdata[[v]])
          colnames(cov_test_list[[v]]) <- paste0(v, "_ns", 1:3)
        } else {
          cov_test_list[[v]] <- matrix(as.numeric(newdata[[v]]), ncol = 1, dimnames = list(NULL, v))
        }
      }
      if (length(object$cov_vars) >= 2L) {
        combos <- utils::combn(object$cov_vars, 2, simplify = FALSE)
        for (cb in combos) {
          col_nm <- paste0(cb[1], "_x_", cb[2])
          cov_test_list[[col_nm]] <- matrix(as.numeric(newdata[[cb[1]]]) * as.numeric(newdata[[cb[2]]]), ncol = 1, dimnames = list(NULL, col_nm))
        }
      }
      X_cov_te <- do.call(cbind, cov_test_list)
    }
    
    if (isTRUE(object$has_trt)) {
      trt_te <- as.numeric(newdata$is_treated)
      step_int_te <- step_te * trt_te
      colnames(step_int_te) <- paste0("trt_t_ge_", t_vals)
      if (!is.null(X_cov_te)) {
        x_te <- cbind(X_cov_te, is_treated = trt_te, step_te, step_int_te)
      } else {
        x_te <- cbind(is_treated = trt_te, step_te, step_int_te)
      }
    } else {
      if (!is.null(X_cov_te)) {
        x_te <- cbind(X_cov_te, step_te)
      } else {
        x_te <- step_te
      }
    }
    return(as.numeric(stats::predict(object$object, newx = x_te, s = "lambda.min", type = "response")))
  } else {
    return(suppressWarnings(stats::predict(object$object, newdata = newdata, type = "response")))
  }
}

# Ensure SuperLearner wrappers exist in global environment for SuperLearner workers
register_sl_earth_int <- function() {
  if (!exists("SL.earth.int", mode = "function", envir = .GlobalEnv)) {
    assign("SL.earth.int", SL.earth.int, envir = .GlobalEnv)
  }
  if (!exists("SL.glm.time.int", mode = "function", envir = .GlobalEnv)) {
    assign("SL.glm.time.int", SL.glm.time.int, envir = .GlobalEnv)
  }
  if (!exists("SL.hal.time", mode = "function", envir = .GlobalEnv)) {
    assign("SL.hal.time", SL.hal.time, envir = .GlobalEnv)
  }
  if (!exists("predict.SL.hal.time", mode = "function", envir = .GlobalEnv)) {
    assign("predict.SL.hal.time", predict.SL.hal.time, envir = .GlobalEnv)
  }
  if (!exists("SL.bart", mode = "function", envir = .GlobalEnv)) {
    assign("SL.bart", SL.bart, envir = .GlobalEnv)
  }
  if (!exists("predict.SL.bart", mode = "function", envir = .GlobalEnv)) {
    assign("predict.SL.bart", predict.SL.bart, envir = .GlobalEnv)
  }
}

