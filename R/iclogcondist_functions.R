#' Find Unique Rows in a Matrix and Their Weights
#'
#' This function finds the unique rows of a given matrix and calculates the frequency (weight) of each unique row. 
#' It returns both the unique rows and the weights (the number of occurrences of each row).
#'
#' @param X A matrix. The matrix whose unique rows are to be found.
#' @return A list containing two components:
#' \describe{
#'   \item{unique_X}{A matrix of the unique rows from the input matrix.}
#'   \item{weight}{An integer vector containing the frequency (weight) of each unique row.}
#' }
#'
unique_X_weight <-  function(X) {
  X_sorted <- X[do.call(order, as.data.frame(X)), , drop = FALSE]
  unique_indices <- !duplicated(X_sorted)
  row_counts <- diff(c(which(unique_indices), nrow(X_sorted) + 1))
  unique_rows <- X_sorted[unique_indices, , drop = FALSE]
  
  return(list(unique_X = unique_rows, weight = row_counts))
}


#' Prepare Data for Interval-Censored Model
#'
#' This function processes interval-censored data and prepares various components needed for model fitting, including unique time points, censoring intervals, and weights.
#'
#' @param X A matrix or data frame of interval-censored data where each row contains the lower and upper bounds of the interval for each observation.
#' @return A list containing:
#' \describe{
#'   \item{tau}{Unique time points.}
#'   \item{m}{The number of unique time points (excluding infinity if present).}
#'   \item{L_Rc}{Indices of observations where the event is in the intersection of L group and the complement of R group. 
#' The L group consists of samples with left intervals time <= min(all right intervals time). 
#' The R group consists of samples with infinity right interval time.}
#'   \item{Lc_R}{Indices of observations where the event is in the intersection of the complement of L group and R group.}
#'   \item{Lc_Rc}{Indices of observations where the event is in the intersection of the complement of L group and the complement of R group.}
#'   \item{ri}{Indices corresponding to the right bounds of the intervals in \code{tau}.}
#'   \item{li}{Indices corresponding to the left bounds of the intervals in \code{tau}.}
#'   \item{tau_no_Inf}{Unique time points excluding infinity.}
#'   \item{weight}{Weights for each unique interval.}
#'   \item{X}{Processed matrix of interval-censored data with unique rows.}
#' }
data_prep <- function(X) {
  temp <- unique_X_weight(X)
  X <- temp$unique_X
  weight <- temp$weight
  n <- dim(X)[1]
  tau <- unique(sort(c(X[X[, 1] > min(X[, 2]), 1], X[, 2])))
  
  if (max(X) == Inf) {
    m <- length(tau) - 1
    tau_no_Inf <- tau[-length(tau)]
  } else {
    m <- length(tau)
    tau_no_Inf <- tau
  }
  
  L_set <- which(X[, 1] < min(X[, 2]))
  R_set <- which(X[, 2] == Inf)
  all <- 1:n
  L_Rc <- intersect(L_set, setdiff(all, R_set))
  Lc_R <- intersect(setdiff(all, L_set), R_set)
  Lc_Rc <- intersect(setdiff(all, L_set), setdiff(all, R_set))
  ri <- match(X[, 2], tau, nomatch = 0)
  li <- match(X[, 1], tau, nomatch = 0)
  
  list(
    tau = tau, m = m, L_Rc = L_Rc, Lc_R = Lc_R, Lc_Rc = Lc_Rc,
    ri = ri, li = li, tau_no_Inf = tau_no_Inf, weight = weight, X = X
  )
}


#' Initial Values for Estimation under Log-concavity with Interval-Censored Data 
#'
#' This function obtains initial values for the maximum likelihood estimation under log-concavity with 
#' interval-censored data based on the unconstrained maximum likelihood estimate (MLE) 
#' or its least concave majorant (LCM). 
#' Alternatively, the user can provide a numeric vector of initial values.
#'
#' @param X A matrix or data frame of interval-censored data, where each row contains the lower and upper bounds of the interval for each observation.
#' @param initial A character string specifying the method for generating initial values. 
#' The default is \code{"LCM"}, which uses the least concave majorant of the log of the unconstrained MLE. 
#' Other options are \code{"MLE"} for the unconstrained maximum likelihood estimate 
#' or a numeric vector provided by the user.
#' 
#' @return A list containing:
#' \describe{
#'   \item{phi_hat}{The initial values of the \code{phi} parameter based on the specified method.}
#'   \item{phi_hat_MLE}{Initial values based on the unconstrained MLE.}
#'   \item{phi_hat_LCM}{Initial values based on the least concave majorant.}
#' }
#' @import icenReg
#'
initial_values <- function(X, initial = "LCM") {
  fit_MLE <- ic_UMLE(X)
  fit_LCM <- ic_LCM_UMLE(X)
  
  if (is.numeric(initial) && is.vector(initial)) {
    phi_hat <- initial
    if (any(phi_hat > 0) || any(diff(phi_hat) < 0)) {
      stop("Please provide negative and increasing initial values")
    }
  } else if (initial == "LCM") {
    phi_hat <- fit_LCM$est$phi_hat  
  } else if (initial == "MLE") {
    phi_hat <- fit_MLE$est$phi_hat
    if (all(phi_hat == 0)) {
      stop("The unconstrained MLE is all 1. Please provide valid data")
    }
  } else {
    stop("Please provide a valid initialization method")
  }
  
  return(list(
    phi_hat = phi_hat,
    phi_hat_MLE = fit_MLE$est$phi_hat,
    phi_hat_LCM = fit_LCM$est$phi_hat 
  ))
}


#' Compute the Negative Log-Likelihood for the Interval-Censored Model
#'
#' This function computes the negative log-likelihood of an interval-censored model based on the specified parameterization.
#'
#' @param x A numeric vector of parameter estimates (can be in terms of \code{phi}, or \code{F}).
#' @param weight A numeric vector of weights for the observations.
#' @param li A numeric vector of indices corresponding to the left bounds of the intervals in \code{tau_no_Inf}.
#' @param ri A numeric vector of indices corresponding to the right bounds of the intervals in \code{tau_no_Inf}.
#' @param L_Rc Indices of observations where the event is in the intersection of L group and the complement of R group. 
#' The L group consists of samples with left intervals time <= min(all right intervals time). 
#' The R group consists of samples with infinity right interval time.
#' @param Lc_R Indices of observations where the event is in the intersection of the complement of L group and R group.
#' @param Lc_Rc Indices of observations where the event is in the intersection of the complement of L group and the complement of R group.
#' @param type A character string indicating the parameterization of \code{x}. Options are \code{"phi"} (log of F), or \code{"F"}.
#' @param tau_no_Inf A numeric vector of unique time points excluding infinity.
#' @return The negative log-likelihood value.
#' 
neg_log_like <- function(x, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "", tau_no_Inf) {
  if (type == "phi") {
    phi <- x
    F <- exp(phi)
  } else if (type == "F") {
    F <- x
    phi <- log(F)
  }
  
  value <- sum(weight[L_Rc] * phi[ri[L_Rc]]) +
    sum(weight[Lc_R] * log(1 - F[li[Lc_R]])) +
    sum(weight[Lc_Rc] * log(F[ri[Lc_Rc]] - F[li[Lc_Rc]]))
  
  return(-value)
}



#' Compute Directional Derivatives for Active Set Algorithm
#'
#' This function computes the directional derivatives for the active set algorithm used in the estimation 
#' of the distribution function under log-concavity with interval-censored data.
#' The calculation takes advantage of the specific structure of the basis matrix, making it efficient to compute in \code{O(n)} time complexity.
#'
#' @param diff_tau A numeric vector containing the differences between consecutive time points (tau).
#' @param first_order A numeric vector representing the first-order derivatives at each time point.
#' @return A numeric vector of length \code{length(diff_tau) + 1} representing the directional derivatives for the active set algorithm.
#'
find_dir_deriv <- function(diff_tau, first_order) {
  m <- length(diff_tau) + 1
  dir_deriv <- rep(0, m)
  dir_deriv[2:m] <- cumsum(-diff_tau * cumsum(first_order)[-m])
  return(dir_deriv)
}


#' Compute qsi Matrix
#'
#' This function computes the qsi matrix for specified indices and time points.
#'
#' @param is Indices of nodes.
#' @param tau_no_Inf Unique time points excluding infinity.
#' @return qsi matrix.
find_qsi <- function(is, tau_no_Inf) {
  k <- length(is)
  m <- length(tau_no_Inf)
  qsi <- matrix(0, nrow = k - 1, ncol = m)
  for (s in 1:(k - 1)) {
    qsi[s, is[s]:is[s + 1]] <- (tau_no_Inf[is[s + 1]] - tau_no_Inf[is[s]:is[s + 1]]) /
      (tau_no_Inf[is[s + 1]] - tau_no_Inf[is[s]])
  }
  return(qsi)
}



#' Evaluate F(x) for Objects of Class 'iclogcondist'
#'
#' Computes the value of the function \eqn{F(x)} for a given object of class \code{iclogcondist}.
#'
#' @param object An object of class \code{iclogcondist}. Must also belong to one of the subclasses:
#'   \code{"ic_LCMLE"}, \code{"ic_LCM_UMLE"}, or \code{"ic_UMLE"}.
#' @param x A numeric vector of values at which \eqn{F(x)} is evaluated. If not specified,
#'   the \code{tau_no_Inf} attribute of the \code{object} object is used.
#' @param log Logical; if \code{TRUE}, returns the result in log-transformed form.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments (not currently used).
#'
#' @return A numeric vector of values, either \eqn{F(x)} or \eqn{log(F(x))}.
#'
#' @examples
#' # Example usage:
#' data(lgnm)
#' 
#' # Evaluate for LCMLE object
#' fit_LCMLE <- ic_LCMLE(lgnm)
#' get_F_at_x(fit_LCMLE)
#' 
#' # Evaluate for UMLE object
#' fit_UMLE <- ic_UMLE(lgnm)
#' x = seq(0.001, 6, length.out = 1000)
#' get_F_at_x(fit_UMLE, x = x)
#' @export
get_F_at_x.iclogcondist <- function(object, x = NA, log = FALSE, ...) {
  
  if (!inherits(object, "iclogcondist")) stop("Object must be of class 'iclogcondist'")
  if (all(is.na(x))) x = object$est$tau_no_Inf
  
  tau_no_Inf <- object$est$tau_no_Inf
  phi <- object$est$phi_hat
  m <- length(phi)
  k <- length(x)
  
  out <- numeric(k)
  out[x < tau_no_Inf[1]] <- -Inf 
  out[x >= tau_no_Inf[m]] <- phi[m]
  in_range <- which(x >= tau_no_Inf[1] & x < tau_no_Inf[m])
  
  idx <- findInterval(x[in_range], tau_no_Inf)
  tau1 <- tau_no_Inf[idx]
  tau2 <- tau_no_Inf[pmin(idx + 1, m)]
  phi1 <- phi[idx]
  phi2 <- phi[pmin(idx + 1, m)]
  
  if (class(object)[2] %in% c("ic_LCMLE", "ic_LCM_UMLE")) {
    out[in_range] <- phi1 + (phi2 - phi1) * (x[in_range] - tau1) / (tau2 - tau1)
  } else if (class(object)[2] %in% c("ic_UMLE")){
    out[in_range] <- phi1 
  } else {
    stop("Object must be from class 'ic_LCMLE', 'ic_LCM_UMLE', 'ic_UMLE'")
  }
  
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}

#' Generic Function to compute F at X
#'
#' Computes the value of the function \eqn{F(x)} for a given object of class \code{iclogcondist}. 
#' This is a generic function to compute \eqn{F(x)} for object class \code{iclogcondist}. 
#' For usage details, please refer to function \code{get_F_at_x.iclogcondist}
#'
#' @param object An object for which the method is defined.
#' @param ... Additional arguments passed to the method.
#' @return A numeric vector of values, either \eqn{F(x)} or \eqn{log(F(x))}.
#'
#' @examples
#' # Example usage:
#' data(lgnm)
#' 
#' # Evaluate for LCMLE object
#' fit_LCMLE <- ic_LCMLE(lgnm)
#' get_F_at_x(fit_LCMLE)
#' 
#' # Evaluate for UMLE object
#' fit_UMLE <- ic_UMLE(lgnm)
#' x = seq(0.001, 6, length.out = 1000)
#' get_F_at_x(fit_UMLE, x = x)
#' @export 
get_F_at_x <- function(object, ...) {
  UseMethod("get_F_at_x")
}
