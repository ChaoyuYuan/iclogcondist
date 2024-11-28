#' Compute Unconstrained Maximum Likelihood Estimate for Interval-Censored Data
#'
#' This function computes the unconstrained maximum likelihood estimate (UMLE) for interval-censored data. 
#' It utilizes the non-parametric MLE from the \code{ic_np} function in the \code{icenReg} package as a starting point 
#' and prepares key components such as cumulative probabilities, log-transformed values, and knot information.
#'
#' @details
#' The \code{ic_np} function from the \code{icenReg} package is used to compute the non-parametric MLE for 
#' interval-censored data. This provides initial estimates of probabilities (\code{p_hat}) and jump points 
#' (\code{knot}) in the cumulative distribution function. These are then processed to compute the 
#' cumulative probabilities (\code{F_hat}) and log-transformed values (\code{phi_hat}) at unique time points.
#' 
#' @param X A matrix with two columns, where each row represents an interval (L, R] for interval-censored data.
#'          \code{L} and \code{R} are left and right endpoints, respectively, with \code{R = Inf} indicating right-censoring.
#' @return A list containing:
#' \item{est}{A list with \code{tau_no_Inf} (finite values of \code{tau}), \code{phi_hat} (log of \code{F_hat} values), 
#'           and \code{F_hat} (MLE cumulative distribution function values).}
#' \item{knot_info}{A list with \code{knot_index} (indices of knots in \code{tau_no_Inf}), \code{tau_on_knot} (values of \code{tau} at knots), 
#'                 \code{F_on_knot} (MLE at knots), and \code{phi_on_knot} (log of \code{F} at knots).}
#' \item{neg_log_likelihood}{The negative log-likelihood of the MLE fit.}
#' \item{weight}{Vector of weights corresponding to each interval in the data.}
#' \item{X}{The original interval-censored data matrix input.}
#' @importFrom icenReg ic_np
#' @examples
#' data(lgnm)
#' result <- ic_UMLE(X = lgnm)
#' @export
#' 
#' @references 
#' Anderson-Bergman, C. (2016) An efficient implementation of the EMICM algorithm for the interval censored NPMLE
#' \emph{Journal of Computational and Graphical Statistics}.
#' 
ic_UMLE <- function(X) {
  
  ## Step 1: Compute non-parametric MLE for interval-censored data and data prepare
  fit <- ic_np(X)
  temp <- data_prep(X)
  tau_no_Inf <- temp$tau_no_Inf
  m <- temp$m
  L_Rc <- temp$L_Rc
  Lc_R <- temp$Lc_R
  Lc_Rc <- temp$Lc_Rc
  ri <- temp$ri
  li <- temp$li
  weight <- temp$weight
  
  # Define F(t) as constant in intervals, with jumps at right interval points (F(R) - F(L) likelihood)
  index <- fit$p_hat > 0
  knot <- fit$T_bull_Intervals[2, index]
  p_hat <- fit$p_hat[index]
  
  # Initialize F_hat for cumulative probabilities over tau_no_Inf
  F_hat <- rep(0, length(tau_no_Inf))
  knot_index = which(tau_no_Inf %in% knot)
  if (tail(knot, 1) == Inf) {
    F_hat[knot_index] <- p_hat[-length(p_hat)]
  } else {
    F_hat[knot_index] <- p_hat
  }
  F_hat <- cumsum(F_hat)
  phi_hat <- log(F_hat)
  knot_index = union(knot_index, m)# Include the last point if the tail is flat
  
  ## Step 2: Package outputs for MLE
  est <- list(
    tau_no_Inf = tau_no_Inf,
    phi_hat = phi_hat,
    F_hat = exp(phi_hat)
  )
  knot_info <- list(
    knot_index = knot_index,
    tau_on_knot = tau_no_Inf[knot_index],
    F_on_knot = F_hat[knot_index],
    phi_on_knot = log(F_hat[knot_index])
  )
  output <- list(
    est = est,
    knot_info = knot_info,
    neg_log_likelihood = neg_log_like(phi_hat, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf), 
    weight = weight, X = X
  )
  
  class(output) <- c("iclogcondist", "ic_UMLE")
  return(output)
}

#' Compute Least Concave Majorant (LCM) of the log of the Unconstrained MLE for Interval-Censored Data
#'
#' This function computes the Least Concave Majorant (LCM) of the log of the unconstrained MLE for interval-censored data. 
#'
#' @param X A matrix with two columns, where each row represents an interval (L, R] for interval-censored data.
#'          \code{L} and \code{R} are the left and right endpoints, respectively, with \code{R = Inf} indicating right-censoring.
#' @return A list containing:
#' \item{est}{A list with \code{tau_no_Inf} (finite values of \code{tau}), \code{phi_hat} (LCM of log of \code{F_hat}), 
#'           and \code{F_hat} (exp of \code{phi_hat}).}
#' \item{knot_info}{A list with \code{knot_index} (indices of knots in \code{tau_no_Inf}), 
#'                  \code{tau_on_knot} (values of \code{tau} at knots), 
#'                 \code{F_on_knot} (\code{F_hat} at knots), 
#'                 and \code{phi_on_knot} (\code{phi_hat} at knots).}
#' \item{neg_log_likelihood}{The negative log-likelihood of the LCM fit.}
#' \item{weight}{Vector of weights corresponding to each interval in the data.}
#' \item{X}{The original interval-censored data matrix input.}
#' @importFrom fdrtool gcmlcm
#' @examples
#' data(lgnm)
#' result <- ic_LCM_UMLE(X = lgnm)
#' @export
ic_LCM_UMLE <- function(X) {
  
  ## Step 1: Compute MLE from interval-censored data
  fit <- ic_UMLE(X)
  tau_no_Inf <- fit$est$tau_no_Inf
  phi_hat <- fit$est$phi_hat
  
  temp <- data_prep(X)
  m <- temp$m
  L_Rc <- temp$L_Rc
  Lc_R <- temp$Lc_R
  Lc_Rc <- temp$Lc_Rc
  ri <- temp$ri
  li <- temp$li
  weight <- temp$weight
  
  ## Step 2: Compute Least Concave Majorant (LCM) of phi_hat
  phi_hat_LCM_knot <- gcmlcm(tau_no_Inf, phi_hat, type = "lcm")
  tau_tilde <- phi_hat_LCM_knot$x.knots
  phi_tilde <- phi_hat_LCM_knot$y.knots
  
  # Include the last point if the tail is flat
  is <- union(match(tau_tilde, tau_no_Inf), length(tau_no_Inf))
  qsi <- find_qsi(is, tau_no_Inf)
  phi_hat_LCM <- find_re_phi_cpp(phi_tilde, is, qsi)
  
  ## Step 3: Package outputs for LCM
  est <- list(
    tau_no_Inf = tau_no_Inf,
    phi_hat = phi_hat_LCM,
    F_hat = exp(phi_hat_LCM)
  )
  knot_info <- list(
    knot_index = is,
    tau_on_knot = tau_no_Inf[is],
    F_on_knot = exp(phi_tilde),
    phi_on_knot = phi_tilde
  )
  output <- list(
    est = est,
    knot_info = knot_info,
    is = is,
    neg_log_likelihood = neg_log_like(phi_hat_LCM, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf), 
    weight = weight, X = X
  )
  
  class(output) <- c("iclogcondist", "ic_LCM_UMLE")
  return(output)
}







