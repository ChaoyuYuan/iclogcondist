#' Compute the log-concave MLE for Interval-censored Data using an Active Set Algorithm
#'
#' This function computes the log-concave MLE of the cumulative distribution function for interval-censored data under log-concavity
#' on the underlying distribution function based on an active set algorithm. The active set algorithm adjusts the knots set based on
#' certain directional derivatives.
#'
#' @param X A matrix with two columns, where each row represents an interval (L, R] for interval-censored data.
#'          \code{L} and \code{R} are the left and right endpoints, with \code{R = Inf} indicating right-censoring.
#' @param initial A character string specifying the method of obtaining an initial value ("LCM" or "MLE") for the estimation process. Default is \code{"LCM"}.
#' @param print Logical. If \code{TRUE}, prints the iterative process details. Default is \code{FALSE}.
#' @param max_iter An integer specifying the maximum number of iterations for the algorithm. Default is 500.
#' @param tol_conv A numeric tolerance level for convergence based on the directional derivatives. Default is \code{1e-7}.
#' @param tol_conv_like A numeric tolerance level for convergence based on log-likelihood difference. Default is \code{1e-10}.
#' @param tol_K A numeric tolerance for checking if v^T phi is in K (constraint set) in active constraint set \code{A}. Default is \code{1e-5}.
#' @return A list with the following components:
#' \item{est}{A list containing \code{tau_no_Inf} (unique finite \code{tau} values), \code{phi_hat} (estimate of \code{log F}), and \code{F_hat} (estimate of \code{F}).}
#' \item{knot_info}{A list with \code{knot_index} (indices of active knots in \code{tau_no_Inf}), \code{tau_on_knot} (tau values at active knots), 
#'                 \code{F_on_knot} (MLE cumulative distribution function values at active knots), and \code{phi_on_knot} (logarithmic estimates of \code{F} at active knots).}
#' \item{neg_log_likelihood}{Vector of negative log-likelihood values for each iteration of the algorithm.}
#' \item{dir_derivs}{Vector of directional derivatives for each iteration.}
#' \item{iter_no}{Integer representing the total number of iterations.}
#' \item{weight}{Vector of weights corresponding to each interval in the data.}
#' \item{X}{The original interval-censored data matrix input.}
#' @export
#' @useDynLib iclogcondist, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' # Example usage:
#' data(lgnm)
#' result <- ic_LCMLE(X = lgnm, initial = "LCM", print = TRUE, max_iter = 500)
#' print(result$est)
ic_LCMLE <- function(X,
                      initial = "LCM",
                      print = FALSE,
                      max_iter = 500,
                      tol_conv = 1e-7,
                      tol_conv_like = 1e-10,
                      tol_K = 1e-5) {
  # Initialize values
  temp <- initial_values(X, initial)
  phi_initial <- temp$phi_hat
  
  temp <- data_prep(X)
  tau_no_Inf <- temp$tau_no_Inf
  m <- temp$m
  L_Rc <- temp$L_Rc
  Lc_R <- temp$Lc_R
  Lc_Rc <- temp$Lc_Rc
  ri <- temp$ri
  li <- temp$li
  weight <- temp$weight
  diff_tau <- diff(tau_no_Inf)

  neg_log_like_initial <- neg_log_like(phi_initial, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf)
  
  # Initialize tracking variables
  dir_deriv_save <- rep(0, max_iter)
  neg_likeli_save <- rep(0, max_iter)
  diff_log_like <- 100
  
  # Step 1: Update initial value
  if (print) cat("Step 1 (update initial value) Start\n")
  v_phi <- compute_constraint(diff_tau, phi_initial)
  A <- which(v_phi >= -tol_K) + 1
  is <- setdiff(1:m, A)
  
  temp <- convert_phi_cand_in_K(phi_initial, v_phi, diff_tau, tol_K, print, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight, m)
  phi_hat <- temp$phi_cand
  v_phi <- temp$v_phi_cand
  A <- which(v_phi >= -tol_K) + 1
  is <- setdiff(1:m, A)
  if (print) cat("Step 1 (update initial value) Done\n=======================\n")
  
  # Step 2: Altering the set of active constraints
  if (print) cat("Step 2 (Altering the set of active constraints) Start\n")
  first_order <- compute_first_second_diag_order_phi(weight, phi_hat, li, ri, L_Rc, Lc_R, Lc_Rc)$first_order
  dir_deriv <- find_dir_deriv(diff_tau, first_order)
  
  iter <- 0
  number_adjust_step2 = 0
  while (max(dir_deriv[A]) > tol_conv && abs(diff_log_like) > tol_conv_like && iter < max_iter) {
    iter <- iter + 1
    is <- sort(union(is, which.max(dir_deriv)))
    
    temp <- convert_phi_cand_in_K(phi_hat, v_phi, diff_tau, tol_K, FALSE, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight, m)
    phi_hat <- temp$phi_cand
    v_phi <- temp$v_phi_cand
    A <- which(v_phi >= -tol_K) + 1
    is <- setdiff(1:m, A)
    
    first_order <- compute_first_second_diag_order_phi(weight, phi_hat, li, ri, L_Rc, Lc_R, Lc_Rc)$first_order
    dir_deriv <- find_dir_deriv(diff_tau, first_order)
    
    dir_deriv_save[iter] <- max(dir_deriv[A])
    neg_likeli_save[iter] <- neg_log_like(phi_hat, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf)
    if (iter > 1) diff_log_like <- neg_likeli_save[iter - 1] - neg_likeli_save[iter]
    
    if (print) {
      cat("Iter:", iter, "\nKnots:", paste(is, collapse = ", "), "\nDD:", dir_deriv_save[iter], "\nneg log-like:", neg_likeli_save[iter], "\n=======================\n")
    }
  }
  
  if (iter == max_iter) print("Max iterations reached; check if log-likelihood has stabilized.")
  
  est <- list(tau_no_Inf = tau_no_Inf, phi_hat = phi_hat, F_hat = exp(phi_hat))
  knot_info <- list(
    knot_index = is,
    tau_on_knot = tau_no_Inf[is],
    F_on_knot = exp(phi_hat[is]),
    phi_on_knot = phi_hat[is]
  )
  output <- list(
    est = est, 
    knot_info = knot_info, 
    neg_log_likelihood_list = neg_likeli_save[1:iter],
    neg_log_likelihood = neg_likeli_save[iter], 
    dir_derivs = dir_deriv_save[1:iter], 
    iter_no = iter, weight = weight, X = X)
  class(output) <- c("iclogcondist", "ic_LCMLE")
  return(output)
}


# Get optimized phi with constraint A(phi), where A + is = 1:m
iclogcondist_subset_algorithm <- function(phi, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight) {
  qsi <- find_qsi(is, tau_no_Inf)
  phi_tilde <- phi[is]
  
  refit <- icm_subset_cpp(phi_tilde, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight)
  phi_new <- find_re_phi_cpp(refit$phi_tilde_hat, is, qsi)
  return(phi_new)
}


# With phi_cur, generate A = A(phi_cur) and get phi_cand = tilde(A) and then make phi_cand belong to domain K
convert_phi_cand_in_K <- function(phi_cur, v_phi, diff_tau, tol_K, print, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight, m) {
  phi_cand <- iclogcondist_subset_algorithm(phi_cur, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight)
  v_phi_cand <- compute_constraint(diff_tau, phi_cand)
  iter_1 <- 0
  
  while (max(v_phi_cand) > tol_K) {
    iter_1 <- iter_1 + 1
    if (print) cat("Iter", iter_1, "\n")
    
    index <- which(v_phi_cand > tol_K)
    t_value <- min(-v_phi[index] / (v_phi_cand[index] - v_phi[index]))
    temp_t_value = -v_phi[index] /(v_phi_cand[index] - v_phi[index])
    temp_t_value_01 = temp_t_value[temp_t_value < 1 & temp_t_value > 0]
    if (length(temp_t_value_01) > 0) {
      t_value = min(temp_t_value_01)  
    } else {
      break
    }
    phi_cur <- (1 - t_value) * phi_cur + t_value * phi_cand
    
    v_phi <- compute_constraint(diff_tau, phi_cur)
    A <- which(v_phi >= -tol_K) + 1
    is <- setdiff(1:m, A)
    phi_cand <- iclogcondist_subset_algorithm(phi_cur, is, tau_no_Inf, L_Rc, Lc_R, Lc_Rc, ri, li, weight)
    
    v_phi_cand <- compute_constraint(diff_tau, phi_cand)
    if (t_value < 1e-10) {
      break
    }
  }
  
  return(list(phi_cand = phi_cand, v_phi_cand = v_phi_cand))
}


#' Iterative Convex Minorant (ICM) Subset Algorithm
#'
#' This function implements the ICM algorithm for solving the sub-problem in the active set algorithm.
#' This is a support of the active set algorithm, computing the optimal values \code{phi_tilde} with reduced number of knots in the sub-problem.
#' It uses backtracking to ensure convergence (Jongbloed, 1998). 
#'
#' @param phi_tilde_initial A numeric vector representing the initial values of the reduced variables \code{phi_tilde}.
#' @param is A numeric vector indicating the nodes with unequal left-hand slope and right-hand slope.
#' @param tau_no_Inf A numeric vector containing the unique time points, excluding infinity.
#' @param L_Rc Indices of observations where the event is in the intersection of L group and the complement of R group. 
#' The L group consists of samples with left intervals time <= min(all right intervals time). 
#' The R group consists of samples with infinity right interval time.
#' @param Lc_R Indices of observations where the event is in the intersection of the complement of L group and R group.
#' @param Lc_Rc Indices of observations where the event is in the intersection of the complement of L group and the complement of R group.
#' @param ri A numeric vector of indices corresponding to the right bounds of the intervals in \code{tau_no_Inf}.
#' @param li A numeric vector of indices corresponding to the left bounds of the intervals in \code{tau_no_Inf}.
#' @param weight A numeric vector representing the weights for each observation.
#' @param tol A numeric value specifying the tolerance for convergence. Default is \code{1e-10}.
#' @param max_iter An integer specifying the maximum number of iterations. Default is \code{500}.
#' @return A list containing:
#' \describe{
#'   \item{phi_tilde_hat}{The estimated values of the reduced variable \code{phi_tilde} at the end of the ICM iterations.}
#' }
#' @import monotone
#'
#' 
#' @references 
#' Jongbloed, G.: The iterative convex minorant algorithm for nonparametric estimation. J. Comput. Gr. Stat. 7(3), 310–321 (1998)
#' 
icm_subset_cpp = function(phi_tilde_initial,
                          is,
                          tau_no_Inf,
                          L_Rc,
                          Lc_R,
                          Lc_Rc,
                          ri,
                          li,
                          weight,
                          tol = 1e-10,
                          max_iter = 500) {
  m = length(tau_no_Inf)
  tol = 1e-7
  
  
  conv_criterion = Inf
  iter = 0
  likelihood_save = rep(0, max_iter)
  delta = 0.1
  line_search = 0.25
  
  phi_tilde_cur = phi_tilde_initial
  qsi = find_qsi(is, tau_no_Inf)
  
  
  while (conv_criterion  > tol & iter < max_iter) {
    iter = iter + 1
    phi_cur = find_re_phi_cpp(phi_tilde_cur, is, qsi)
    
    ## Step 1: Isotonic Regression
    # compute 1st and 2nd derivatives w.r.t. phi, then derivative w.r.t phi_tilde (reduced variables)
    temp = compute_first_second_diag_order_phi(weight, phi_cur, li, ri, L_Rc, Lc_R, Lc_Rc)
    phi_first_order = temp$first_order
    phi_second_order = temp$second_order
    
    temp = compute_1st_2nd_order_tilde_cpp(phi_first_order, phi_second_order, is, qsi)
    phi_tilde_first_order = temp$first_order_tilde
    phi_tilde_second_order = temp$second_order_tilde
    phi_tilde_second_order[phi_tilde_second_order == 0] = -delta
    
    # change to convex function
    phi_tilde_first_order = -phi_tilde_first_order
    phi_tilde_second_order = -phi_tilde_second_order
    
    # isotonic regression
    phi_tilde_cand = pmin(
      monotone(
        phi_tilde_cur - 1 / phi_tilde_second_order * phi_tilde_first_order,
        w = phi_tilde_second_order
      ),
      0
    )
    
    
    ## Step 2: Damping
    
    l_cur = neg_log_like(phi_cur, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf)
    phi_cand = find_re_phi_cpp(phi_tilde_cand, is, qsi)
    l_new = neg_log_like(phi_cand, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf)
    
    if (l_new < l_cur + (1 - line_search) * sum(phi_tilde_first_order *
                                                (phi_tilde_cand - phi_tilde_cur))) {
      phi_tilde_cand = phi_tilde_cand
    } else {
      lambda = 1
      s = 1 / 2
      phi_tilde_s = phi_tilde_cand
      
      phi_s = find_re_phi_cpp(phi_tilde_s, is, qsi)
      l_s = neg_log_like(phi_s, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf)
      cond1 = l_cur + (1 - line_search) * sum(phi_tilde_first_order * (phi_tilde_s - phi_tilde_cur))
      cond2 = l_cur + line_search * sum(phi_tilde_first_order * (phi_tilde_s - phi_tilde_cur))
      
      
      while (l_s < cond1 | l_s > cond2) {
        if (l_s < cond1) {
          lambda = lambda + s
        }
        if (l_s > cond2) {
          lambda = lambda - s
        }
        #
        phi_tilde_s = phi_tilde_cur + lambda * (phi_tilde_cand - phi_tilde_cur)
        s = s / 2
        
        phi_s = find_re_phi_cpp(phi_tilde_s, is, qsi)
        l_s = neg_log_like(phi_s,
                           weight,
                           li,
                           ri,
                           L_Rc,
                           Lc_R,
                           Lc_Rc,
                           type = "phi",
                           tau_no_Inf)
        cond1 = l_cur + (1 - line_search) * sum(phi_tilde_first_order *
                                                  (phi_tilde_s - phi_tilde_cur))
        cond2 = l_cur + line_search * sum(phi_tilde_first_order * (phi_tilde_s - phi_tilde_cur))
      }
      phi_tilde_cand = phi_tilde_s
    }
    
    phi_cand = find_re_phi_cpp(phi_tilde_cand, is, qsi)
    l_new = neg_log_like(phi_cand, weight, li, ri, L_Rc, Lc_R, Lc_Rc, type = "phi", tau_no_Inf)
    conv_criterion = l_cur - l_new
    
    phi_tilde_cur = phi_tilde_cand
    
    conv_set = c(
      sum(phi_tilde_cur * phi_tilde_first_order),
      sum(phi_tilde_first_order),
      min(
        sum(phi_tilde_first_order) - cumsum(phi_tilde_first_order)
      )
    )
    
    likelihood_save[iter] = l_new
  }
  
  
  return(list(phi_tilde_hat = phi_tilde_cur))
}