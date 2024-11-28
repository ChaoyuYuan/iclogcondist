#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/2
//

// for active set algo, to check whether phi satisfies the constraint set
// [[Rcpp::export]]
NumericVector compute_constraint(const NumericVector& delta_tau, const NumericVector& phi) {
// NumericVector compute_constraint(NumericVector delta_tau, NumericVector phi) {
  int m = delta_tau.size() + 1; // delta_tau has m-1 elements, so m = delta_tau.size() + 1
  NumericVector result(m - 1);  // Result vector for indices i = 2, ..., m
  
  for (int i = 1; i < m - 1; ++i) {
    double term1 = 1.0 / delta_tau[i-1] * phi[i-1];
    double term2 = (1.0 / delta_tau[i] + 1.0 / delta_tau[i-1]) * phi[i];
    double term3 = 1.0 / delta_tau[i] * phi[i+1];
    
    result[i-1] = term1 - term2 + term3;
  }
  result[m-2] = 1000*(phi[m-2] - phi[m-1]);
  
  return result;
}




// [[Rcpp::export]]
NumericVector find_re_phi_cpp(NumericVector tilde_phi, IntegerVector is, NumericMatrix qsi) {
  int k = is.size();
  int m = qsi.ncol();
  NumericVector re_phi(m);
  
  // linear interpolation to get back the values between the knots
  for (int s = 0; s < (k - 1); ++s) {
    for (int j = is[s] - 1; j < is[s + 1]; ++j) {
      re_phi[j] = qsi(s, j) * tilde_phi[s] + (1 - qsi(s, j)) * tilde_phi[s + 1];
    }
  }
  
  // handle the case when tau_m is not a knot, then phi is flat from tau_{is[k]+1} to tau_m
  if (is[k-1] < m) {
    for (int j = is[k-1]; j < m; ++j) {
      re_phi[j] = tilde_phi[k-1];
    }
  }
  
  return re_phi;
}

// [[Rcpp::export]]
Rcpp::List compute_first_second_diag_order_phi(
    Rcpp::NumericVector weight,
    Rcpp::NumericVector phi,
    Rcpp::IntegerVector li,
    Rcpp::IntegerVector ri,
    Rcpp::IntegerVector L_Rc,
    Rcpp::IntegerVector Lc_R,
    Rcpp::IntegerVector Lc_Rc) {
  
  // This is to calculate the first order and diagonal second order
  // For icm_subset in main algorithm
  
  int m = phi.size();
  Rcpp::NumericVector first_order(m);
  Rcpp::NumericVector second_order(m);
  
  // Create local copies of the input index vectors and adjust them to be zero-based
  IntegerVector L_Rc_copy = clone(L_Rc) - 1;
  IntegerVector Lc_R_copy = clone(Lc_R) - 1;
  IntegerVector Lc_Rc_copy = clone(Lc_Rc) - 1;
  IntegerVector ri_copy = clone(ri) - 1;
  IntegerVector li_copy = clone(li) - 1;
  
  NumericVector F = exp(phi);
  NumericVector Fn = exp(-phi);
  
  // First loop: L_Rc
  for (int i : L_Rc_copy) {
    first_order[ri_copy[i]] += weight[i];
    second_order[ri_copy[i]] -= 0;
  }
  
  // Second loop: Lc_R
  for (int i : Lc_R_copy) {
    double denom = Fn[li_copy[i]] - 1;
    first_order[li_copy[i]] -= weight[i] / denom;
    second_order[li_copy[i]] -= (weight[i] / (denom * denom)) * Fn[li_copy[i]];
  }
  
  // Third loop: Lc_Rc, r[i] = j
  for (int i : Lc_Rc_copy) {
    double denom = 1 - F[li_copy[i]] * Fn[ri_copy[i]];
    first_order[ri_copy[i]] += weight[i] / denom;
    second_order[ri_copy[i]] -= (weight[i] / (denom * denom)) * F[li_copy[i]] * Fn[ri_copy[i]];
  }
  
  // Fourth loop: Lc_Rc, l[i] = j
  for (int i : Lc_Rc_copy) {
    double denom = F[ri_copy[i]] * Fn[li_copy[i]] - 1;
    first_order[li_copy[i]] -= weight[i] / denom;
    second_order[li_copy[i]] -= (weight[i] / (denom * denom)) * F[ri_copy[i]] * Fn[li_copy[i]];
  }
  
  return Rcpp::List::create(
    Rcpp::Named("first_order") = first_order,
    Rcpp::Named("second_order") = second_order
  );
}


 

// [[Rcpp::export]]
List compute_1st_2nd_order_tilde_cpp(NumericVector first_order, NumericVector second_order, 
                                     IntegerVector is, NumericMatrix qsi) {
  int k = is.size();
  int m = qsi.ncol();
  
  NumericVector first_order_tilde(k, 0.0);
  NumericVector second_order_tilde(k, 0.0);
  
  // First element for first and second order
  for (int j = is[0] - 1; j < is[1] - 1; ++j) {
    first_order_tilde[0] += first_order[j] * qsi(0, j);
    second_order_tilde[0] += second_order[j] * std::pow(qsi(0, j), 2);
  }
  
  // Loop over s = 2:(k-1)
  for (int s = 1; s < (k - 1); ++s) {
    // First order tilde
    for (int j = is[s-1]; j < is[s] - 1; ++j) {
      first_order_tilde[s] += first_order[j] * (1 - qsi(s-1, j));
      second_order_tilde[s] += second_order[j] * std::pow(1 - qsi(s-1, j), 2);
    }
    for (int j = is[s] - 1; j < is[s+1] - 1; ++j) {
      first_order_tilde[s] += first_order[j] * qsi(s, j);
      second_order_tilde[s] += second_order[j] * std::pow(qsi(s, j), 2);
    }
  }
  
  // Last element for first and second order
  for (int j = is[k-2]; j < is[k-1]; ++j) {
    first_order_tilde[k-1] += first_order[j] * (1 - qsi(k-2, j));
    second_order_tilde[k-1] += second_order[j] * std::pow(1 - qsi(k-2, j), 2);
  }
  
  // handle the case when i(k) < m
  if (is[k-1] < m) {
    for (int j = is[k-1]; j < m; ++j) {
      first_order_tilde[k-1] += first_order[j];
      second_order_tilde[k-1] += second_order[j];
    }
  }
  
  return List::create(Named("first_order_tilde") = first_order_tilde,
                      Named("second_order_tilde") = second_order_tilde);
}
