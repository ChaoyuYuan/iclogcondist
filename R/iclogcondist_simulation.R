#' Construct Case I Interval Censoring Data (Current Status Data)
#' 
#' This function constructs case I interval-censored data (current status data)
#' using the provided event times and censoring (survey) times.
#' Each individual's event time is either left-censored or right-censored at their survey time,
#' depending on whether the event has occurred by the survey time.
#'
#' @param event_times A numeric vector of event times for each individual.
#' @param survey_times A numeric vector of censoring (survey) times for each individual.
#' @return A matrix with two columns, where each row represents an individual's interval-censored data.
#'         The first column is the left endpoint, and the second column is the right endpoint.
#'         If the event time is before the survey time, the interval is \code{(0, survey_time]}.
#'         If the event time is after the survey time, the interval is \code{(survey_time, Inf)}.
current_status_X <- function(event_times, survey_times) {
  # Input validation
  if (length(event_times) != length(survey_times)) stop("event_times and survey_times must be of equal length")
  if (any(event_times < 0) || any(survey_times < 0)) stop("event_times and survey_times must be non-negative")
  
  n <- length(event_times)
  censor <- matrix(0, nrow = n, ncol = 2)
  
  # Assign intervals based on event times relative to survey times
  index <- which(event_times < survey_times)
  censor[index, 2] <- survey_times[index]  # Right-censored for events before survey time
  censor[-index, 1] <- survey_times[-index]  # Left-censored for events after survey time
  censor[-index, 2] <- Inf  # Right endpoint is Inf for right-censored intervals
  
  return(censor)
}




#' Construct Case II Interval Censoring Data
#'
#' This function constructs case II interval-censored data using the provided event times and censoring (survey) times.
#' Each individual's event time is either left-censored, right-censored, or interval-censored
#' based on two survey times: the left and right bounds of the interval.
#'
#' @param event_times A numeric vector of event times for each individual.
#' @param survey_times A numeric matrix with two columns, where each row contains the left and right
#'        censoring (survey) times for each individual.
#' @return A matrix with two columns, where each row represents an individual's interval-censored data.
#'         The first column is the left endpoint, and the second column is the right endpoint.
#'         If the event time is before the left survey time, the interval is \code{(0, left survey time]}.
#'         If the event time is after the right survey time, the interval is \code{(right survey time, Inf)}.
#'         If the event time falls between the left and right survey times, the interval is \code{(left survey time, right survey time]}.
case_II_X <- function(event_times, survey_times) {
  # Input validation
  if (length(event_times) != nrow(survey_times)) stop("event_times and survey_times must have the same number of rows")
  if (ncol(survey_times) != 2) stop("survey_times must have two columns (left and right censoring times)")
  if (any(event_times < 0) || any(survey_times < 0)) stop("event_times and survey_times must be non-negative")
  if (any(survey_times[, 1] > survey_times[, 2])) stop("Left survey times must be less than or equal to right survey times")
  
  n <- length(event_times)
  censor <- matrix(0, nrow = n, ncol = 2)
  
  # Case: Event before the left survey time (right-censored)
  index <- event_times <= survey_times[, 1]
  censor[index, 2] <- survey_times[index, 1]
  
  # Case: Event after the right survey time (left-censored)
  index2 <- event_times > survey_times[, 2]
  censor[index2, 1] <- survey_times[index2, 2]
  censor[index2, 2] <- Inf
  
  # Case: Event within the interval (interval-censored)
  index3 <- (event_times > survey_times[, 1]) & (event_times <= survey_times[, 2])
  censor[index3, ] <- survey_times[index3, ]
  
  return(censor)  
}




#' Simulate Interval-Censored Data
#'
#' This function generates interval-censored data, where the event times
#' are generated from one of the following distributions: Weibull, log-normal and log-logistic.
#' It supports both case 1 and case 2 interval censoring.
#'
#' @param n An integer specifying the number of observations to generate.
#' @param dist A character string indicating the distribution to use for event times. 
#'             Options are \code{"lognormal"}, \code{"weibull"}, or \code{"loglogistic"}.
#' @param para1 A numeric value representing the first parameter of the distribution:
#'              \itemize{
#'                \item \code{"lognormal"}: Mean of the log-normal distribution (meanlog).
#'                \item \code{"weibull"} and \code{"loglogistic"}: Shape parameter.
#'              }
#' @param para2 A numeric value representing the second parameter of the distribution:
#'              \itemize{
#'                \item \code{"lognormal"}: Standard deviation of the log-normal distribution (sdlog).
#'                \item \code{"weibull"} and \code{"loglogistic"}: Scale parameter.
#'              }
#' @param upper_bound A numeric value specifying the upper bound for event times, 
#' corresponding to a truncated distribution. Default is \code{Inf}.
#' @param C1_upper A numeric value specifying the upper limit for the first censoring time \code{C1}. 
#' Default is 1.
#' @param case An integer specifying the censoring case to simulate:
#'             \itemize{
#'               \item \code{1}: Current status (case 1 interval censoring)
#'               \item \code{2}: Case 2 Interval censoring
#'             }
#' @param rounding A logical value. If \code{TRUE}, generated times are rounded to a specified number of decimal places. Default is \code{FALSE}.
#' @param round_digit An integer specifying the number of digits for rounding when \code{rounding = TRUE}. Default is 4.
#' 
#' @return A matrix of interval-censored data where each row represents an interval (L, R] containing the unobserved event time.
#'
#' @details
#' \itemize{
#'   \item **Censoring Times**:
#'         \itemize{
#'           \item In \code{case = 1} (current status), one censoring time is generated, where it follows \code{U(0, C1_upper)}.
#'           \item In \code{case = 2} (case 2 interval censoring), two censoring times are generated:
#'                 \itemize{
#'                   \item \code{C1}: sampled from \code{U(0, C1_upper)}.
#'                   \item \code{C2}: sampled from \code{U(C1, min(upper_bound, 20))}.
#'                 }
#'         }
#'   \item **Distributions**:
#'         \itemize{
#'           \item **Weibull**: Parameterized by shape (\code{para1}) and scale (\code{para2}).
#'           \item **Log-logistic**: Parameterized by shape (\code{para1}) and scale (\code{para2}).
#'           \item **Log-normal**: Parameterized by mean (\code{para1}) and standard deviation (\code{para2}).
#'         }
#' }
#' 
#' @examples
#' # Simulate data with a truncated Weibull distribution and case II interval censoring
#' simulate_ic_data(n = 100, dist = "weibull", para1 = 2, para2 = 1, upper_bound = 5, case = 2)
#' @export
#' @importFrom stats runif
#' 
simulate_ic_data <- function(n, dist, para1, para2, upper_bound = Inf, C1_upper = 1, case = 2, 
                     rounding = FALSE, round_digit = 4) {
  
  # Input validation
  if (dist == "lognormal") {
    if (para2 <= 0) stop("For 'lognormal' distribution, 'para2' (sdlog) must be positive.")
  } else if (dist %in% c("weibull", "loglogistic")) {
    if (para1 <= 0 || para2 <= 0) stop("For 'weibull' and 'loglogistic' distributions, both 'para1' (shape) and 'para2' (scale) must be positive.")
  } else {
    stop("Invalid distribution. Choose 'lognormal', 'weibull', or 'loglogistic'.")
  }
  if (any(c(n, upper_bound) <= 0)) stop("n, and upper_bound must be positive")
  if (!case %in% c(1, 2)) stop("case must be 1 (current status) or 2 (interval censoring)")
  if (!is.logical(rounding)) stop("rounding must be a logical value (TRUE or FALSE)")
  
  # Generate event times based on specified distribution
  event_times <- switch(dist,
                        weibull = rtweibull(n, shape = para1, scale = para2, upper_bound = upper_bound),
                        loglogistic = rtllogis(n, shape = para1, scale = para2, upper_bound = upper_bound),
                        lognormal = rtlnorm(n, meanlog = para1, sdlog = para2, upper_bound = upper_bound))
  
  # Case I: Current status interval censoring
  if (case == 1) {
    survey_times <- if (rounding) round(runif(n, 0, C1_upper), round_digit) else runif(n, 0, C1_upper)
    X <- current_status_X(event_times, survey_times)
    return(X)
  }
  
  # Case II: Interval censoring
  if (case == 2) {
    C1 <- if (rounding) round(runif(n, 0, C1_upper), round_digit) else runif(n, 0, C1_upper)
    C2 <- if (rounding) round(runif(n, C1, min(upper_bound, 20)), round_digit) else runif(n, C1, min(upper_bound, 20))
    
    survey_times <- cbind(C1, C2)
    X <- case_II_X(event_times, survey_times)
    return(X)
  }
}
