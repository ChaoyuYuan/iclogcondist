#' LGNM Data: Case II Interval Censoring Example
#'
#' This dataset, \code{lgnm}, provides an example of case II interval censoring data for illustrating the functions in this package. 
#' The event time is simulated from a log-normal distribution with parameters mean = 0 and standard deviation = 1.
#' The left censoring time is drawn from a uniform distribution between 0 and 2, and the right censoring time is drawn 
#' from a uniform distribution between the left censoring time and 20. Both the left and right censoring times are rounded to four decimal places.
#'
#' @name lgnm
#' @docType data
#' @format A data frame with 100 observations on the following 2 variables:
#' \describe{
#'     \item{left}{The left censoring time.}
#'     \item{right}{The right censoring time.}
#' }
#' @keywords datasets
#' @source Synthetic data generated for illustration purposes.
#' @examples
#' data(lgnm)
#' head(lgnm)
NULL
