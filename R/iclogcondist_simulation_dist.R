#' Simulate from a Truncated Weibull Distribution
#'
#' This function generates random samples from a truncated Weibull distribution
#' using inverse transform sampling. When \code{shape = 1}, it reduces to a truncated exponential distribution.
#' 
#'
#' @param n An integer specifying the number of random samples to generate.
#' @param shape A positive numeric value representing the shape parameter of the Weibull distribution. Default is \code{1}.
#' @param scale A positive numeric value representing the scale parameter of the Weibull distribution. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of \code{n} random samples from the truncated Weibull distribution.
#' @importFrom stats rweibull runif
#' @export
#' @examples
#' # Generate 10 random samples from a truncated Weibull distribution
#' rtweibull(10, shape = 2, scale = 1, upper_bound = 5)
#' 
rtweibull <- function(n, shape = 1, scale = 1, upper_bound = Inf) {
  if (any(c(n, shape, scale) <= 0)) stop("n, shape, and scale must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated Weibull distribution
    return(rweibull(n, shape, scale))
  } else {
    # Truncated Weibull distribution
    F_b <- 1 - exp(-(upper_bound / scale)^shape)  # Truncated CDF value at upper bound
    U <- runif(n)
    X <- scale * (-log(1 - U * F_b))^(1 / shape)
    return(X)
  }
}

#' Cumulative Distribution Function of a Truncated Weibull Distribution
#'
#' This function computes the cumulative distribution function (CDF) of a truncated Weibull distribution
#' at a given point \code{x}.
#'
#' @param x A numeric vector at which to evaluate the CDF.
#' @param shape A positive numeric value representing the shape parameter of the Weibull distribution. Default is \code{1}.
#' @param scale A positive numeric value representing the scale parameter of the Weibull distribution. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of the CDF values of the truncated Weibull distribution at \code{x}.
#' @importFrom stats pweibull
#' @export
#' @examples
#' # Evaluate the CDF at x = 2 for a truncated Weibull distribution
#' ptweibull(2, shape = 2, scale = 1, upper_bound = 5)
#' 
ptweibull <- function(x, shape = 1, scale = 1, upper_bound = Inf) {
  if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated Weibull CDF
    return(pweibull(x, shape, scale))
  } else {
    # Truncated Weibull CDF
    trunc_factor <- pweibull(upper_bound, shape, scale)
    p_vals <- pweibull(x, shape, scale) / trunc_factor
    p_vals[x > upper_bound] <- 1  # CDF is 1 for x >= upper_bound
    return(pmin(p_vals, 1))       # Ensure values don't exceed 1
  }
}

#' Quantile Function of a Truncated Weibull Distribution
#'
#' This function computes the quantiles of a truncated Weibull distribution for a given probability vector.
#'
#' @param q A numeric vector of probabilities for which to calculate the quantiles.
#' @param shape A positive numeric value representing the shape parameter of the Weibull distribution. Default is \code{1}.
#' @param scale A positive numeric value representing the scale parameter of the Weibull distribution. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of quantiles corresponding to the given probabilities in \code{q}.
#' @importFrom stats qweibull pweibull uniroot
#' @export
#' @examples
#' # Calculate the 0.5 quantile of a truncated Weibull distribution
#' qtweibull(0.5, shape = 2, scale = 1, upper_bound = 5)
#' 
qtweibull <- function(q, shape = 1, scale = 1, upper_bound = Inf) {
  if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  if (any(q < 0 | q > 1)) stop("Probabilities must be between 0 and 1")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated Weibull quantile
    return(qweibull(q, shape, scale))
  } else {
    # Truncated Weibull quantile
    trunc_factor <- pweibull(upper_bound, shape, scale)
    output <- sapply(q, function(prob) {
      uniroot(function(x) pweibull(x, shape, scale) / trunc_factor - prob, 
              interval = c(0, upper_bound))$root
    })
    return(output)
  }
}




#' Simulate from a Truncated Log-Logistic Distribution
#'
#' This function generates random samples from a truncated log-logistic distribution
#' using an acceptance-rejection method.
#'
#' @param n An integer specifying the number of random samples to generate.
#' @param shape A positive numeric value representing the shape parameter of the log-logistic distribution. Default is \code{1}.
#' @param scale A positive numeric value representing the scale parameter of the log-logistic distribution. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of \code{n} random samples from the truncated log-logistic distribution.
#' @importFrom flexsurv rllogis
#' @export
#' @examples
#' # Generate 10 random samples from a truncated log-logistic distribution
#' rtllogis(10, shape = 2, scale = 1, upper_bound = 5)
#' 
rtllogis <- function(n, shape = 1, scale = 1, upper_bound = Inf) {
  if (any(c(n, shape, scale) <= 0)) stop("n, shape, and scale must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated log-logistic distribution
    return(flexsurv::rllogis(n, shape, scale))
  } else {
    # Truncated log-logistic distribution
    temp <- flexsurv::rllogis(10 * n, shape, scale)
    while (sum(temp < upper_bound) < n) {
      temp <- c(temp, flexsurv::rllogis(10 * n, shape, scale))
    }
    temp <- temp[temp < upper_bound]
    return(temp[1:n])
  }
}

#' Cumulative Distribution Function of a Truncated Log-Logistic Distribution
#'
#' This function computes the cumulative distribution function (CDF) of a truncated log-logistic distribution
#' at a given point \code{x}.
#'
#' @param x A numeric vector at which to evaluate the CDF.
#' @param shape A positive numeric value representing the shape parameter of the log-logistic distribution. Default is \code{1}.
#' @param scale A positive numeric value representing the scale parameter of the log-logistic distribution. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of the CDF values of the truncated log-logistic distribution at \code{x}.
#' @importFrom flexsurv pllogis
#' @export
#' @examples
#' # Evaluate the CDF at x = 2 for a truncated log-logistic distribution
#' ptllogis(2, shape = 2, scale = 1, upper_bound = 5)
#' 
ptllogis <- function(x, shape = 1, scale = 1, upper_bound = Inf) {
  if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated log-logistic CDF
    return(flexsurv::pllogis(x, shape, scale))
  } else {
    # Truncated log-logistic CDF
    trunc_factor <- flexsurv::pllogis(upper_bound, shape, scale)
    p_vals <- flexsurv::pllogis(x, shape, scale) / trunc_factor
    p_vals[x > upper_bound] <- 1  # CDF is 1 for x >= upper_bound
    return(pmin(p_vals, 1))       # Ensure values don't exceed 1
  }
}

#' Quantile Function of a Truncated Log-Logistic Distribution
#'
#' This function computes the quantiles of a truncated log-logistic distribution for a given probability vector.
#'
#' @param q A numeric vector of probabilities for which to calculate the quantiles.
#' @param shape A positive numeric value representing the shape parameter of the log-logistic distribution. Default is \code{1}.
#' @param scale A positive numeric value representing the scale parameter of the log-logistic distribution. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of quantiles corresponding to the given probabilities in \code{q}.
#' @importFrom flexsurv pllogis
#' @importFrom stats uniroot
#' @export
#' @examples
#' # Calculate the 0.5 quantile of a truncated log-logistic distribution
#' qtllogis(0.5, shape = 2, scale = 1, upper_bound = 5)
#' 
qtllogis <- function(q, shape = 1, scale = 1, upper_bound = Inf) {
  if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  if (any(q < 0 | q > 1)) stop("Probabilities must be between 0 and 1")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated log-logistic quantile
    return(flexsurv::qllogis(q, shape, scale))
  } else {
    # Truncated log-logistic quantile
    trunc_factor <- flexsurv::pllogis(upper_bound, shape, scale)
    output <- sapply(q, function(prob) {
      uniroot(function(x) flexsurv::pllogis(x, shape, scale) / trunc_factor - prob, 
              interval = c(0, upper_bound))$root
    })
    return(output)
  }
}





#' Simulate from a Truncated Log-Normal Distribution
#'
#' This function generates random samples from a truncated log-normal distribution
#' using an acceptance-rejection method.
#'
#' @param n An integer specifying the number of random samples to generate.
#' @param meanlog A numeric value representing the mean of the log-normal distribution on the log scale. Default is \code{0}.
#' @param sdlog A positive numeric value representing the standard deviation of the log-normal distribution on the log scale. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of \code{n} random samples from the truncated log-normal distribution.
#' @importFrom stats  rlnorm
#' @export
#' @examples
#' # Generate 10 random samples from a truncated log-normal distribution
#' rtlnorm(10, meanlog = 0, sdlog = 1, upper_bound = 5)
#' 
rtlnorm <- function(n, meanlog = 0, sdlog = 1, upper_bound = Inf) {
  if (any(c(n, sdlog) <= 0)) stop("n and sdlog must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated log-normal distribution
    return(rlnorm(n, meanlog, sdlog))
  } else {
    # Truncated log-normal distribution
    temp <- rlnorm(10 * n, meanlog, sdlog)
    while (sum(temp < upper_bound) < n) {
      temp <- c(temp, rlnorm(10 * n, meanlog, sdlog))
    }
    temp <- temp[temp < upper_bound]
    return(temp[1:n])
  }
}

#' Cumulative Distribution Function of a Truncated Log-Normal Distribution
#'
#' This function computes the cumulative distribution function (CDF) of a truncated log-normal distribution
#' at a given point \code{x}.
#'
#' @param x A numeric vector at which to evaluate the CDF.
#' @param meanlog A numeric value representing the mean of the log-normal distribution on the log scale. Default is \code{0}.
#' @param sdlog A positive numeric value representing the standard deviation of the log-normal distribution on the log scale. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of the CDF values of the truncated log-normal distribution at \code{x}.
#' @importFrom stats plnorm
#' @export
#' @examples
#' # Evaluate the CDF at x = 2 for a truncated log-normal distribution
#' ptlnorm(2, meanlog = 0, sdlog = 1, upper_bound = 5)
#' 
ptlnorm <- function(x, meanlog = 0, sdlog = 1, upper_bound = Inf) {
  if (sdlog <= 0) stop("sdlog must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated log-normal CDF
    return(plnorm(x, meanlog, sdlog))
  } else {
    # Truncated log-normal CDF
    trunc_factor <- plnorm(upper_bound, meanlog, sdlog)
    p_vals <- plnorm(x, meanlog, sdlog) / trunc_factor
    p_vals[x > upper_bound] <- 1  # CDF is 1 for x >= upper_bound
    return(pmin(p_vals, 1))       # Ensure values don't exceed 1
  }
}

#' Quantile Function of a Truncated Log-Normal Distribution
#'
#' This function computes the quantiles of a truncated log-normal distribution for a given probability vector.
#'
#' @param q A numeric vector of probabilities for which to calculate the quantiles.
#' @param meanlog A numeric value representing the mean of the log-normal distribution on the log scale. Default is \code{0}.
#' @param sdlog A positive numeric value representing the standard deviation of the log-normal distribution on the log scale. Default is \code{1}.
#' @param upper_bound A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).
#' @return A numeric vector of quantiles corresponding to the given probabilities in \code{q}.
#' @importFrom stats qlnorm plnorm uniroot
#' @export
#' @examples
#' # Calculate the 0.5 quantile of a truncated log-normal distribution
#' qtlnorm(0.5, meanlog = 0, sdlog = 1, upper_bound = 5)
#' 
qtlnorm <- function(q, meanlog = 0, sdlog = 1, upper_bound = Inf) {
  if (sdlog <= 0) stop("sdlog must be positive")
  if (upper_bound <= 0) stop("upper_bound must be positive")
  if (any(q < 0 | q > 1)) stop("Probabilities must be between 0 and 1")
  
  if (is.infinite(upper_bound)) {
    # Non-truncated log-normal quantile
    return(qlnorm(q, meanlog, sdlog))
  } else {
    # Truncated log-normal quantile
    trunc_factor <- plnorm(upper_bound, meanlog, sdlog)
    output <- sapply(q, function(prob) {
      uniroot(function(x) plnorm(x, meanlog, sdlog) / trunc_factor - prob, 
              interval = c(0, upper_bound))$root
    })
    return(output)
  }
}


