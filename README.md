# iclogcondist: Nonparametric Estimation for a Log-concave Distribution Function with Interval-censored Data in R

## Overview

The iclogcondist package provides an efficient algorithm to compute the nonparametric maximum likelihood estimator (NPMLE) of a log-concave distribution function for the underlying event time in mixed-case interval-censored data. The computational algorithm combines an active set method with an iterative convex minorant algorithm.

## Features

- **Distribution Function Estimation**: The package provides a function `ic_LCMLE` to estimate the distribution function assuming its log-concavity, as well as functions to compute the unconstrained MLE and the distribution function corresponding to the least concave majorant (LCM) of the logarithm of the unconstrained MLE.
  
- **Simulation Capabilities**: Users can simulate interval-censored datasets where the underlying event times follows some common parametric families with a log-concave distribution function.

- **Visualization**: The package includes plotting functions to compare the estimated distribution functions obtained from the log-concave MLE, the unconstrained MLE, and the distribution function corresponding to the LCM of the logarithm of the unconstrained MLE.


## Installation

You can clone the repository [iclogcondist](https://github.com/ChaoyuYuan/iclogcondist) locally, by running this command in your terminal

```sh
git clone https://github.com/ChaoyuYuan/iclogcondist.git
```

### Option 1: Build and Install Locally

If you have cloned the repository locally, follow these steps to build and install the package manually:

1. **Build the package:** 

  In R, navigate to the package directory and build the package:

```r
setwd("path/to/iclogcondist")  # Change to your package directory
devtools::document()     # Generate documentation
devtools::build()        # Build the package
```
2. **Install the package:**

  After building, you will get the package object `iclogcondist_1.0.0.tar.gz`. You can install the package locally:

```r
install.packages("path/to/iclogcondist_1.0.0.tar.gz", repos = NULL, type = "source")
```
or

```r
devtools::install_local("path/to/iclogcondist_1.0.0.tar.gz")
```
The second command will install the prerequisite packages, such as Rcpp and ggplots, automatically during the process, while the first command will not.


### Option 2: Install Directly from Local Files
You can also install the package directly from the local directory without building manually:

```r
devtools::install("path/to/iclogcondist")
```
Similarly, this command will install the prerequisite packages automatically.

## Loading the Package
Once the package is installed, load it with:

```r
library(iclogcondist)
```


## Usage Example
Here is a simple example of how to use the `iclogcondist` package:

```r
data(lgnm)
X <- lgnm

# Run the LC MLE algorithm
result <- ic_LCMLE(X)

# Print the estimated distribution function
print(result$est$F_hat)
```
For more details in usage examples, please refer to the [iclogcondist_Example.pdf](https://github.com/ChaoyuYuan/iclogcondist/blob/main/iclogcondist_Example.pdf) file.


## License
This package is licensed under the GPL-3 License.

## Contributing
Contributions are welcome! If you would like to report issues or contribute to the development, please feel free to open an issue or submit a pull request.