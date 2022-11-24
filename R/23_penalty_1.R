#' Penalty Function I
#'
#' Test function 23 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n + 1}.
#'   \item Minima: \code{f = 2.24997...e-5} if \code{n = 4};
#'   \code{f = 7.08765...e-5} if \code{n = 10}.
#' }
#'
#' The number of parameters, \code{n}, in the objective function is not
#' specified when invoking this function. It is implicitly set by the length of
#' the parameter vector passed to the objective and gradient functions that this
#' function creates. See the 'Examples' section.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Objective function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Gradient function which calculates the gradient vector
#'   given input parameter vector.
#'   \item \code{he} If available, the hessian matrix (second derivatives)
#'   of the function w.r.t. the parameters at the given values.
#'   \item \code{fg} A function which, given the parameter vector, calculates
#'   both the objective value and gradient, returning a list with members
#'   \code{fn} and \code{gr}, respectively.
#'   \item \code{x0} Function returning the standard starting point, given
#'   \code{n}, the number of variables desired.
#'   \item \code{fmin} reported minimum
#'   \item \code{xmin} parameters at reported minimum
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \doi{doi.org/10.1145/355934.355936}
#'
#' Gill, P. E., Murray, W., & Pitfield, R. A. (1972).
#' \emph{The implementation of two revised quasi-Newton algorithms for
#' unconstrained optimization} (Report NAC-11).
#' Teddington, London: National Physical Laboratory, Division of Numerical
#' Analysis and Computing.
#' @examples
#' pen1 <- penalty_1()
#' # 6 variable problem using the standard starting point
#' x0_6 <- pen1$x0(6)
#' res_6 <- stats::optim(x0_6, pen1$fn, pen1$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(pen1$x0(8), pen1$fn, pen1$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), pen1$fn, pen1$gr,
#'                       method = "L-BFGS-B")
#' @export
penalty_1 <- function() {
  a <- 1e-5
  sqrta <- sqrt(a)
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }
      fsum <- 0
      fn1 <- 0
      for (i in 1:n) {
        fi <- sqrta * (par[i] - 1)
        fsum <- fsum + fi * fi
        fn1 <- fn1 + par[i] * par[i]
      }

      fn1 <- fn1 - 0.25
      fsum <- fsum + fn1 * fn1
      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }
      grad <- rep(0, n)
      fn1 <- 0

      for (i in 1:n) {
        fi <- sqrta * (par[i] - 1)
        grad[i] <- grad[i] + 2 * sqrta * fi
        fn1 <- fn1 + par[i] * par[i]
      }
      fn1 <- fn1 - 0.25
      grad <- grad + 4 * par * fn1
      grad
    },
    he = function(x) {
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)
       t1 <- -0.25
       for (j in 1:n) {
          t1 <- t1 + x[j] ^ 2
       }
       d1 <- 2.0e-5
       th <- 4.0*t1
       for (j in 1:n) {
         for (k in 1:(j-1)) {
             h[k,j] <- 8.0*x[j]*x[k]
         }
          h[j,j] <- d1 + th + 8.0*x[j] ^ 2
          ## ! h[j,j) <- th + 8.0*x(j) ^ 2 - 1.0
       }
       for (j in 1:(n-1)) { # symmetrize
         for (k in (j+1):n) {
           h[k,j] <- h[j,k]        
         }
       }
       h
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }

      fn1 <- 0
      grad <- rep(0, n)
      fsum <- 0
      for (i in 1:n) {
        fi <- sqrta * (par[i] - 1)
        fsum <- fsum + fi * fi
        grad[i] <- grad[i] + 2 * sqrta * fi
        fn1 <- fn1 + par[i] * par[i]
      }
      fn1 <- fn1 - 0.25
      fsum <- fsum + fn1 * fn1
      grad <- grad + 4 * par * fn1

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 25) {
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }
      1:n
    },
    fmin = 2.24997e-5, 
    xmin = c(0.2500075, 0.2500075, 0.2500075, 0.2500075) # n=4 case
  )
}
