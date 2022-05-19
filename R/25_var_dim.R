#' Variably Dimensioned Function
#'
#' Test function 25 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n + 2}.
#'   \item Minima: \code{f = 0} at \code{rep(1, n)}.
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
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' @examples
#' vdim <- var_dim()
#' # 6 variable problem using the standard starting point
#' x0_6 <- vdim$x0(6)
#' res_6 <- stats::optim(x0_6, vdim$fn, vdim$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(vdim$x0(8), vdim$fn, vdim$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), vdim$fn, vdim$gr,
#'                       method = "L-BFGS-B")
#' @export
var_dim <- function() {
  list(
    m = length(par) + 1,
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Variably Dimensioned: n must be positive")
      }
      fsum <- 0
      fn1 <- 0
      for (j in 1:n) {
        fj <- par[j] - 1
        fsum <- fsum + fj * fj
        fn1 <- fn1 + j * fj
      }
      fn1_fn1 <- fn1 * fn1
      # f_n+1 and f_n+2
      fsum <- fsum + fn1_fn1 + fn1_fn1 * fn1_fn1

      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Variably Dimensioned: n must be positive")
      }

      fsum <- 0
      grad <- rep(0, n)
      fn1 <- 0
      for (j in 1:n) {
        fj <- par[j] - 1
        fn1 <- fn1 + j * fj
        grad[j] <- grad[j] + 2 * fj
      }

      fn1_2 <- fn1 * 2
      fn13_4 <- fn1_2 * fn1_2 * fn1
      grad <- grad + 1:n * (fn1_2 + fn13_4)

      grad
    },
    he = function(x) { # rather large discrepancy from approximations in n, n
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)
       t1 <- 0.0
       for (j in 1:n) {
          t1 <- t1 + j*( x[j]-1.0 )
       }
       t <- 1.0 + 6.0*t1^2 # ?? why?
       for (j in 1:n) {
          h[j,j] <- 2.0 + 2.0 * t * j^2
          if (j > 1) {
            for (k in 1:(j-1)) {
               h[k,j] <- 2.0*t*j*k
            }
          }
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
        stop("Variably Dimensioned: n must be positive")
      }

      fsum <- 0
      grad <- rep(0, n)
      fn1 <- 0
      for (j in 1:n) {
        fj <- par[j] - 1
        fn1 <- fn1 + j * fj
        fsum <- fsum + fj * fj
        grad[j] <- grad[j] + 2 * fj
      }

      fn1_fn1 <- fn1 * fn1
      # f_n+1 and f_n+2
      fsum <- fsum + fn1_fn1 + fn1_fn1 * fn1_fn1

      fn1_2 <- fn1 * 2
      fn13_4 <- fn1_fn1 * fn1_2 * 2
      grad <- grad + 1:n * (fn1_2 + fn13_4)

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 30) {
      if (n < 1) {
        stop("Variably Dimensioned: n must be positive")
      }
      1 - (1:n) / n
    }
  )
}
