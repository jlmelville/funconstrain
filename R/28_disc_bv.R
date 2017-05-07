#' Discrete Boundary Value Function
#'
#' Test function 28 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n}.
#'   \item Minima: \code{f = 0}.
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
#' More', J. J., & Cosnard, M. Y. (1979).
#' Numerical solution of nonlinear equations.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{5}(1), 64-85.
#' \url{https://doi.org/10.1145/355815.355820}
#'
#' @examples
#' dbv <- disc_bv()
#' # 6 variable problem using the standard starting point
#' x0_6 <- dbv$x0(6)
#' res_6 <- stats::optim(x0_6, dbv$fn, dbv$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(dbv$x0(8), dbv$fn, dbv$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), dbv$fn, dbv$gr,
#'                       method = "L-BFGS-B")
#' @export
disc_bv <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Discrete Boundary Value: n must be positive")
      }
      h <- 1 / (n + 1)
      hsq <- h * h
      ti <- 1:n * h
      pt1 <- par + ti + 1

      fi <- 2 * par + hsq * 0.5 * pt1 * pt1 * pt1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      fsum <- sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Discrete Boundary Value: n must be positive")
      }
      grad <- rep(0, n)

      h <- 1 / (n + 1)
      hsq <- h * h
      ti <- 1:n * h
      pt1 <- par + ti + 1

      fi <- 2 * par + hsq * 0.5 * pt1 * pt1 * pt1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]

      grad <- 2 * fi * (1.5 * hsq * pt1 * pt1 + 2)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 2 * fi[1:(n - 1)]
      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Discrete Boundary Value: n must be positive")
      }

      h <- 1 / (n + 1)
      hsq <- h * h
      ti <- 1:n * h
      pt1 <- par + ti + 1

      fi <- 2 * par + hsq * 0.5 * pt1 * pt1 * pt1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      fsum <- sum(fi * fi)

      grad <- 2 * fi * (1.5 * hsq * pt1 * pt1 + 2)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 2 * fi[1:(n - 1)]

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 35) {
      if (n < 1) {
        stop("Discrete Boundary Value: n must be positive")
      }
      (1:n / n + 1) * ((1:n / n + 1) - 1)
    }
  )
}
