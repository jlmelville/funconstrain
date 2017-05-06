#' Broyden Tridiagonal Function
#'
#' Test function 30 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable, \code{m = n}.
#'   \item Minima: \code{f = 0}.
#' }
#'
#' The number of variables in the function, \code{n}, is determined by the
#' length of the vector passed to the function and gradient routines. See
#' the 'Examples' section.
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
#' Broyden, C. G. (1965).
#' A class of methods for solving nonlinear simultaneous equations.
#' \emph{Mathematics of computation}, \emph{19}(92), 577-593.
#' \url{https://doi.org/10.2307/2003941}
#'
#' @examples
#' btri <- broyden_tri()
#' # 6 variable problem using the standard starting point
#' x0_6 <- btri$x0(6)
#' res_6 <- stats::optim(x0_6, btri$fn, btri$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(btri$x0(8), btri$fn, btri$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), btri$fn, btri$gr,
#'                       method = "L-BFGS-B")
#' @export
broyden_tri <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]

      grad <- 2 * fi * (3 - 4 * par)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 4 * fi[1:(n - 1)]

      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      fsum <- sum(fi * fi)

      grad <- 2 * fi * (3 - 4 * par)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 4 * fi[1:(n - 1)]

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n) {
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }
      rep(-1, n)
    }
  )
}
