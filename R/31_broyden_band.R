#' Broyden Banded Function
#'
#' Test function 31 from the More', Garbow and Hillstrom paper.
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
#'   \item \code{fn} Function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Function which calculates the gradient vector given input
#'   parameter vector.
#'   \item \code{fg} Function which calculates both the value and gradient,
#'   given the parameter vector, returning a list.
#'   \item \code{x0} Function returning the standard starting point, given
#'   \code{n}, the number of variables desired.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Broyden, C. G. (1971).
#' The convergence of an algorithm for solving sparse nonlinear systems.
#' \emph{Mathematics of Computation}, \emph{25}(114), 285-294.
#' \url{https://doi.org/10.1090/S0025-5718-1971-0297122-5}
#'
#' @examples
#' btri <- broyden_band()
#' # 6 variable problem using the standard starting point
#' x0_6 <- btri$x0(6)
#' res_6 <- stats::optim(x0_6, btri$fn, btri$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(btri$x0(8), btri$fn, btri$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), btri$fn, btri$gr,
#'                       method = "L-BFGS-B")
#' @export
broyden_band <- function() {

  ml <- 5
  mu <- 1
  # return the valid indexes in a band around i
  #  the range i-mlow:i+mupp inclusive, excluding i, and j < 1 or j > n
  band_j <- function(i, n, mlow = ml, mupp = mu) {
    j <- max(1, i - mlow):min(n, i + mupp)
    j[j != i]
  }

  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }

      xx <- par * par
      fi <- (2 + 5 * xx) * par + 1
      fj <- par + xx
      for (i in 1:n) {
        fi[i] <- fi[i] - sum(fj[band_j(i, n)])
      }
      sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }

      xx <- par * par
      fi <- (2 + 5 * xx) * par + 1
      fj <- par + xx
      for (i in 1:n) {
        fi[i] <- fi[i] - sum(fj[band_j(i, n)])
      }

      grad <- 2 * fi * (15 * xx + 2)
      gj <- 2 * par + 1
      for (i in 1:n) {
        # reversing the limits returns the bands that i is part of
        grad[i] <- grad[i] - 2 * gj[i] * sum(fi[band_j(i, n, mu, ml)])
      }

      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }

      xx <- par * par
      fi <- (2 + 5 * xx) * par + 1
      fj <- par + xx
      for (i in 1:n) {
        fi[i] <- fi[i] - sum(fj[band_j(i, n)])
      }
      fsum <- sum(fi * fi)

      grad <- 2 * fi * (15 * xx + 2)
      gj <- 2 * par + 1
      for (i in 1:n) {
        grad[i] <- grad[i] - 2 * gj[i] * sum(fi[band_j(i, n, mu, ml)])
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n) {
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }
      rep(-1, n)
    }
  )
}
