#' Trigonometric Function
#'
#' Test function 26 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable, \code{m = n}.
#'   \item Minima: \code{f = 0} at \code{rep(0, n)}. But there are many others.
#' }
#'
#' The number of variables in the function, \code{n}, is determined by the
#' length of the vector passed to the function and gradient routines. See
#' the 'Examples' section.
#'
#' @note It seems to be extremely difficult to reach the global minimum for
#' any \code{n} from the starting location using the typical gradient-based
#' methods (e.g. conjugate gradient, BFGS, L-BFGS).
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
#' Spedicato, E. (1975).
#' \emph{Computational experience with quasi-Newton algorithms for minimization
#' problems of moderately large size} (Report CISE-N-175).
#' Segrate, Milano: Computing Center, CISE.
#'
#' @examples
#' trig <- trigon()
#' # 6 variable problem using the standard starting point
#' x0_6 <- trig$x0(6)
#' res_6 <- stats::optim(x0_6, trig$fn, trig$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(trig$x0(8), trig$fn, trig$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), trig$fn, trig$gr,
#'                       method = "L-BFGS-B")
#' @export
trigon <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }

      cos_sum <- sum(cos(par))
      fi <- n - cos_sum + 1:n * (1 - cos(par)) - sin(par)
      sum(fi * fi)

    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }
      cosx <- cos(par)
      sinx <- sin(par)
      cos_sum <- sum(cosx)
      fi <- n - cos_sum + 1:n * (1 - cosx) - sinx

      2 * (fi * (1:n * sinx - cosx) + sinx * sum(fi))
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }
      cosx <- cos(par)
      sinx <- sin(par)
      cos_sum <- sum(cosx)
      fi <- n - cos_sum + 1:n * (1 - cosx) - sinx

      fsum <- sum(fi * fi)
      grad <- 2 * (fi * (1:n * sinx - cosx) + sinx * sum(fi))

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n) {
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }
      rep(1 / n, n)
    }
  )
}

