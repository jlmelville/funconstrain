#' Brown Almost-Linear Function
#'
#' Test function 27 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable, \code{m = n}.
#'   \item Minima: \code{f = 0} at \code{(a, a, a, ..., a ^ (1 - n))}, where
#'   \code{a} satisfies \code{n * a ^ n - (n + 1) * a ^ (n - 1) + 1 = 0};
#'   \code{f = 1} at \code{c(0, 0, ..., n + 1)}.
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
#' Brown, K. M. (1969).
#' A quadratically convergent Newton-like method based upon Gaussian
#' elimination.
#' \emph{SIAM Journal on Numerical Analysis}, \emph{6}(4), 560-569.
#' \url{http://dx.doi.org/10.1137/0706051}
#'
#' @examples
#' bal <- brown_al()
#' # 6 variable problem using the standard starting point
#' x0_6 <- bal$x0(6)
#' res_6 <- stats::optim(x0_6, bal$fn, bal$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(bal$x0(8), bal$fn, bal$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), bal$fn, bal$gr,
#'                       method = "L-BFGS-B")
#' @export
brown_al <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      fi <- par + sum(par) - (n + 1)
      fi[n] <- prod(par) - 1
      sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      fi <- par + sum(par) - (n + 1)
      grad <- rep(2 * sum(fi[1:(n - 1)]), n)
      grad[1:(n - 1)] <- grad[1:(n - 1)] + 2 * fi[1:(n - 1)]
      prod_x <- prod(par)
      if (prod_x > 0) {
        fn <- prod_x - 1
        grad <- grad + 2 * fn * (prod_x / par)
      }
      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      fi <- par + sum(par) - (n + 1)
      grad <- rep(2 * sum(fi[1:(n - 1)]), n)
      grad[1:(n - 1)] <- grad[1:(n - 1)] + 2 * fi[1:(n - 1)]
      prod_x <- prod(par)
      fn <- prod_x - 1

      fi[n] <- prod_x - 1
      fsum <- sum(fi * fi)
      if (prod_x > 0) {
        grad <- grad + 2 * fn * (prod_x / par)
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n) {
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      rep(0.5, n)
    }
  )
}
