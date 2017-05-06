#' Extended Rosenbrock Function
#'
#' Test function 21 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable but even, \code{m = n}.
#'   \item Minima: \code{f = 0} at \code{rep(1, n)}
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
#' Spedicato, E. (1975).
#' \emph{Computational experience with quasi-Newton algorithms for minimization
#' problems of moderately large size} (Report CISE-N-175).
#' Segrate, Milano: Computing Center, CISE.
#' @examples
#' exros <- ex_rosen()
#' # 6 variable problem using the standard starting point
#' x0_6 <- exros$x0(6)
#' res_6 <- stats::optim(x0_6, exros$fn, exros$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(exros$x0(8), exros$fn, exros$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), exros$fn, exros$gr,
#'                       method = "L-BFGS-B")
#' @export
ex_rosen <- function() {

  list(
    fn = function(par) {
      n <- length(par)
      if (n %% 2 != 0) {
        stop("Extended Rosenbrock: n must be even")
      }

      fsum <- 0
      for (i in 1:(n / 2)) {
        p2 <- 2 * i
        p1 <- p2 - 1

        f_p1 <- 10 * (par[p2] - par[p1] ^ 2)
        f_p2 <- 1 - par[p1]
        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
      }

      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n %% 2 != 0) {
        stop("Extended Rosenbrock: n must be even")
      }

      grad <- rep(0, n)
      for (i in 1:(n / 2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]

        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]

        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
      }

      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n %% 2 != 0) {
        stop("Extended Rosenbrock: n must be even")
      }

      fsum <- 0
      grad <- rep(0, n)
      for (i in 1:(n / 2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]

        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]

        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2

        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n) {
      if (n %% 2 != 0) {
        stop("Extended Rosenbrock: n must be even")
      }
      rep(c(-1.2, 1), n / 2)
    }
  )
}
