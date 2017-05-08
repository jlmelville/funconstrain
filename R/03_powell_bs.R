#' Powell Badly Scaled Function
#'
#' Test function 3 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 2}, number of summand
#'   functions \code{m = 2}.
#'   \item Minima: \code{f = 0} at \code{(1.098... e-5, 9.106...) },
#' }
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
#'   \item \code{x0} Standard starting point.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Powell, M. J. D. (1970).
#' A hybrid method for nonlinear equations.
#' In P. Rabinowitz (Ed.),
#' \emph{Numerical Methods for Nonlinear Algebraic Equations} (pp87-114).
#' New York: Gordon & Breach.
#'
#' @examples
#' fun <- powell_bs()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
powell_bs <- function() {
  list(
    fn = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- 1e4 * x * y - 1
      f2 <- exp(-x) + exp(-y) - 1.0001

      f1 * f1 + f2 * f2
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]
      a <- 2e4 * (1e4 * x * y - 1)
      b <- 2 * (exp(-x) + exp(-y) - 1.0001)

      c(
          y * a - exp(-x) * b,
          x * a - exp(-y) * b
      )
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- 1e4 * x * y - 1
      f2 <- exp(-x) + exp(-y) - 1.0001

      a <- 2e4 * f1
      b <- 2 * f2

      list(
        fn = f1 * f1 + f2 * f2,
        gr = c(
          y * a - exp(-x) * b,
          x * a - exp(-y) * b
        )
      )
    },
    x0 = c(0, 1)
  )
}
