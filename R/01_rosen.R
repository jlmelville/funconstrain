#' Rosenbrock Function
#'
#' Test function 1 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 2}, number of summand
#'   functions \code{m = 2}.
#'   \item Minima: \code{f = 0} at \code{(1, 1)}
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
#' Rosenbrock, H. (1960).
#' An automatic method for finding the greatest or least value of a function.
#' \emph{The Computer Journal}, \emph{3}(3), 175-184.
#' \url{https://doi.org/10.1093/comjnl/3.3.175}
#'
#' @examples
#' fun <- rosen()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
rosen <- function() {
  list(
    fn = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
    },
    gr = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      c(
        -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
        200 *      (x2 - x1 * x1))
    },
    fg = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x2x11 <- x2 - x1 * x1
      list(
        fn = 100 * x2x11 ^ 2 + (1 - x1) ^ 2,
        gr = c(
          -400 * x1 * x2x11 - 2 * (1 - x1),
           200 * x2x11
        )
      )
    },
    x0 = c(-1.2, 1)
  )
}
