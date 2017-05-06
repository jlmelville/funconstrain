#' Brown Badly Scaled Function
#'
#' Test function 4 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 2}, number of summand
#'   functions \code{m = 3}.
#'   \item Minima: \code{f = 0} at \code{(1e6, 2e-6) },
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
#' @examples
#' fun <- brown_bs()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
brown_bs <- function() {
  list(
    fn = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- x - 1e6
      f2 <- y - 2e-6
      f3 <- x * y - 2

      f1 * f1 + f2 * f2 + f3 * f3
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]

      f3 <- x * y - 2

      c(
        2 * y * f3 + 2 * (x - 1e6),
        2 * x * f3 + 2 * (y - 2e-6)
      )
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- x - 1e6
      f2 <- y - 2e-6
      f3 <- x * y - 2

      list(
        fn = f1 * f1 + f2 * f2 + f3 * f3,
        gr = c(
          2 * y * f3 + 2 * f1,
          2 * x * f3 + 2 * f2
        )
      )
    },
    x0 = c(1, 1)
  )
}
