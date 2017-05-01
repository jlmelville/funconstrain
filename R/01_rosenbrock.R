#' Rosenbrock Function
#'
#' Test function 1 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 2}, \code{m = 2}.
#'   \item Minima: \code{f = 0} at \code{(1, 1)}
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Function which calculates the gradient vector given input
#'   parameter vector.
#'   \item \code{fg} A function which calculates both the value and gradient,
#'   given the parameter vector, returning a list.
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
#' @export
rosenbrock <- function() {
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
