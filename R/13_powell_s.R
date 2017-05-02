#' Powell Singular Function
#'
#' Test function 13 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 4}, \code{m = 4}.
#'   \item Minima: \code{f = 0} at \code{rep(0, 4)}.
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
#' Powell, M. J. D. (1962).
#' An iterative method for finding stationary values of a function of several
#' variables.
#' \emph{The Computer Journal}, \emph{5}(2), 147-151.
#' \url{https://doi.org/10.1093/comjnl/5.2.147}
#' @export
powell_s <- function() {
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      f1 <- x1 + 10 * x2
      f2 <- sqrt(5) * (x3 - x4)
      f3 <- (x2 - 2 * x3) ^ 2
      f4 <- sqrt(10) * (x1 - x4) ^ 2

      f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      x12 <- x1 + 10 * x2
      x14 <- x1 - x4
      x14_3 <- x14 * x14 * x14
      x34 <- x3 - x4
      x23 <- x2 - 2 * x3
      x23_3 <- x23 * x23 * x23

      c(
        40 * x14_3 + 2 * x12,
        4 * x23_3 + 20 * x12,
        10 * x34 - 8 * x23_3,
        -10 * x34 - 40 * x14_3
      )
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      x14 <- x1 - x4
      x14_2 <- x14 * x14
      x14_3 <- x14_2 * x14
      x34 <- x3 - x4
      x23 <- x2 - 2 * x3
      x23_2 <- x23 * x23
      x23_3 <- x23_2 * x23
      x12 <- x1 + 10 * x2

      f1 <- x12
      f2 <- sqrt(5) * x34
      f3 <- x23_2
      f4 <- sqrt(10) * x14_2

      list(
        fn = f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4,
        gr = c(
          40 * x14_3 + 2 * x12,
          4 * x23_3 + 20 * x12,
          10 * x34 - 8 * x23_3,
          -10 * x34 - 40 * x14_3
        )
      )
    },
    x0 = c(3, -1, 0, 1)
  )
}