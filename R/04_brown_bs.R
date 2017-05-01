#' Brown Badly Scaled Function
#'
#' Test function 4 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 2}, \code{m = 3}.
#'   \item Minima: \code{f = 0} at \code{(1e6, 2e-6) },
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
