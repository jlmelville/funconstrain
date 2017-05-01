#' Beale Function
#'
#' Test function 5 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 2}, \code{m = 3}.
#'   \item Minima: \code{f = 0} at \code{(3, 0.5) },
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} The function.
#'   \item \code{gr} The gradient.
#'   \item \code{x0} Standard starting point.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Beale, E. M. L. (1958).
#' \emph{On an iterative method for finding a local minimum of a function of more
#' than one variable} (No. 25).
#' Statistical Techniques Research Group, Section of Mathematical Statistics,
#' Department of Mathematics, Princeton University.
#' @export
beale <- function() {
  list(
    fn = function(par) {
      x <- par[1]
      y <- par[2]
      yy <- y * y
      yyy <- yy * y

      f1 <- 1.5 - x * (1 - y)
      f2 <- 2.25 - x * (1 - yy)
      f3 <- 2.625 - x * (1 - yyy)

      f1 * f1 + f2 * f2 + f3 * f3
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]

      yy <- y * y
      yyy <- yy * y

      f1 <- 1.5 - x * (1 - y)
      f2 <- 2.25 - x * (1 - yy)
      f3 <- 2.625 - x * (1 - yyy)

      c(
        2 * (y - 1) * f1 + 2 * (yy - 1) * f2 + 2 * (yyy - 1) * f3,
        2 * x * f1 + 4 * x * y * f2 + 6 * x * yy * f3
      )
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]
      yy <- y * y
      yyy <- yy * y

      f1 <- 1.5 - x * (1 - y)
      f2 <- 2.25 - x * (1 - yy)
      f3 <- 2.625 - x * (1 - yyy)

      list(
        fn = f1 * f1 + f2 * f2 + f3 * f3,
        gr = c(
          2 * (y - 1) * f1 + 2 * (yy - 1) * f2 + 2 * (yyy - 1) * f3,
          2 * x * f1 + 4 * x * y * f2 + 6 * x * yy * f3
        )
      )
    },
    x0 = c(1, 1)
  )
}
