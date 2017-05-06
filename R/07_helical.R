#' Helical Valley Function
#'
#' Test function 7 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 3}, \code{m = 3}.
#'   \item Minima: \code{f = 0} at \code{(1, 0, 0)}.
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
#' Fletcher, R., & Powell, M. J. (1963).
#' A rapidly convergent descent method for minimization.
#' \emph{The Computer Journal}, \emph{6}(2), 163-168.
#' @export
helical <- function() {
  theta <- function(x1, x2) {
    res <- 0.5 * atan(x2 / x1) / pi
    if (x1 < 0) {
      res <- res + 0.5
    }
    res
  }

  list(
    fn = function(par) {
      x <- par[1]
      y <- par[2]
      z <- par[3]

      f1 <- 10 * (z - 10 * theta(x, y))
      f2 <- 10 * (sqrt(x * x + y * y) - 1)
      f3 <- z

      f1 * f1 + f2 * f2 + f3 * f3
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]
      z <- par[3]

      xx <- x * x
      yy <- y * y
      sxy <- sqrt(xx + yy)
      sxy1 <- sxy - 1
      txy <- theta(x, y)
      ztxy <- z - 10 * txy
      pyyxx <- pi * (yy / xx  + 1)

      c(
        1000 * y * 10 * ztxy / (pyyxx * xx) + 200 * x * sxy1 / sxy,
        200 * y * sxy1 / sxy - 1000 * ztxy / (pyyxx * x),
        200 * ztxy + 2 * z
      )
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]
      z <- par[3]

      xx <- x * x
      yy <- y * y
      sxy <- sqrt(xx + yy)
      sxy1 <- sxy - 1
      txy <- theta(x, y)
      ztxy <- z - 10 * txy
      pyyxx <- pi * (yy / xx  + 1)

      f1 <- 10 * ztxy
      f2 <- 10 * sxy1
      f3 <- z

      fsum <- f1 * f1 + f2 * f2 + f3 * f3
      grad <- c(
        1000 * y * 10 * ztxy / (pyyxx * xx) + 200 * x * sxy1 / sxy,
        200 * y * sxy1 / sxy - 1000 * ztxy / (pyyxx * x),
        200 * ztxy + 2 * z
      )

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(-1, 0, 0)
  )
}
