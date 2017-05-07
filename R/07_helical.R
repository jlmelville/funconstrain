#' Helical Valley Function
#'
#' Test function 7 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 3}, number of summand
#'   functions \code{m = 3}.
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
#'
#' @examples
#' fun <- helical()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
helical <- function() {
  one_div_2pi <- 0.5 / pi

  theta <- function(x1, x2) {
    res <- one_div_2pi * atan(x2 / x1)
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
      pyyxx <- pi * (yy / xx  + 1)

      fx <- 10 * (z - 10 * theta(x, y))
      fy <- 10 * (sxy - 1)
      fz <- z

      dx <- (100 * y * fx) / (pyyxx * xx) + (20 * x * fy) / sxy
      dy <- (-100 * fx) / (pyyxx * x) + (20 * y * fy) / sxy
      dz <- 20 * fx + 2 * z

      c(dx, dy, dz)
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]
      z <- par[3]

      xx <- x * x
      yy <- y * y
      sxy <- sqrt(xx + yy)
      pyyxx <- pi * (yy / xx  + 1)

      fx <- 10 * (z - 10 * theta(x, y))
      fy <- 10 * (sxy - 1)
      fz <- z

      dx <- (100 * y * fx) / (pyyxx * xx) + (20 * x * fy) / sxy
      dy <- (-100 * fx) / (pyyxx * x) + (20 * y * fy) / sxy
      dz <- 20 * fx + 2 * z

      fsum <- fx * fx + fy * fy + fz * fz
      grad <- c(dx, dy, dz)

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(-1, 0, 0)
  )
}
