#' Gaussian Function
#'
#' Test function 9 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 3}, \code{m = 15}.
#'   \item Minima: \code{f = 1.12793...e-8}.
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
gauss <- function() {
  y <- c(0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521, 0.3989,
         0.3521, 0.2420, 0.1295, 0.0540, 0.0175, 0.0044, 0.0009)
  m <- 15
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      for (i in 1:m) {
        ti <- (8 - i) * 0.5
        f <- x1 * exp(-0.5 * x2 * (ti - x3) ^ 2) - y[i]
        fsum <- fsum + f * f
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      grad <- c(0, 0, 0)
      for (i in 1:m) {
        ti <- (8 - i) * 0.5
        tx3 <- ti - x3
        tx3s <- tx3 * tx3
        g <- exp(-0.5 * x2 * tx3s)
        x1g <- x1 * g
        f <- x1g - y[i]

        grad[1] <- grad[1] + 2 * g * f
        grad[2] <- grad[2] - x1g * tx3s * f
        grad[3] <- grad[3] + 2 * x1g * x2 * tx3 * f
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      grad <- c(0, 0, 0)
      for (i in 1:m) {
        ti <- (8 - i) * 0.5
        tx3 <- ti - x3
        tx3s <- tx3 * tx3
        g <- exp(-0.5 * x2 * tx3s)
        x1g <- x1 * g
        f <- x1g - y[i]

        fsum <- fsum + f * f
        grad[1] <- grad[1] + 2 * g * f
        grad[2] <- grad[2] - x1g * tx3s * f
        grad[3] <- grad[3] + 2 * x1g * x2 * tx3 * f
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.4, 1, 0)
  )
}
