#' Bard Function
#'
#' Test function 8 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 3}, \code{m = 15}.
#'   \item Minima: \code{f = 8.21487...e-3} and
#'   \code{f = 17.4286} at \code{(0.8406, -Inf, -Inf)}.
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
#' Bard, Y. (1970).
#' Comparison of gradient methods for the solution of nonlinear parameter
#' estimation problems.
#' \emph{SIAM Journal on Numerical Analysis}, \emph{7}(1), 157-186.
#' \url{http://dx.doi.org/10.1137/0707011}
#' @export
bard <- function() {
  y <- c(0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58,
         0.73, 0.96, 1.34, 2.10, 4.39)
  m <- 15

  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      for (u in 1:m) {
        v <- 16 - u
        w <- min(u, v)
        f <- y[u] - (x1 + u / (v * x2 + w * x3))
        fsum <- fsum + f * f
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      grad <- c(0, 0, 0)
      for (u in 1:m) {
        v <- 16 - u
        w <- min(u, v)
        a <- v * x2 + w * x3
        aa <- a * a

        f <- y[u] - (x1 + u / a)
        f2 <- 2 * f

        grad[1] <- grad[1] - f2
        grad[2] <- grad[2] + f2 * u * v / aa
        grad[3] <- grad[3] + f2 * u * w / aa
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      grad <- c(0, 0, 0)
      for (u in 1:m) {
        v <- 16 - u
        w <- min(u, v)
        a <- v * x2 + w * x3
        aa <- a * a

        f <- y[u] - (x1 + u / a)
        fsum <- fsum + f * f

        f2 <- 2 * f

        grad[1] <- grad[1] - f2
        grad[2] <- grad[2] + f2 * u * v / aa
        grad[3] <- grad[3] + f2 * u * w / aa
      }
      grad

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(1, 1, 1)
  )
}
