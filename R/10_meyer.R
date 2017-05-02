#' Meyer Function
#'
#' Test function 10 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 3}, \code{m = 16}.
#'   \item Minima: \code{f = 87.9458...}.
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
#' Meyer, R. R. (1970).
#' Theoretical and computational aspects of nonlinear regression.
#' In J. B. Rosen, O. L. Mangasarian, and K. Ritter (Eds.)
#' \emph{Nonlinear programming} (pp465-496).
#' New York: Academic Press.
#' @export
meyer <- function() {
  y <- c(34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030,
         6005, 5147, 4427, 3820, 3307, 2872)
  m <- 16
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      for (i in 1:m) {
        ti <- 45 + 5 * i
        f <- x1 * exp(x2 / (ti + x3)) - y[i]
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
        ti <- 45 + 5 * i
        tix3 <- ti + x3
        tix3s <- tix3 * tix3
        g <- exp(x2 / tix3)
        f <- x1 * g - y[i]
        gf2 <- g * f * 2
        grad[1] <- grad[1] + gf2
        grad[2] <- grad[2] + (x1 * gf2) / tix3
        grad[3] <- grad[3] - (x1 * x2 * gf2) / tix3s
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
        ti <- 45 + 5 * i
        tix3 <- ti + x3
        g <- exp(x2 / tix3)
        f <- x1 * g - y[i]
        fsum <- fsum + f * f

        tix3s <- tix3 * tix3
        gf2 <- g * f * 2
        grad[1] <- grad[1] + gf2
        grad[2] <- grad[2] + (x1 * gf2) / tix3
        grad[3] <- grad[3] - (x1 * x2 * gf2) / tix3s
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.02, 4000, 250)
  )
}