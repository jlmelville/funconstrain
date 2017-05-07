#' Meyer Function
#'
#' Test function 10 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 3}, number of summand
#'   functions \code{m = 16}.
#'   \item Minima: \code{f = 87.9458...}.
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
#' Meyer, R. R. (1970).
#' Theoretical and computational aspects of nonlinear regression.
#' In J. B. Rosen, O. L. Mangasarian, and K. Ritter (Eds.)
#' \emph{Nonlinear programming} (pp465-496).
#' New York: Academic Press.
#'
#' @examples
#' fun <- meyer()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
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

      ti <- 45 + 5 * (1:m)
      fi <- x1 * exp(x2 / (ti + x3)) - y

      sum(fi * fi)
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      ti <- 45 + 5 * (1:m)
      tix3 <- ti + x3
      tix3s <- tix3 * tix3
      e <- exp(x2 / tix3)
      fi <- x1 * e - y
      ef2 <- e * fi * 2

      dx <- sum(ef2)
      dy <- sum((x1 * ef2) / tix3)
      dz <- -sum((x1 * x2 * ef2) / tix3s)
      c(dx, dy, dz)
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      ti <- 45 + 5 * (1:m)
      tix3 <- ti + x3
      tix3s <- tix3 * tix3
      e <- exp(x2 / tix3)
      fi <- x1 * e - y
      fsum <- sum(fi * fi)

      ef2 <- e * fi * 2

      dx <- sum(ef2)
      dy <- sum((x1 * ef2) / tix3)
      dz <- -sum((x1 * x2 * ef2) / tix3s)
      grad <- c(dx, dy, dz)

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.02, 4000, 250)
  )
}
