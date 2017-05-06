#' Osborne 1 Function
#'
#' Test function 17 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 5}, number of summand
#'   functions \code{m = 33}.
#'   \item Minima: \code{f = 5.46489...e-5}.
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
#' Osborne, M. R. (1972).
#' Some aspects of nonlinear least squares calculations.
#' In F. A. Lootsma (Ed.),
#' \emph{Numerical methods for nonlinear optimization} (pp. 171-189).
#' New York, NY: Academic Press.
#'
#' @examples
#' fun <- osborne_1()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(1:5, fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
osborne_1 <- function() {
  m <- 33
  y <- c(0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850, 0.818, 0.784,
         0.751, 0.718, 0.685, 0.658, 0.628, 0.603, 0.580, 0.558, 0.538, 0.522,
         0.506, 0.490, 0.478, 0.467, 0.457, 0.448, 0.438, 0.431, 0.424, 0.420,
         0.414, 0.411, 0.406)
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]

      fsum <- 0
      for (i in 1:m) {
        ti <- 10 * (i - 1)
        f <- y[i] - (x1 + x2 * exp(-ti * x4) + x3 * exp(-ti * x5))
        fsum <- fsum + f * f
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]

      grad <- c(0, 0, 0, 0, 0)
      for (i in 1:m) {
        ti <- 10 * (i - 1)
        a <- exp(-ti * x4)
        b <- exp(-ti * x5)
        f <- y[i] - (x1 + x2 * a + x3 * b)
        f2 <- 2 * f
        fa2 <- f2 * a
        fb2 <- f2 * b

        grad[1] <- grad[1] - f2
        grad[2] <- grad[2] - fa2
        grad[3] <- grad[3] - fb2
        grad[4] <- grad[4] + fa2 * ti * x2
        grad[5] <- grad[5] + fb2 * ti * x3
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]

      fsum <- 0
      grad <- c(0, 0, 0, 0, 0)
      for (i in 1:m) {
        ti <- 10 * (i - 1)
        a <- exp(-ti * x4)
        b <- exp(-ti * x5)
        f <- y[i] - (x1 + x2 * a + x3 * b)

        fsum <- fsum + f * f

        f2 <- 2 * f
        fa2 <- f2 * a
        fb2 <- f2 * b

        grad[1] <- grad[1] - f2
        grad[2] <- grad[2] - fa2
        grad[3] <- grad[3] - fb2
        grad[4] <- grad[4] + fa2 * ti * x2
        grad[5] <- grad[5] + fb2 * ti * x3
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.5, 1.5, -1, 0.01, 0.02)
  )
}
