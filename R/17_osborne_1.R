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
#'   \item \code{he} If available, the hessian matrix (second derivatives)
#'   of the function w.r.t. the parameters at the given values.
#'   \item \code{fg} A function which, given the parameter vector, calculates
#'   both the objective value and gradient, returning a list with members
#'   \code{fn} and \code{gr}, respectively.
#'   \item \code{x0} Standard starting point.
#'   \item \code{fmin} reported minimum
#'   \item \code{xmin} parameters at reported minimum
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \doi{doi.org/10.1145/355934.355936}
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

      ti <- 10 * 0:(m - 1)
      fi <- y - (x1 + x2 * exp(-ti * x4) + x3 * exp(-ti * x5))
      sum(fi * fi)

    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]

      ti <- 10 * 0:(m - 1)
      a <- exp(-ti * x4)
      b <- exp(-ti * x5)
      f <- y - (x1 + x2 * a + x3 * b)
      f2 <- 2 * f
      fa2 <- f2 * a
      fb2 <- f2 * b

      c(
        -sum(f2),
        -sum(fa2),
        -sum(fb2),
        sum(fa2 * ti * x2),
        sum(fb2 * ti * x3)
      )

    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      h <- matrix(0.0, ncol=5, nrow=5)
      for (i in 1:m) {
         t1 <- 10.0*(i-1)
         d1 <- exp( - t1*x4 )
         d2 <- exp( - t1*x5 )
         s1 <- y[i] - ( x1 + x2*d1 + x3*d2 )
          h[1,1] <- h[1,1] + 2.0
          h[1,2] <- h[1,2] + 2.0*d1
          h[2,2] <- h[2,2] + 2.0*d1 ^ 2
          h[1,3] <- h[1,3] + 2.0*d2
          h[2,3] <- h[2,3] + 2.0*d1*d2
          h[3,3] <- h[3,3] + 2.0*d2 ^ 2
          h[1,4] <- h[1,4] - 2.0*t1*x2*d1
          h[2,4] <- h[2,4] + 2.0*t1*d1*( s1 - x2*d1 )
          h[3,4] <- h[3,4] - 2.0*t1*x2*d1*d2
          h[4,4] <- h[4,4] + 2.0*t1 ^ 2*x2*d1*( x2*d1 - s1 )
          h[1,5] <- h[1,5] - 2.0*t1*x3*d2
          h[2,5] <- h[2,5] - 2.0*t1*x3*d1*d2
          h[3,5] <- h[3,5] + 2.0*t1*d2*( s1 - x3*d2 )
          h[4,5] <- h[4,5] + 2.0*t1 ^ 2*d1*d2*x2*x3
          h[5,5] <- h[5,5] + 2.0*t1 ^ 2*x3*d2*( x3*d2 - s1 )
      }
      h[2,1] <- h[1,2]
      h[3,1] <- h[1,3]
      h[3,2] <- h[2,3]
      h[4,1] <- h[1,4]
      h[4,2] <- h[2,4]
      h[4,3] <- h[3,4]
      h[5,1] <- h[1,5]
      h[5,2] <- h[2,5]
      h[5,3] <- h[3,5]
      h[5,4] <- h[4,5]
      h
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]

      ti <- 10 * 0:(m - 1)
      a <- exp(-ti * x4)
      b <- exp(-ti * x5)
      f <- y - (x1 + x2 * a + x3 * b)
      f2 <- 2 * f
      fa2 <- f2 * a
      fb2 <- f2 * b

      fsum <- sum(f * f)
      grad <- c(
        -sum(f2),
        -sum(fa2),
        -sum(fb2),
        sum(fa2 * ti * x2),
        sum(fb2 * ti * x3)
      )

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.5, 1.5, -1, 0.01, 0.02),
    fmin = 5.464895e-5,
    xmin = c(0.3754101, 1.935847, -1.4646871, 0.01286753, 0.02212270)
  )
}
