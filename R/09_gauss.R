#' Gaussian Function
#'
#' Test function 9 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 3}, number of summand
#'   functions \code{m = 15}.
#'   \item Minima: \code{f = 1.12793...e-8}.
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
#' @examples
#' fun <- gauss()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
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
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

       y9 <- c(0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521,
         0.3989, 0.3521, 0.2420, 0.1295, 0.0540, 0.0175, 0.0044, 0.0009)      
       h <- matrix(0.0, ncol=3, nrow=3)
       for (i in 1:m) {
          d1 <- 0.5*(i-1)
          d2 <- 3.5 - d1 - x3
          arg <- 0.5*x2*d2 ^ 2
          r <- exp( - arg )
          t <- x1*r - y9[i]
          t1 <- 2.0*x1*r - y9[i]
          h[1,1] <- h[1,1] + r ^ 2
          h[1,2] <- h[1,2] - r*t1*d2 ^ 2
          h[2,2] <- h[2,2] + r*t1*d2 ^ 4
          h[1,3] <- h[1,3] + d2*r*t1
          h[2,3] <- h[2,3] + d2*r*( t - arg*t1 )
          h[3,3] <- h[3,3] + r*( x2*t1*d2 ^ 2 - t )
       }
       h[1,1] <- 2.0*h[1,1]
       h[2,2] <- 0.5*x1*h[2,2]
       h[1,3] <- 2.0*x2*h[1,3]
       h[2,3] <- 2.0*x1*h[2,3]
       h[3,3] <- 2.0*x1*x2*h[3,3]
       h[2,1] <- h[1,2]
       h[3,1] <- h[1,3]
       h[3,2] <- h[2,3]
       h
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
    x0 = c(0.4, 1, 0),
    fmin = 1.12793e-8,
    xmin = c( 0.3989561, 1.0000191, 2.787451e-20) # APPROX!
  )
}
