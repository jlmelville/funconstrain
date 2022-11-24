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
#'   \item Minima: the MGH (1981) only provides the optimal value, with \code{f
#'   = 87.9458...}. Meyer and Roth (1972) give the optimal parameter values as
#'   \code{(0.0056, 6181.4, 345.2)}, with \code{f = 88}.
#' }
#'
#' @note
#' The gradient is large even at the optimal value, and enormous if using the
#' level of precision given by Meyer and Roth (smallest gradient component is at
#' least 1e4). It is not recommended to rely on the typical gradient norm
#' termination conditions if using this test function.
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
#' Meyer, R. R. (1970).
#' Theoretical and computational aspects of nonlinear regression.
#' In J. B. Rosen, O. L. Mangasarian, and K. Ritter (Eds.)
#' \emph{Nonlinear programming} (pp465-496).
#' New York: Academic Press.
#'
#' Meyer, R. R., & Roth, P. M. (1972).
#' Modified damped least squares: an algorithm for non-linear estimation.
#' \emph{IMA Journal of Applied Mathematics}, \emph{9}(2), 218-233.
#' \doi{doi.org/10.1093/imamat/9.2.218}
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
    he = function(par) {
#      y10 <- c(34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0,
#           8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0)
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      h <- matrix(0.0, ncol=3, nrow=3)
      for (i in 1:m) {
        d1 <- 4.5e+01 + 5.0*i
        d2 <- d1 + x3
        if ( d2 != 0.0 ) {
          t1 <- exp( x2 / d2 )
          t2 <- x1*t1 - y[i]
          s1 <- t2 + t1*x1
          h[1,1] <- h[1,1] + 2.0*t1 ^ 2
          h[1,2] <- h[1,2] + 2.0*s1*t1/d2
          h[2,2] <- h[2,2] + 2.0*t1*s1*x1/d2 ^ 2
          h[1,3] <- h[1,3] - 2.0*t1*s1*x2/d2 ^ 2
          h[2,3] <- h[2,3] - 2.0*t1*x1/d2 ^ 2*( t2 + s1*x2/d2 )
          h[3,3] <- h[3,3] + 2.0*t1*x1*x2/d2 ^ 3*( 2.0*t2 + s1*x2/d2 )
        } else {
          h <- matrix(0.0, ncol=3, nrow=3)
          # flag <- - 3
          return()
        } 
      } # end loop
      h[2,1] <- h[1,2]
      h[3,1] <- h[1,3]
      h[3,2] <- h[2,3]
      h
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
    x0 = c(0.02, 4000, 250),
    fmin = 87.9458,
    xmin = c(0.0056, 6181.4, 345.2) # Meyer and Roth (1972) APPROX!!!
  )
}
