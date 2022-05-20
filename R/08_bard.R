#' Bard Function
#'
#' Test function 8 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 3}, number of summand
#'   functions \code{m = 15}.
#'   \item Minima: \code{f = 8.21487...e-3} and
#'   \code{f = 17.4286} at \code{(0.8406, -Inf, -Inf)}.
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
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \doi{doi.org/10.1145/355934.355936}
#'
#' Bard, Y. (1970).
#' Comparison of gradient methods for the solution of nonlinear parameter
#' estimation problems.
#' \emph{SIAM Journal on Numerical Analysis}, \emph{7}(1), 157-186.
#' \doi{dx.doi.org/10.1137/0707011}
#'
#' @examples
#' fun <- bard()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
bard <- function() {
  y <- c(0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58,
         0.73, 0.96, 1.34, 2.10, 4.39)
  m <- 15
  n <- 3
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
    he = function(par) {
       y8 <- c(0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 
          0.58, 0.73, 0.96, 1.34, 2.10, 4.39 )

      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      h <- matrix(0.0, nrow=3, ncol=3)
      for (i in 1:m) {
         d1 <- i
         d2 <- 16.0 - i
         d3 <- min( d1, d2 )

         s1 <- d2*x2 + d3*x3
         if ( s1 != 0.0 ) {
             t1 <- y8[i] - ( x1 + ( d1 / s1 ) )
             t2 <- d1 / s1 ^ 2
             t3 <- - 2.0*d1 / s1 ^ 3
             s2 <- t1*t3 + t2 ^ 2
             h[1,1] <- h[1,1] + 2.0
             h[1,2] <- h[1,2] - 2.0*d2*t2
             h[2,2] <- h[2,2] + 2.0*s2*d2 ^ 2
             h[1,3] <- h[1,3] - 2.0*t2*d3
             h[2,3] <- h[2,3] + 2.0*s2*d2*d3
             h[3,3] <- h[3,3] + 2.0*s2*d3 ^ 2
             
         } else {
             h[1:n, 1,] <- .Machine$double.xmax
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
