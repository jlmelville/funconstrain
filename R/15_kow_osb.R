#' Kowalik and Osborne Function
#'
#' Test function 15 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 4}, number of summand
#'   functions \code{m = 11}.
#'   \item Minima: \code{f = 3.07505...e-4};
#'   and \code{f = 1.02734...e-3} at \code{(Inf, -14.07..., -Inf, -Inf)}.
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
#' Kowalik, J. S., & Osborne, M. R. (1968).
#' \emph{Methods for unconstrained optimization problems.}
#' New York, NY: Elsevier North-Holland.
#'
#' @examples
#' fun <- kow_osb()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
#' "L-BFGS-B")
#' @export
kow_osb <- function() {
  m <- 11
  y <- c(0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342,
         0.0323, 0.0235, 0.0246)
  u <- c(4, 2, 1, 0.5, 0.25, 0.167, 0.125, 0.1, 0.0833, 0.0714, 0.0625)
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      fsum <- 0
      for (i in 1:m) {
        ui2 <- u[i] * u[i]
        num <- (ui2 + u[i] * x2)
        den <- (ui2 + u[i] * x3 + x4)
        fi <- y[i] - (x1 * num) / den
        fsum <- fsum + fi * fi
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      grad <- c(0, 0, 0, 0)
      for (i in 1:m) {
        ui2 <- u[i] * u[i]
        num <- ui2 + u[i] * x2
        den <- (ui2 + u[i] * x3 + x4)
        den2 <- den * den
        fi <- y[i] - (x1 * num) / den

        f2 <- 2 * fi
        nf2 <- f2 * num
        nfx2 <- nf2 * x1

        grad[1] <- grad[1] - nf2 / den
        grad[2] <- grad[2] - (f2 * u[i] * x1) / den
        grad[3] <- grad[3] + (nfx2 * u[i]) / den2
        grad[4] <- grad[4] + nfx2 / den2
      }
      grad
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      h <- matrix(0.0, ncol=4, nrow=4)
      for (i in 1:m) {
         s1 <- u[i] + x2
          t1 <- u[i] ^ 2 + u[i]*x3 + x4

          if ( t1 != 0.0 ) {
             d1 <- y[i] - x1*u[i]*s1/t1
             h[1,1] <- h[1,1] + 2.0*( u[i]*s1/t1 ) ^ 2
             h[1,2] <- h[1,2] + 2.0*u[i]/t1*( x1*u[i]*s1/t1 - d1 )
             h[2,2] <- h[2,2] + 2.0*( x1*u[i]/t1 ) ^ 2
             h[1,3] <- h[1,3] + 2.0*s1*u[i] ^ 2/t1 ^ 2*( d1 - u[i]*s1*x1/t1 )
             h[2,3] <- h[2,3] + 2.0*x1*u[i] ^ 2/t1 ^ 2*( d1 - x1*u[i]*s1/t1 )
             h[3,3] <- h[3,3] + 2.0*u[i] ^ 3*x1*s1/t1 ^ 3*( x1*u[i]*s1/t1 - 2.0*d1 )
             h[1,4] <- h[1,4] + 2.0*u[i]*s1/t1 ^ 2*( d1 - u[i]*s1*x1/t1 )
             h[2,4] <- h[2,4] + 2.0*x1*u[i]/t1 ^ 2*( d1 - x1*u[i]*s1/t1 )
             h[3,4] <- h[3,4] + 2.0*x1*u[i] ^ 2*s1/t1 ^ 3*( x1*u[i]*s1/t1 - 2.0*d1 )
             h[4,4] <- h[4,4] + 2.0*x1*u[i]*s1/t1 ^ 3*( x1*u[i]*s1/t1 - 2.0*d1 )
          } else {
             h<-matrix(.Machine$double.eps, nrow=4, ncol=4)
             # flag <- - 3
             return()
          } 
      }
      h[2,1] <- h[1,2]
      h[3,1] <- h[1,3]
      h[3,2] <- h[2,3]
      h[4,1] <- h[1,4]
      h[4,2] <- h[2,4]
      h[4,3] <- h[3,4]
      h
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      fsum <- 0
      grad <- c(0, 0, 0, 0)
      for (i in 1:m) {
        ui2 <- u[i] * u[i]
        num <- ui2 + u[i] * x2
        den <- (ui2 + u[i] * x3 + x4)
        den2 <- den * den
        fi <- y[i] - (x1 * num) / den
        f2 <- 2 * fi
        nf2 <- f2 * num
        nfx2 <- nf2 * x1

        fsum <- fsum + fi * fi

        grad[1] <- grad[1] - nf2 / den
        grad[2] <- grad[2] - (f2 * u[i] * x1) / den
        grad[3] <- grad[3] + (nfx2 * u[i]) / den2
        grad[4] <- grad[4] + nfx2 / den2
      }
      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.25, 0.39, 0.415, 0.39)
  )
}
