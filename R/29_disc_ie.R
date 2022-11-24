#' Discrete Integral Equation Function
#'
#' Test function 29 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n}.
#'   \item Minima: \code{f = 0} (at the same location as \code{\link{disc_bv}}).
#' }
#'
#' The number of parameters, \code{n}, in the objective function is not
#' specified when invoking this function. It is implicitly set by the length of
#' the parameter vector passed to the objective and gradient functions that this
#' function creates. See the 'Examples' section.
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
#'   \item \code{x0} Function returning the standard starting point, given
#'   \code{n}, the number of variables desired.
#'   \item \code{fmin} reported minimum
#'   \item \code{xmin} parameters at reported minimum
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \doi{doi.org/10.1145/355934.355936}
#'
#' More', J. J., & Cosnard, M. Y. (1979).
#' Numerical solution of nonlinear equations.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{5}(1), 64-85.
#' \doi{doi.org/10.1145/355815.355820}
#'
#' @examples
#' d_ie <- disc_ie()
#' # 6 variable problem using the standard starting point
#' x0_6 <- d_ie$x0(6)
#' res_6 <- stats::optim(x0_6, d_ie$fn, d_ie$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(d_ie$x0(8), d_ie$fn, d_ie$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), d_ie$fn, d_ie$gr,
#'                       method = "L-BFGS-B")
#' @export
disc_ie <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Discrete Integral Equation: n must be positive")
      }
      h <- 1 / (n + 1)

      pt1 <- par + (1:n * h) + 1
      pt13 <- pt1 * pt1 * pt1
      tii <- 1:n * h

      # for each i, sum1 loops from 1:i, so we start at 0 and accumulate by
      # adding the new i value at each iteration
      sum1 <- 0
      # for each i, sum2 loops from i+1:n, so here we start with the sum of 1:n
      # and subtract the ith value at each iteration
      sum2 <- sum((1 - tii) * pt13)

      fsum <- 0
      for (i in 1:n) {
        ti <- tii[i]
        ti1 <- 1 - ti

        sum1 <- sum1 + ti * pt13[i]
        sum2 <- sum2 - ti1 * pt13[i]

        fi <- par[i] + 0.5 * h * (ti1 * sum1 + ti * sum2)
        fsum <- fsum + fi * fi
      }
      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Discrete Integral Equation: n must be positive")
      }
      grad <- rep(0, n)

      h <- 1 / (n + 1)

      pt1 <- par + (1:n * h) + 1
      pt12 <- pt1 * pt1
      pt13 <- pt12 * pt1
      hp3 <- 3 * h * pt12
      tii <- 1:n * h

      sum1 <- 0
      sum2 <- sum((1 - tii) * pt13)

      for (i in 1:n) {
        ti <- tii[i]
        ti1 <- 1 - ti

        sum1 <- sum1 + ti * pt13[i]
        sum2 <- sum2 - ti1 * pt13[i]

        fi <- par[i] + 0.5 * h * (ti1 * sum1 + ti * sum2)

        fhp3 <- fi * hp3
        grad[1:i] <- grad[1:i] + fhp3[1:i] * (1 - tii[i]) * tii[1:i]
        grad[i] <- grad[i] + 2 * fi
        if (i < n) {
          r2 <- (i + 1):n
          grad[r2] <- grad[r2] + fhp3[r2] * (1 - tii[r2]) * tii[i]
        }
      }
      grad
    },
    he = function(x) { 
       n <- length(x)
       w1 <- rep(0, n)
       w2 <- rep(0, n+1)
       gvec <- rep(0,n)
       h <- matrix(0.0, nrow=n, ncol=n)
#       for (j in 1:n){
#          for (i in 1:j) {
#             h[i,j] <- 0.0
#          }
#       }
       
       d1 <- 1.0/(n + 1.0)
       w1[1]          <- d1*( x[1] + d1 + 1.0 )^3
       w2[n]   <- ( 1.0 - n*d1 )*( x[n] + n*d1 + 1.0 )^3
       w2[n+1] <- 0.0
       for (i in 2:n) {
          t1 <- i*d1
          t2 <- (n-i+1)*d1
          w1[i]            <- w1[i-1] + t1*( x[i] + t1 + 1.0)^3
          w2[n-i+1] <- w2[n-i+2] + ( 1.0 - t2 )*( x[n-i+1] + t2 + 1.0 )^3
       }

       for (i in 1:n){
#          cat("i:",i)
          t1 <- i*d1
          t  <- x[i] + 0.5*d1*( ( 1.0 - t1 )*w1[i] + t1*w2[i+1] )

          for (j in 1:i) {
#             cat(" j:",j)
             t2 <- j*d1
             gvec[j] <- 1.5*d1*( 1.0 - t1 )*t2*( x[j] + t2 + 1.0 )^2

             h[j,j] <- h[j,j] + 2.0*t*( 3.0*d1*( 1.0 - t1 )*t2*( x[j] + t2 + 1.0 ) )
          }

          gvec[i] <- gvec[i] + 1.0

#          cat("new j loop\n")
          if (i < 8) {
            for (j in (i+1):n) {
#               cat(" j:",j)
               t2 <- j*d1
               gvec[j] <- 1.5*d1*t1*( 1.0 - t2 )*( x[j] + t2 + 1.0 )^2
               h[j,j] <- h[j,j] + 2.0*t*( 3.0*d1*t1*(1.0 - t2)*(x[j] + t2 + 1.0))
            }
          }
#          cat("new k loop\n")
          for (k in 1:n) {
#             cat(" k:",j)
             for (j in 1:k) {
#             cat(" j:",j)
                h[j,k] <- h[j,k] + 2.0*gvec[j]*gvec[k]
             }
#             cat("\n")
          }
       }

       for (j in 1:(n-1)) { # symmetrize
         for (k in (j+1):n) {
           h[k,j] <- h[j,k]        
         }
       }
       h
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Discrete Integral Equation: n must be positive")
      }

      h <- 1 / (n + 1)

      pt1 <- par + (1:n * h) + 1
      pt12 <- pt1 * pt1
      pt13 <- pt12 * pt1
      hp3 <- 3 * h * pt12
      tii <- 1:n * h

      sum1 <- 0
      sum2 <- sum((1 - tii) * pt13)

      fsum <- 0
      grad <- rep(0, n)
      for (i in 1:n) {
        ti <- tii[i]
        ti1 <- 1 - ti

        sum1 <- sum1 + ti * pt13[i]
        sum2 <- sum2 - ti1 * pt13[i]

        fi <- par[i] + 0.5 * h * (ti1 * sum1 + ti * sum2)
        fsum <- fsum + fi * fi

        fhp3 <- fi * hp3
        grad[1:i] <- grad[1:i] + fhp3[1:i] * (1 - tii[i]) * tii[1:i]
        grad[i] <- grad[i] + 2 * fi
        if (i < n) {
          r2 <- (i + 1):n
          grad[r2] <- grad[r2] + fhp3[r2] * (1 - tii[r2]) * tii[i]
        }
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 35) {
      if (n < 1) {
        stop("Discrete Integral Equation: n must be positive")
      }
      (1:n / n + 1) * ((1:n / n + 1) - 1)
    },
    fmin = 0,
    xmin =     xmin = c(-0.07502213, -0.1319762, -0.1648488, -0.1646647, -0.1174177) # n=5 case
  )
}
