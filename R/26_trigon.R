#' Trigonometric Function
#'
#' Test function 26 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n}.
#'   \item Minima: \code{f = 0} at \code{rep(0, n)}. But there are many others.
#' }
#'
#' The number of parameters, \code{n}, in the objective function is not
#' specified when invoking this function. It is implicitly set by the length of
#' the parameter vector passed to the objective and gradient functions that this
#' function creates. See the 'Examples' section.
#'
#' @note It seems to be extremely difficult to reach the global minimum for
#' any \code{n} from the starting location using the typical gradient-based
#' methods (e.g. conjugate gradient, BFGS, L-BFGS).
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
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Spedicato, E. (1975).
#' \emph{Computational experience with quasi-Newton algorithms for minimization
#' problems of moderately large size} (Report CISE-N-175).
#' Segrate, Milano: Computing Center, CISE.
#'
#' @examples
#' trig <- trigon()
#' # 6 variable problem using the standard starting point
#' x0_6 <- trig$x0(6)
#' res_6 <- stats::optim(x0_6, trig$fn, trig$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(trig$x0(8), trig$fn, trig$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), trig$fn, trig$gr,
#'                       method = "L-BFGS-B")
#' @export
trigon <- function() {
  list(
    m = length(par),
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }

      cos_sum <- sum(cos(par))
      fi <- n - cos_sum + 1:n * (1 - cos(par)) - sin(par)
      sum(fi * fi)

    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }
      cosx <- cos(par)
      sinx <- sin(par)
      cos_sum <- sum(cosx)
      fi <- n - cos_sum + 1:n * (1 - cosx) - sinx

      2 * (fi * (1:n * sinx - cosx) + sinx * sum(fi))
    },
    he = function(x) { 
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)

       s1 <- 0.0
       for (j in 1:n) {
          h[j,j] <- sin( x[j] )
          s1 <- s1 + cos( x[j] )
       }
       s2 <- 0.0
       for (j in 1:n) {
          th <- cos( x[j] )
          t <- ( n+j ) - h[j,j] - s1 - j*th
          s2 <- s2 + t
          if (j > 1) {
            for (k in (1:(j-1))){
              h[k,j] <- 2.0*(sin(x[k])*(( n+j+k )*h[j,j]-th) - h[j,j]*cos(x[k]) )
            }
          }
          h[j,j] <- (j*(j+2)+n)*h[j,j]^2 + 
               th*(th-(2*j+2)*h[j,j]) + t*(j*th + h[j,j] )
       }

       for (j in 1:n) {
          h[j,j] <- 2.0*( h[j,j] + cos(x[j])*s2 )
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
        stop("Trigonometric: n must be positive")
      }
      cosx <- cos(par)
      sinx <- sin(par)
      cos_sum <- sum(cosx)
      fi <- n - cos_sum + 1:n * (1 - cosx) - sinx

      fsum <- sum(fi * fi)
      grad <- 2 * (fi * (1:n * sinx - cosx) + sinx * sum(fi))

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 30) {
      if (n < 1) {
        stop("Trigonometric: n must be positive")
      }
      rep(1 / n, n)
    }
  )
}

