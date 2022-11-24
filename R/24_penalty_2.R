#' Penalty Function II
#'
#' Test function 24 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n + 1}.
#'   \item Minima: \code{f = 9.37629...e-6} if \code{n = 4};
#'   \code{f = 2.93660...e-4} if \code{n = 10}.
#' }
#'
#' The number of parameters, \code{n}, in the objective function is not
#' specified when invoking this function. It is implicitly set by the length of
#' the parameter vector passed to the objective and gradient functions that this
#' function creates. See the 'Examples' section..
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
#' Gill, P. E., Murray, W., & Pitfield, R. A. (1972).
#' \emph{The implementation of two revised quasi-Newton algorithms for
#' unconstrained optimization} (Report NAC-11).
#' Teddington, London: National Physical Laboratory, Division of Numerical
#' Analysis and Computing.
#' @examples
#' pen2 <- penalty_2()
#' # 6 variable problem using the standard starting point
#' x0_6 <- pen2$x0(6)
#' res_6 <- stats::optim(x0_6, pen2$fn, pen2$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(pen2$x0(8), pen2$fn, pen2$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), pen2$fn, pen2$gr,
#'                       method = "L-BFGS-B")
#' @export
penalty_2 <- function() {
  a <- 1e-5
  e <- exp(-0.1)

  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function II: n must be positive")
      }
      ex01 <- exp(par * 0.1)
      ei <- exp(1:n * 0.1)

      fsum <- 0

      f1 <- par[1] - 0.2
      fsum <- fsum + f1 * f1

      f2n <- n * par[1] * par[1]
      for (i in 2:n) {
        yi <- ei[i] + ei[i - 1]
        fi <- ex01[i] + ex01[i - 1] - yi

        fsum <- fsum + a * fi * fi
        f2n <- f2n + (n - i + 1) * par[i] * par[i]
      }
      for (i in (n + 1):(2 * n - 1)) {
        fi <- ex01[i - n + 1] - e
        fsum <- fsum + a * fi * fi
      }

      f2n <- f2n - 1
      fsum <- fsum + f2n * f2n
      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function II: n must be positive")
      }
      grad <- rep(0, n)
      grad[1] <- grad[1] + 2 * (par[1] - 0.2)
      ex01 <- exp(par * 0.1)
      ei <- exp(1:n * 0.1)

      f2n <- n * par[1] * par[1]
      for (i in 2:n) {
        yi <- ei[i] + ei[i - 1]
        fi <- (ex01[i] + ex01[i - 1] - yi)

        grad[i] <- grad[i] + 0.2 * fi * a * ex01[i]
        grad[i - 1] <- grad[i - 1] + 0.2 * fi * a * ex01[i - 1]

        f2n <- f2n + (n - i + 1) * par[i] * par[i]
      }

      for (i in (n + 1):(2 * n - 1)) {
        fi <- ex01[i - n + 1] - e
        grad[i - n + 1] <- grad[i - n + 1] + 0.2 * fi * a * ex01[i - n + 1]
      }

      f2n <- f2n - 1

      for (i in 1:n) {
        grad[i] <- grad[i] + 4 * (n - i + 1) * par[i] * f2n
      }

      grad
    },
    he = function(x) { # ?? failing? Why?
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)
       t1 <- -1.0
       for (j in 1:n){
          t1 <- t1 + ( n-j+1 ) * x[j]^2
       }
       d1 <- exp( 0.1 )
       d2 <- 1.0
       th <- 4.0*t1
       s2 <- 0.0
       for (j in 1:n){
          h[j,j] <- 8.0*( ( n-j+1 )*x[j] )^2 + ( n-j+1 )*th
          s1 <- exp( x[j]/10.0 )
          if ( j > 1 ) {
             s3 <- s1 + s2 - d2*( d1 + 1.0 )
             h[j  ,j  ] <- h[j  ,j  ] + 1.0e-5*s1*( s3 + s1 - 1.0/d1 + 2.0*s1 ) / 50.0
             h[j-1,j-1] <- h[j-1,j-1] + 1.0e-5*s2*( s2 + s3 ) / 50.0
             for (k in 1:(j-1)){
                t1 <- exp( k/10.0 )
                h[k,j] <- 8.0*( n-j+1 )*( n-k+1 )*x[j]*x[k]
             }
             h[j-1,j] <- h[j-1,j] + 1.0e-5*s1*s2/50.0
          }
          s2 <- s1
          d2 <- d1*d2
       }
       h[1,1] <- h[1,1] + 2.0
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
        stop("Penalty Function II: n must be positive")
      }
      ex01 <- exp(par * 0.1)
      ei <- exp(1:n * 0.1)

      fsum <- 0

      grad <- rep(0, n)
      grad[1] <- grad[1] + 2 * (par[1] - 0.2)

      f1 <- par[1] - 0.2
      fsum <- fsum + f1 * f1
      f2n <- n * par[1] * par[1]

      for (i in 2:n) {
        yi <- ei[i] + ei[i - 1]
        fi <- (ex01[i] + ex01[i - 1] - yi)
        fsum <- fsum + a * fi * fi

        grad[i] <- grad[i] + 0.2 * fi * a * ex01[i]
        grad[i - 1] <- grad[i - 1] + 0.2 * fi * a * ex01[i - 1]
        f2n <- f2n + (n - i + 1) * par[i] * par[i]
      }

      for (i in (n + 1):(2 * n - 1)) {
        fi <- ex01[i - n + 1] - e
        fsum <- fsum + a * fi * fi
        grad[i - n + 1] <- grad[i - n + 1] + 0.2 * fi * a * ex01[i - n + 1]
      }

      f2n <- f2n - 1
      fsum <- fsum + f2n * f2n

      for (i in 1:n) {
        grad[i] <- grad[i] + 4 * (n - i + 1) * par[i] * f2n
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 25) {
      if (n < 1) {
        stop("Penalty Function II: n must be positive")
      }
      rep(0.5, n)
    },
    fmin = 9.376293e-6,
    xmin = c(0.1999993, 0.19131669, 0.48010149, 0.5188454) # n=4 case
  )
}
