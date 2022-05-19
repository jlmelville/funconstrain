#' Brown Almost-Linear Function
#'
#' Test function 27 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of summand
#'   functions \code{m = n}.
#'   \item Minima: \code{f = 0} at \code{(a, a, a, ..., a ^ (1 - n))}, where
#'   \code{a} satisfies \code{n * a ^ n - (n + 1) * a ^ (n - 1) + 1 = 0};
#'   \code{f = 1} at \code{c(0, 0, ..., n + 1)}.
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
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Brown, K. M. (1969).
#' A quadratically convergent Newton-like method based upon Gaussian
#' elimination.
#' \emph{SIAM Journal on Numerical Analysis}, \emph{6}(4), 560-569.
#' \url{http://dx.doi.org/10.1137/0706051}
#'
#' @examples
#' bal <- brown_al()
#' # 6 variable problem using the standard starting point
#' x0_6 <- bal$x0(6)
#' res_6 <- stats::optim(x0_6, bal$fn, bal$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(bal$x0(8), bal$fn, bal$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), bal$fn, bal$gr,
#'                       method = "L-BFGS-B")
#' @export
brown_al <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      fi <- par + sum(par) - (n + 1)
      fi[n] <- prod(par) - 1
      sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      fi <- par + sum(par) - (n + 1)
      grad <- rep(2 * sum(fi[1:(n - 1)]), n)
      grad[1:(n - 1)] <- grad[1:(n - 1)] + 2 * fi[1:(n - 1)]
      prod_x <- prod(par)
      if (prod_x > 0) {
        fn <- prod_x - 1
        grad <- grad + 2 * fn * (prod_x / par)
      }
      grad
    },
    he = function(x) { 
       n <- length(x)
       m <- n
       h <- matrix(0.0, nrow=n, ncol=n)
       pp <- matrix(0.0, nrow=(n+1), ncol=(n+1))
       # pp[i,j] = prod(i-1, j)  i,j in 1,(n+1)
       # prod(i,j) = pp[i+1,j] i in 0:n  j in 1:(n+1)
       for (j in (1:n)) { #        do j = 1, global_n
          pp[j,j] <- 1  #          prod(j-1,j) = 1.0_rk
          for (k in (j:n)){ #      do k = j, global_n
            pp[k+1,j] <- pp[k,j]*x[k] #   prod(k,j) = prod(k-1,j)*x(k)
          } #          end do
          pp[j+1,n+1] <- 1.0 #          prod(j,global_n+1) = 1.0_rk
       } #       end do

#       cat("pp\n")
#       print(pp)

       for (i in 1:m) {#       do i = 1, global_m
          if ( i == n ) {#          if ( i .eq. global_n ) then
             for (j in 1:n) {#             do j = 1, global_n
                t1 <- pp[n+1,1] #          t1 = prod(global_n,1)
                t <- t1 - 1.0 #                t = t1 - 1.0_rk
                h[j,j] <- h[j,j] + 2.0*( pp[j,1]*pp[n+1,j+1])^2
#                h(j,j) = h(j,j) + 2.0_rk*( prod(j-1,1)*prod(global_n,j+1) )**2
                for (l in 1:(j-1)) { #                do l = 1, j-1
                   t2 <- pp[l,1]*pp[j,l+1]*pp[n+1,j+1]
#                  t2 = prod(l-1,1)*prod(j-1,l+1)*prod(global_n,j+1)
                   h[l,j] <- h[l,j] + 2.0*t2*( 2.0*t1 - 1.0 )
#                  h(l,j) = h(l,j) + 2.0_rk*t2*( 2.0_rk*t1 - 1.0_rk )
                }  #                end do
             }
          } else {
             for (j in 1:n) {
                for (k in 1:j) {
                   if ( (j == i) && (k == i) ) {
                      h[k,j] <- h[k,j] + 8.0
                   }
                   else if ( (j == i) || (k == i) ){
                      h[k,j] <- h[k,j] + 4.0
                   } else {
                      h[k,j] <- h[k,j] + 2.0
                   }
                }
             }
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
        stop("Brown Almost-Linear: n must be positive")
      }
      fi <- par + sum(par) - (n + 1)
      grad <- rep(2 * sum(fi[1:(n - 1)]), n)
      grad[1:(n - 1)] <- grad[1:(n - 1)] + 2 * fi[1:(n - 1)]
      prod_x <- prod(par)
      fn <- prod_x - 1

      fi[n] <- prod_x - 1
      fsum <- sum(fi * fi)
      if (prod_x > 0) {
        grad <- grad + 2 * fn * (prod_x / par)
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 30) {
      if (n < 1) {
        stop("Brown Almost-Linear: n must be positive")
      }
      rep(0.5, n)
    }
  )
}
