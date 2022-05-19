#' Chebyquad Function
#'
#' Test function 35 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, number of
#'   summand functions \code{m = n}.
#'   \item Minima: \code{f = 0} for \code{m = n}, \code{1 <= n <= 7} and \code{n
#'   = 9};
#'   \code{f = 3.51687...e-3} for \code{m = n = 8};
#'   \code{f = 6.50395...e-3} for \code{m = n = 10}.
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
#' Fletcher, R. (1965).
#' Function minimization without evaluating derivatives - a review.
#' \emph{The Computer Journal}, \emph{8}(1), 33-41.
#' \url{https://doi.org/10.1093/comjnl/8.1.33}
#'
#' @examples
#' cheb <- chebyquad()
#' # 6 variable problem using the standard starting point
#' x0_6 <- cheb$x0(n = 6)
#' res_6 <- stats::optim(x0_6, cheb$fn, cheb$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(cheb$x0(8), cheb$fn, cheb$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), cheb$fn, cheb$gr, method =
#' "L-BFGS-B")
#' @export
chebyquad <- function() {

  list(
    m = NA,
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Chebyquad: n must be positive")
      }

      # y is the shifted x
      y <- 2 * par - 1

      t0 <- rep(1, n)
      t1 <- y
      ti <- t1
      fsum <- 0
      for (i in 1:n) {
        if (i > 1) {
          ti <- 2 * y * t1 - t0
          t0 <- t1
          t1 <- ti
        }

        fi <- sum(ti) / n
        if (i %% 2 == 0) {
          fi <- fi + 1 / (i * i - 1)
        }
        fsum <- fsum + fi * fi
      }

      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Chebyquad: n must be positive")
      }

      y <- 2 * par - 1

      g0 <- rep(0, n)
      # grad of T1 wrt to shift is 2, but we'll account for that later
      g1 <- rep(1, n)

      t0 <- rep(1, n)
      t1 <- y
      grad <- rep(0, n)

      gi <- g1
      ti <- t1
      for (i in 1:n) {
        if (i > 1) {
          # Gn = 2T_old + 2xG_old - G_old_old
          gi <- 2 * t1 + 2 * y * g1 - g0
          g0 <- g1
          g1 <- gi

          # Update T
          ti <- 2 * y * t1 - t0
          t0 <- t1
          t1 <- ti
        }

        # Correct for the fact we want want deriv wrt x, not y
        # G wrt x is 2 * G wrt shifted (y)
        gi <- 2 * gi / n

        fi <- sum(ti) / n
        if (i %% 2 == 0) {
          fi <- fi + 1 / (i * i - 1)
        }

        grad <- grad + 2 * fi * gi
      }
      grad
    },
    he = function(x) { # ?? seems to have error in 1,1 element
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)
       fvec <- rep(0.0,n)
       gvec <- rep(0.0,n)
       for (j in 1:n) {
          t1 <- 1.0
          t2 <- 2.0*x[j] - 1.0
          t <- 2.0*t2
          for (i in 1:n){ 
             fvec[i] <- fvec[i] + t2
             th <- t*t2 - t1
             t1 <- t2
             t2 <- th
          }
       }
       d1 <- 1.0/n
       for (i in 1:n) {
          fvec[i] <- d1*fvec[i]
          if ( (i %% 2) == 0 ) {
             fvec[i] <- fvec[i] + 1.0/( i^2 - 1.0 )
          }
       }
       d2 <- 2.0*d1
       for (j in 1:n) { 
          h[j,j] <- 4.0*d1
          t1 <- 1.0
          t2 <- 2.0*x[j] - 1.0
          t <- 2.0*t2
          s1 <- 0.0
          s2 <- 2.0
          p1 <- 0.0
          p2 <- 0.0
          gvec[1] <- s2
          for (i in 2:n) {
             th <- 4.0*t2 + t*s2 - s1
             s1 <- s2
             s2 <- th
             th <- t*t2 - t1
             t1 <- t2
             t2 <- th
             th <- 8.0*s1 + t*p2 - p1
             p1 <- p2
             p2 <- th
             gvec[i] <- s2
             h[j,j] <- h[j,j] + fvec[i]*th + d1*s2 ^ 2
          }
          h[j,j] <- d2*h[j,j]
          if (j > 1) { # needed for R
           for (k in 1:(j-1)) {
             h[k,j] <- 0.0
             tt1 <- 1.0
             tt2 <- 2.0*x[k] - 1.0
             tt <- 2.0*tt2
             ss1 <- 0.0
             ss2 <- 2.0
             for (i in 1:n) {
                h[k,j] <- h[k,j] + ss2*gvec[i]
                tth <- 4.0*tt2 + tt*ss2 - ss1
                ss1 <- ss2
                ss2 <- tth
                tth <- tt*tt2 - tt1
                tt1 <- tt2
                tt2 <- tth
             }
             h[k,j] <- d2*d1*h[k,j]
           }
          } # end if
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
        stop("Chebyquad: n must be positive")
      }

      y <- 2 * par - 1

      g0 <- rep(0, n)
      # grad of T1 wrt to shift is 2, but we'll account for that later
      g1 <- rep(1, n)

      t0 <- rep(1, n)
      t1 <- y

      fsum <- 0
      grad <- rep(0, n)

      gi <- g1
      ti <- t1
      for (i in 1:n) {
        if (i > 1) {
          # Gn = 2T_old + 2xG_old - G_old_old
          gi <- 2 * t1 + 2 * y * g1 - g0
          g0 <- g1
          g1 <- gi

          # Update T
          ti <- 2 * y * t1 - t0
          t0 <- t1
          t1 <- ti
        }

        # Correct for the fact we want want deriv wrt x, not y
        # G wrt x is 2 * G wrt shifted (y)
        gi <- 2 * gi / n

        fi <- sum(ti) / n
        if (i %% 2 == 0) {
          fi <- fi + 1 / (i * i - 1)
        }
        fsum <- fsum + fi * fi
        grad <- grad + 2 * fi * gi
      }
      grad

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 50) {
      if (n < 1) {
        stop("Chebyquad: n must be positive")
      }

      1:n / (n + 1)
    }
  )
}
