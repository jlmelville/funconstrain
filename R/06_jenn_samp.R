#' Jennrich and Sampson Function
#'
#' Test function 6 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 2}, number of summand
#'   functions \code{m >= n}.
#'   \item Minima: \code{f = 124.362...} at \code{(x1 = x2 = 0.2578)} for
#'   \code{m = 10},
#' }
#'
#' @note This test problem isn't really unconstrained. \code{x1} must take
#' a value between \code{(-1, 1)}. Included for the sake of completeness.
#'
#' @param m Number of summand functions in the objective function. Should be
#'   equal to or greater than 2.
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
#' Jennrich, R. I., & Sampson, P. F. (1968).
#' Application of stepwise regression to non-linear estimation.
#' \emph{Technometrics}, \emph{10}(1), 63-72.
#'
#' @examples
#' # Use m = 10 summand functions
#' fun <- jenn_samp(m = 10)
#' # Optimize using the standard starting point
#' # Set 'lower' and 'upper' parameter to constrain par[1]. Only works with
#' # L-BFGS-B.
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B", lower = -1, upper = 1)
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B",
#' lower = -1, upper = 1)
#'
#' # Use 20 summand functions
#' fun20 <- jenn_samp(m = 20)
#' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B",
#' lower = -1, upper = 1)
#' @export
jenn_samp <- function(m = 10) {
  if (m < 2) {
    stop("Jennrich-Sampson: m must be >= 2")
  }
  list(
    m = NA,
    fn = function(par) {
      x <- par[1]
      y <- par[2]

      fsum <- 0
      for (i in 1:m) {
        fi <- 2 + 2 * i - (exp(i * x) + exp(i * y))
        fsum <- fsum + fi * fi
      }
      fsum
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]

      grad <- c(0, 0)
      for (i in 1:m) {
        eix <- exp(i * x)
        eiy <- exp(i * y)
        fi <- 2 + 2 * i - (eix + eiy)
        grad[1] <- grad[1] - 2 * i * eix * fi
        grad[2] <- grad[2] - 2 * i * eiy * fi
      }
      grad
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      h <- matrix(0.0, nrow=2, ncol=2)
      for (i in 1:m) {
        d1 <- exp( i*x1 )
        d2 <- exp( i*x2 )
        t1 <- 2.0 + 2.0*i - ( d1 + d2 )
        h[1,1] <- h[1,1] + 2.0*( (i*d1) ^ 2 - t1*i ^ 2*d1 )
        h[1,2] <- h[1,2] + 2.0*i ^ 2*d1*d2
        h[2,2] <- h[2,2] + 2.0*( (i*d2) ^ 2 - t1*i ^ 2*d2 )
      }
      h[2,1] <- h[1,2]
      h
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]

      fsum <- 0
      grad <- c(0, 0)
      for (i in 1:m) {
        eix <- exp(i * x)
        eiy <- exp(i * y)
        fi <- 2 + 2 * i - (eix + eiy)
        fsum <- fsum + fi * fi
        grad[1] <- grad[1] - 2 * i * eix * fi
        grad[2] <- grad[2] - 2 * i * eiy * fi
      }
      grad

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0.3, 0.4),
    fmin = 124.362,
    xmin = c(0.2578, 0.2578) # for m = 10 (Caution!)
  )
}
