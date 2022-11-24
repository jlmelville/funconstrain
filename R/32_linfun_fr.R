#' Linear Function - Full Rank
#'
#' Test function 32 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n} variable, \code{m >= n}.
#'   \item Minima: \code{f = m - n} at \code{rep(-1, n)}.
#' }
#'
#' The number of parameters, \code{n}, in the objective function is not
#' specified when invoking this function. It is implicitly set by the length of
#' the parameter vector passed to the objective and gradient functions that this
#' function creates. See the 'Examples' section.
#'
#' @param m Number of summand functions in the objective function. Should be
#'   equal to or greater than \code{n}.
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
#' @examples
#' linfr <- linfun_fr(m = 10)
#' # 6 variable problem using the standard starting point
#' x0_6 <- linfr$x0(n = 6)
#' res_6 <- stats::optim(x0_6, linfr$fn, linfr$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(linfr$x0(8), linfr$fn, linfr$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), linfr$fn, linfr$gr,
#'                       method = "L-BFGS-B")
#'
#' # Use 20 summand functions
#' linfr_m20 <- linfun_fr(m = 20)
#' # Repeat 4 parameter optimization with new test function
#' res_n4_m20 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), linfr_m20$fn, linfr_m20$gr,
#' method = "L-BFGS-B")
#' @export
linfun_fr <- function(m = 100) {
  if (m < 1) {
    stop("Linear Function - Full Rank: m must be positive")
  }

  m2 <- 2 / m
  m4 <- 2 * m2

  list(
    m = m,
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Linear Function - Full Rank: n must be positive")
      }
      if (m < n) {
        stop("Linear Function - Full Rank: m must be >= n")
      }
      sum_x <- sum(par)
      msx <- m2 * sum_x + 1
      fi <- par - msx
      fnm <- -msx

      sum(fi * fi) + (m - n) * fnm * fnm
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Linear Function - Full Rank: n must be positive")
      }
      if (m < n) {
        stop("Linear Function - Full Rank: m must be >= n")
      }
      sum_x <- sum(par)
      msx <- m2 * sum_x + 1
      fi <- par - msx
      fnm <- -msx

      2 * fi - (m4 * (sum(fi) + (m - n) * fnm))
    },
    he = function(x) {
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)
       for (i in 1:n){
         h[i, i] <- 2.0 # since quadratic
       }
       h
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Linear Function - Full Rank: n must be positive")
      }
      if (m < n) {
        stop("Linear Function - Full Rank: m must be >= n")
      }
      sum_x <- sum(par)
      msx <- m2 * sum_x + 1
      fi <- par - msx
      fnm <- -msx

      fsum <- sum(fi * fi) + (m - n) * fnm * fnm
      grad <- 2 * fi - (m4 * (sum(fi) + (m - n) * fnm))

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 45) {
      if (n < 1) {
        stop("Linear Function - Full Rank: n must be positive")
      }
      if (m < n) {
        stop("Linear Function - Full Rank: m must be >= n")
      }
      rep(1, n)
    }
  )
}
