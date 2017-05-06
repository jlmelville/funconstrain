#' Discrete Integral Equation Function
#'
#' Test function 29 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable, \code{m = n}.
#'   \item Minima: \code{f = 0} (at the same location as
#'   \code{\link{disc_bv}}).
#' }
#'
#' The number of variables in the function, \code{n}, is determined by the
#' length of the vector passed to the function and gradient routines. See
#' the 'Examples' section.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Objective function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Gradient function which calculates the gradient vector
#'   given input parameter vector.
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
#' More', J. J., & Cosnard, M. Y. (1979).
#' Numerical solution of nonlinear equations.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{5}(1), 64-85.
#' \url{https://doi.org/10.1145/355815.355820}
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
    x0 = function(n) {
      if (n < 1) {
        stop("Discrete Integral Equation: n must be positive")
      }
      (1:n / n + 1) * ((1:n / n + 1) - 1)
    }
  )
}
