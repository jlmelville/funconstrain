#' Broyden Tridiagonal Function
#'
#' Test function 30 from the Moré, Garbow and Hillstrom paper.
#'
#' The objective function is the sum of `m` functions, each of `n`
#' parameters.
#'
#' - Dimensions: Number of parameters `n` variable, number of summand
#'   functions `m = n`.
#' - Minima: `f = 0`.
#'
#' The number of parameters, `n`, in the objective function is not
#' specified when invoking this function. It is implicitly set by the length of
#' the parameter vector passed to the objective and gradient functions that this
#' function creates. See the 'Examples' section.
#'
#' @template factory-return
#' @references
#' Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' *ACM Transactions on Mathematical Software (TOMS)*, *7*(1), 17-41.
#' <https://doi.org/10.1145/355934.355936>
#'
#' Broyden, C. G. (1965).
#' A class of methods for solving nonlinear simultaneous equations.
#' *Mathematics of computation*, *19*(92), 577-593.
#' <https://doi.org/10.2307/2003941>
#'
#' @examples
#' btri <- broyden_tri()
#' # 6 variable problem using the standard starting point
#' x0_6 <- btri$x0(6)
#' res_6 <- stats::optim(x0_6, btri$fn, btri$gr, method = "L-BFGS-B")
#' # Standard starting point with 8 variables
#' res_8 <- stats::optim(btri$x0(8), btri$fn, btri$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), btri$fn, btri$gr,
#'                       method = "L-BFGS-B")
#' @export
broyden_tri <- function() {
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]

      grad <- 2 * fi * (3 - 4 * par)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 4 * fi[1:(n - 1)]

      grad
    },
    he = function(x) {
      n <- length(x)
      h <- matrix(0.0, nrow = n, ncol = n)
      #       ! For i <- 1
      t <- (3.0 - 2.0 * x[1]) * x[1] - 2.0 * x[2] + 1.0
      t1 <- 3.0 - 4.0 * x[1]
      h[1, 1] <- 2.0 * (t1^2 - 4.0 * t)
      h[1, 2] <- -4.0 * t1
      h[2, 2] <- 8.0

      if (n > 2) {
        for (i in 2:(n - 1)) {
          t <- (3.0 - 2.0 * x[i]) * x[i] - x[i - 1] - 2.0 * x[i + 1] + 1.0
          t1 <- 3.0 - 4.0 * x[i]
          h[i - 1, i - 1] <- h[i - 1, i - 1] + 2.0
          h[i - 1, i] <- h[i - 1, i] - 2.0 * t1
          h[i, i] <- h[i, i] + 2.0 * (t1^2 - 4.0 * t)
          h[i - 1, i + 1] <- h[i - 1, i + 1] + 4.0
          h[i, i + 1] <- h[i, i + 1] - 4.0 * t1
          h[i + 1, i + 1] <- h[i + 1, i + 1] + 8.0
        }
      }

      #       ! For i <- n
      t <- (3.0 - 2.0 * x[n]) * x[n] - x[n - 1] + 1.0
      t1 <- 3.0 - 4.0 * x[n]
      h[n - 1, n - 1] <- h[n - 1, n - 1] + 2.0
      h[n - 1, n] <- h[n - 1, n] - 2.0 * t1
      h[n, n] <- h[n, n] + 2.0 * (t1^2 - 4.0 * t)

      for (j in 1:(n - 1)) {
        # symmetrize
        for (k in (j + 1):n) {
          h[k, j] <- h[j, k]
        }
      }
      h
    },

    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      fsum <- sum(fi * fi)

      grad <- 2 * fi * (3 - 4 * par)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 4 * fi[1:(n - 1)]

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 40) {
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }
      rep(-1, n)
    },
    fmin = 0,
    xmin = c(-0.5648284, -0.6662737, -0.6609170, -0.5950500, -0.4162011) # n = 5 case
  )
}
