#' Broyden Banded Function
#'
#' Test function 31 from the Moré, Garbow and Hillstrom paper.
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
#' Broyden, C. G. (1971).
#' The convergence of an algorithm for solving sparse nonlinear systems.
#' *Mathematics of Computation*, *25*(114), 285-294.
#' <https://doi.org/10.1090/S0025-5718-1971-0297122-5>
#'
#' @examples
#' btri <- broyden_band()
#' # 6 variable problem using the standard starting point
#' x0_6 <- btri$x0(6)
#' res_6 <- stats::optim(x0_6, btri$fn, btri$gr, method = "L-BFGS-B")
#' # Standard starting point with 8 variables
#' res_8 <- stats::optim(btri$x0(8), btri$fn, btri$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), btri$fn, btri$gr,
#'                       method = "L-BFGS-B")
#' @export
broyden_band <- function() {
  ml <- 5
  mu <- 1
  # return the valid indexes in a band around i
  #  the range i-mlow:i+mupp inclusive, excluding i, and j < 1 or j > n
  band_j <- function(i, n, mlow = ml, mupp = mu) {
    j <- max(1, i - mlow):min(n, i + mupp)
    j[j != i]
  }

  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }

      xx <- par * par
      fi <- (2 + 5 * xx) * par + 1
      fj <- par + xx
      for (i in 1:n) {
        fi[i] <- fi[i] - sum(fj[band_j(i, n)])
      }
      sum(fi * fi)
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }

      xx <- par * par
      fi <- (2 + 5 * xx) * par + 1
      fj <- par + xx
      for (i in 1:n) {
        fi[i] <- fi[i] - sum(fj[band_j(i, n)])
      }

      grad <- 2 * fi * (15 * xx + 2)
      gj <- 2 * par + 1
      for (i in 1:n) {
        # reversing the limits returns the bands that i is part of
        grad[i] <- grad[i] - 2 * gj[i] * sum(fi[band_j(i, n, mu, ml)])
      }

      grad
    },
    he = function(x) {
      n <- length(x)
      h <- matrix(0.0, nrow = n, ncol = n)

      for (i in 1:n) {
        s1 <- 0.0
        if (i > 1) {
          # needed for R
          for (j in max(1, i - 5):(i - 1)) {
            s1 <- s1 + x[j] * (1.0 + x[j])
          }
        }
        if (i < n) {
          s1 <- s1 + x[i + 1] * (1.0 + x[i + 1])
        }

        t <- x[i] * (2.0 + 5.0 * x[i]^2) + 1.0 - s1
        d1 <- 2.0 + 15.0 * x[i]^2

        if (i > 1) {
          # needed for R
          for (j in (max(1, i - 5):(i - 1))) {
            d2 <- -1.0 - 2.0 * x[j]
            h[j, j] <- h[j, j] + 2.0 * (d2^2 - 2.0 * t)
            if ((j + 1) < i) {
              # needed for R
              for (l in ((j + 1):(i - 1))) {
                h[j, l] <- h[j, l] + 2.0 * d2 * (-1.0 - 2.0 * x[l])
              }
            }
            h[j, i] <- h[j, i] + 2.0 * d1 * d2
            if (i < n) {
              h[j, i + 1] <- h[j, i + 1] + 2.0 * d2 * (-1.0 - 2.0 * x[i + 1])
            }
          }
        } # end if i>1

        h[i, i] <- h[i, i] + 2.0 * (30.0 * t * x[i] + d1^2)
        if (i < n) {
          d2 <- -1.0 - 2.0 * x[i + 1]
          h[i, i + 1] <- h[i, i + 1] + 2.0 * d1 * d2
          h[i + 1, i + 1] <- h[i + 1, i + 1] + 2.0 * (d2^2 - 2.0 * t)
        }
      }

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
        stop("Broyden Banded: n must be positive")
      }

      xx <- par * par
      fi <- (2 + 5 * xx) * par + 1
      fj <- par + xx
      for (i in 1:n) {
        fi[i] <- fi[i] - sum(fj[band_j(i, n)])
      }
      fsum <- sum(fi * fi)

      grad <- 2 * fi * (15 * xx + 2)
      gj <- 2 * par + 1
      for (i in 1:n) {
        grad[i] <- grad[i] - 2 * gj[i] * sum(fi[band_j(i, n, mu, ml)])
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 40) {
      if (n < 1) {
        stop("Broyden Banded: n must be positive")
      }
      rep(-1, n)
    },
    fmin = 0,
    xmin = c(-0.4283029, -0.4765965, -0.5196377, -0.5588620, -0.5588620) # n = 5 case
  )
}
