#' Penalty Function I
#'
#' Test function 23 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable, \code{m = n + 1}.
#'   \item Minima: \code{f = 2.24997...e-5} if \code{n = 4};
#'   \code{f = 7.08765...e-5} if \code{n = 10}.
#' }
#'
#' The number of variables in the function, \code{n}, is determined by the
#' length of the vector passed to the function and gradient routines. See
#' the 'Examples' section.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Function which calculates the gradient vector given input
#'   parameter vector.
#'   \item \code{fg} Function which calculates both the value and gradient,
#'   given the parameter vector, returning a list.
#'   \item \code{x0} Function returning the standard starting point, given
#'   \code{n}, the number of variables desired.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Gill, P. E., Murray, W., & Pitfield, R. A. (1972).
#' \emph{The implementation of two revised quasi-Newton algorithms for
#' unconstrained optimization} (Report NAC-11).
#' Teddington, London: National Physical Laboratory, Division of Numerical
#' Analysis and Computing.
#' @examples
#' pen1 <- penalty_1()
#' # 6 variable problem using the standard starting point
#' x0_6 <- pen1$x0(6)
#' res_6 <- stats::optim(x0_6, pen1$fn, pen1$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(pen1$x0(8), pen1$fn, pen1$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), pen1$fn, pen1$gr,
#'                       method = "L-BFGS-B")
#' @export
penalty_1 <- function() {
  a <- 1e-5
  sqrta <- sqrt(a)
  list(
    fn = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }
      fsum <- 0
      fn1 <- 0
      for (i in 1:n) {
        fi <- sqrta * (par[i] - 1)
        fsum <- fsum + fi * fi
        fn1 <- fn1 + par[i] * par[i]
      }

      fn1 <- fn1 - 0.25
      fsum <- fsum + fn1 * fn1
      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }
      grad <- rep(0, n)
      fn1 <- 0

      for (i in 1:n) {
        fi <- sqrta * (par[i] - 1)
        grad[i] <- grad[i] + 2 * sqrta * fi
        fn1 <- fn1 + par[i] * par[i]
      }
      fn1 <- fn1 - 0.25
      grad <- grad + 4 * par * fn1
      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }

      fn1 <- 0
      grad <- rep(0, n)
      fsum <- 0
      for (i in 1:n) {
        fi <- sqrta * (par[i] - 1)
        fsum <- fsum + fi * fi
        grad[i] <- grad[i] + 2 * sqrta * fi
        fn1 <- fn1 + par[i] * par[i]
      }
      fn1 <- fn1 - 0.25
      fsum <- fsum + fn1 * fn1
      grad <- grad + 4 * par * fn1

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 10) {
      if (n < 1) {
        stop("Penalty Function I: n must be positive")
      }
      1:n
    }
  )
}