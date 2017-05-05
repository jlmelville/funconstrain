#' Penalty Function II
#'
#' Test function 24 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable, \code{m = n + 1}.
#'   \item Minima: \code{f = 9.37629...e-6} if \code{n = 4};
#'   \code{f = 2.93660...e-4} if \code{n = 10}.
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
    x0 = function(n = 10) {
      if (n < 1) {
        stop("Penalty Function II: n must be positive")
      }
      rep(0.5, n)
    }
  )
}
