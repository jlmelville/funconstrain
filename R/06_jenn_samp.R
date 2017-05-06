#' Jennrich and Sampson Function
#'
#' Test function 6 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 2}, \code{m >= n}.
#'   \item Minima: \code{f = 124.362...} at \code{(x1 = x2 = 0.2578)} for
#'   \code{m = 10},
#' }
#'
#' @param m Number of terms in the objective function.
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Function which calculates the gradient vector given input
#'   parameter vector.
#'   \item \code{fg} A function which calculates both the value and gradient,
#'   given the parameter vector, returning a list.
#'   \item \code{x0} Standard starting point.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Jennrich, R. I., & Sampson, P. F. (1968).
#' Application of stepwise regression to non-linear estimation.
#' \emph{Technometrics}, \emph{10}(1), 63-72.
#' @export
jenn_samp <- function(m = 10) {
  list(
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
    x0 = c(0.3, 0.4)
  )
}
