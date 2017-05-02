#' Box Three-Dimensional Function
#'
#' Test function 12 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 3}, \code{m >= n}.
#'   \item Minima: \code{f = 0} at \code{(1, 10, 1), (10, 1, -1)} and
#'     where \code{x1 = x2} and \code{x3 = 0}.
#' }
#'
#' @param m Number of terms in the objective value. Should be 3 or greater.
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
#' Box, M. J. (1966).
#' A comparison of several current optimization methods, and the use of
#' transformations in constrained problems.
#' \emph{The Computer Journal}, \emph{9}(1), 67-77.
#' \url{https://doi.org/10.1093/comjnl/9.1.67}
#' @export
box_3d <- function(m = 15) {
  if (m < 3) {
    stop("box3d: m must be >= 3")
  }
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      for (i in 1:m) {
        ti <- 0.1 * i
        fi <- exp(-ti * x1) - exp(-ti * x2) - x3 * (exp(-ti) - exp(-i))
        fsum <- fsum + fi * fi
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      grad <- c(0, 0, 0)
      for (i in 1:m) {
        ti <- 0.1 * i
        etx1 <- exp(-ti * x1)
        etx2 <- exp(-ti * x2)
        et <- exp(-ti)
        fi <- etx1 - etx2 - x3 * (et - exp(-i))
        fi2 <- fi * 2
        grad[1] <- grad[1] - fi2 * ti * etx1
        grad[2] <- grad[2] + fi2 * ti * etx2
        grad[3] <- grad[3] - fi2 * (et - exp(-10 * ti))
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      fsum <- 0
      grad <- c(0, 0, 0)
      for (i in 1:m) {
        ti <- 0.1 * i
        etx1 <- exp(-ti * x1)
        etx2 <- exp(-ti * x2)
        et <- exp(-ti)
        fi <- etx1 - etx2 - x3 * (et - exp(-i))

        fsum <- fsum + fi * fi

        fi2 <- fi * 2
        grad[1] <- grad[1] - fi2 * ti * etx1
        grad[2] <- grad[2] + fi2 * ti * etx2
        grad[3] <- grad[3] - fi2 * (et - exp(-10 * ti))
      }
      grad

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(0, 10, 20)
  )
}