#' Biggs EXP6 Function
#'
#' Test function 18 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 6}, \code{m >= n}.
#'   \item Minima: \code{f = 5.65565...e-3} if \code{m = 13};
#'   not reported in the MGH (1981) paper is \code{(f = 0)} at
#'   \code{c(1, 10, 1, 5, 4, 3)} and \code{c(4, 10, 3, 5, 1, 1)}
#'   (and probably others) for all \code{m} (probably: I stopped testing after
#'   \code{m = 1000}).
#' }
#'
#' @param m Number of terms in the objective value.
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
#' Biggs, M. C. (1971).
#' Minimization algorithms making use of non-quadratic properties of the
#' objective function.
#' \emph{IMA Journal of Applied Mathematics}, \emph{8}(3), 315-327.
#' @export
biggs_exp6 <- function(m = 13) {
  if (m < 6) {
    stop("Biggs EXP6 function must have m >= 6")
  }

  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      fsum <- 0
      for (i in 1:m) {
        ti <- 0.1 * i
        y <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
        f <- x3 * exp(-ti * x1) - x4 * exp(-ti * x2) + x6 * exp(-ti * x5) - y
        fsum <- fsum + f * f
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      grad <- c(0, 0, 0, 0, 0, 0)
      for (i in 1:m) {
        ti <- 0.1 * i
        y <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
        f1 <- exp(-ti * x1)
        f31 <- x3 * f1
        f2 <- exp(-ti * x2)
        f42 <- x4 * f2
        f5 <- exp(-ti * x5)
        f65 <- x6 * exp(-ti * x5)
        f <- f31 - f42 + f65 - y

        grad[1] <- grad[1] - 2 * ti * f31 * f
        grad[2] <- grad[2] + 2 * ti * f42 * f
        grad[3] <- grad[3] + 2 * f1 * f
        grad[4] <- grad[4] - 2 * f2 * f
        grad[5] <- grad[5] - 2 * ti * f65 * f
        grad[6] <- grad[6] + 2 * f5 * f
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      fsum <- 0
      grad <- c(0, 0, 0, 0, 0, 0)
      for (i in 1:m) {
        ti <- 0.1 * i
        y <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
        f1 <- exp(-ti * x1)
        f31 <- x3 * f1
        f2 <- exp(-ti * x2)
        f42 <- x4 * f2
        f5 <- exp(-ti * x5)
        f65 <- x6 * exp(-ti * x5)
        f <- f31 - f42 + f65 - y

        fsum <- fsum + f * f

        grad[1] <- grad[1] - 2 * ti * f31 * f
        grad[2] <- grad[2] + 2 * ti * f42 * f
        grad[3] <- grad[3] + 2 * f1 * f
        grad[4] <- grad[4] - 2 * f2 * f
        grad[5] <- grad[5] - 2 * ti * f65 * f
        grad[6] <- grad[6] + 2 * f5 * f
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(1, 2, 1, 1, 1, 1)
  )
}
