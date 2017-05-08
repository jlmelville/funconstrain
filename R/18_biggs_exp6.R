#' Biggs EXP6 Function
#'
#' Test function 18 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 6}, number of summand
#'   functions \code{m >= n}.
#'   \item Minima: \code{f = 5.65565...e-3} if \code{m = 13};
#'   not reported in the MGH (1981) paper is \code{(f = 0)} at
#'   \code{c(1, 10, 1, 5, 4, 3)} and \code{c(4, 10, 3, 5, 1, 1)}
#'   (and probably others) for all \code{m} (probably: I stopped testing after
#'   \code{m = 1000}).
#' }
#'
#' @param m Number of summand functions in the objective function. Should be
#'   equal to or greater than 6.
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Objective function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Gradient function which calculates the gradient vector
#'   given input parameter vector.
#'   \item \code{fg} A function which, given the parameter vector, calculates
#'   both the objective value and gradient, returning a list with members
#'   \code{fn} and \code{gr}, respectively.
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
#'
#' @examples
#' fun <- biggs_exp6()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(1:6, fun$fn, fun$gr, method = "L-BFGS-B")
#'
#' # Use 20 summand functions
#' fun20 <- biggs_exp6(m = 20)
#' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
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

      ti <- 0.1 * (1:m)
      yi <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
      fi <- x3 * exp(-ti * x1) - x4 * exp(-ti * x2) + x6 * exp(-ti * x5) - yi
      sum(fi * fi)

    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      ti <- 0.1 * (1:m)
      yi <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
      f1 <- exp(-ti * x1)
      f31 <- x3 * f1
      f2 <- exp(-ti * x2)
      f42 <- x4 * f2
      f5 <- exp(-ti * x5)
      f65 <- x6 * exp(-ti * x5)
      fi <- f31 - f42 + f65 - yi
      fi2 <- 2 * fi
      tfi2 <- fi2 * ti

      c(
        -sum(tfi2 * f31),
         sum(tfi2 * f42),
         sum(fi2 * f1),
        -sum(fi2 * f2),
        -sum(tfi2 * f65),
         sum(fi2 * f5)
      )
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      ti <- 0.1 * (1:m)
      yi <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
      f1 <- exp(-ti * x1)
      f31 <- x3 * f1
      f2 <- exp(-ti * x2)
      f42 <- x4 * f2
      f5 <- exp(-ti * x5)
      f65 <- x6 * exp(-ti * x5)
      fi <- f31 - f42 + f65 - yi
      fi2 <- 2 * fi
      tfi2 <- fi2 * ti

      fsum <- sum(fi * fi)
      grad <- c(
        -sum(tfi2 * f31),
        sum(tfi2 * f42),
        sum(fi2 * f1),
        -sum(fi2 * f2),
        -sum(tfi2 * f65),
        sum(fi2 * f5)
      )

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(1, 2, 1, 1, 1, 1)
  )
}
