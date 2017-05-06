#' Gulf Research and Development Function
#'
#' Test function 11 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 3}, number of summand
#'   functions \code{3 <= m <= 100}.
#'   \item Minima: \code{f = 0} at \code{(50, 25, 1.5)}.
#' }
#'
#' Note that the equation as published by More', Garbow and Hillstrom (1981)
#' contains an error, where the symbol 'mi' should be interpreted as a minus
#' sign. The corrected version can be found in Jamil and Xang (2013), and no
#' doubt several other publications. The Jamil and Xang equation unfortunately
#' contains its own minor errors, but you can piece together the correct
#' equation from these two sources without too much trouble.
#'
#' @param m Number of summand functions in the objective function. Should be
#'   between 3 and 100, according to the MGH paper. Default value is 99, which
#'   Jamil and Xang (2013) list as the only valid value.
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
#' Jamil, M., & Yang, X. S. (2013).
#' A literature survey of benchmark functions for global optimisation problems.
#' \emph{International Journal of Mathematical Modelling and Numerical
#' Optimisation}, \emph{4}(2), 150-194.
#' \url{https://doi.org/10.1504/IJMMNO.2013.055204}
#' \url{https://arxiv.org/abs/1308.4008}
#'
#' @examples
#' # Use 10 summand functions
#' fun <- gulf(m = 10)
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
#'
#' # Use 20 summand functions
#' fun20 <- gulf(m = 20)
#' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
#' @export
gulf <- function(m = 99) {
  y <- c(34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030,
         6005, 5147, 4427, 3820, 3307, 2872)
  # can be between 3 and 100
  if (m < 3 || m > 100) {
    stop("Gulf research and development function: m must be between 3 and 100")
  }
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      p66 <- 2 / 3
      fsum <- 0
      for (i in 1:m) {
        ti <- i * 0.01
        y <- 25 + (-50 * log(ti)) ^ p66

        fi <- exp(-(abs(x2 - y) ^ x3) / x1) - ti
        fsum <- fsum + fi * fi
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      p66 <- 2 / 3
      grad <- c(0, 0, 0)
      for (i in 1:m) {
        ti <- i * 0.01
        y <- 25 + (-50 * log(ti)) ^ p66

        x2y <- x2 - y
        ax2y <- abs(x2y)
        x2yz <- ax2y ^ x3
        e <- exp(-x2yz / x1)
        fi <- e - ti

        grad[1] <- grad[1] + 2 * x2yz * e * fi / (x1 * x1)
        grad[2] <- grad[2] - 2 * x3 * e * x2yz * fi / (x1 * x2y)
        grad[3] <- grad[3] - 2 * x2yz * log(ax2y) * e * fi / x1
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      p66 <- 2 / 3
      fsum <- 0
      grad <- c(0, 0, 0)
      for (i in 1:m) {
        ti <- i * 0.01
        y <- 25 + (-50 * log(ti)) ^ p66

        x2y <- x2 - y
        ax2y <- abs(x2y)
        x2yz <- ax2y ^ x3
        e <- exp(-x2yz / x1)
        fi <- e - ti

        fsum <- fsum + fi * fi

        grad[1] <- grad[1] + 2 * x2yz * e * fi / (x1 * x1)
        grad[2] <- grad[2] - 2 * x3 * e * x2yz * fi / (x1 * x2y)
        grad[3] <- grad[3] - 2 * x2yz * log(ax2y) * e * fi / x1
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(5, 2.5, 0.15)
  )
}
