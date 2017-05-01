#' Powell Badly Scaled Function
#'
#' Test function 3 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 2}, \code{m = 2}.
#'   \item Minima: \code{f = 0} at \code{(1.098... e-5, 9.106...) },
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} The function.
#'   \item \code{gr} The gradient.
#'   \item \code{x0} Standard starting point.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Powell, M. J. D. (1970).
#' A hybrid method for nonlinear equations.
#' In P. Rabinowitz (Ed.),
#' \emph{Numerical Methods for Nonlinear Algebraic Equations} (pp87-114).
#' New York: Gordon & Breach.
#' @export
powell_bs <- function() {
  list(
    fn = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- 1e4 * x * y - 1
      f2 <- exp(-x) + exp(-y) - 1.0001

      f1 ^ 2 + f2 ^ 2
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]
      a <- 2e4 * (1e4 * x * y - 1)
      b <- 2 * (exp(-x) + exp(-y) - 1.0001)

      c(
          y * a - exp(-x) * b,
          x * a - exp(-y) * b
      )
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- 1e4 * x * y - 1
      f2 <- exp(-x) + exp(-y) - 1.0001

      a <- 2e4 * f1
      b <- 2 * f2

      list(
        fn = f1 * f1 + f2 * f2,
        gr = c(
          y * a - exp(-x) * b,
          x * a - exp(-y) * b
        )
      )
    },
    x0 = c(0, 1)
  )
}
