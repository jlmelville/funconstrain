#' Freudenstein and Roth Function
#'
#' Test function 2 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 2}, \code{m = 2}.
#'   \item Minima: \code{f = 0} at \code{(5, 4)},
#'   \code{f = 48.9842...} at \code{(11.41..., -0.8968...)}
#' }
#'
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
#' Freudenstein, F., & Roth, B. (1963).
#' Numerical solution of systems of nonlinear equations.
#' \emph{Journal of the ACM (JACM)}, \emph{10}(4), 550-556.
#' \url{https://doi.org/10.1145/321186.321200}
#' @export
freud_roth <- function() {
  list(
    fn = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- -13 + x + ((5 - y) * y - 2) * y
      f2 <- -29 + x + ((y + 1) * y - 14) * y

      f1 * f1 + f2 * f2
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]
      yy <- y * y
      yyyy <- yy * yy
      c(
        4 * x + 12 * y * y - 32 * y - 84,
        12 * yyyy * y - 40 * yyyy + 8 * yy * y - 240 * yy +
          (24 * x + 24) * y - 32 * x + 864
      )
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]
      f1 <- -13 + x + ((5 - y) * y - 2) * y
      f2 <- -29 + x + ((y + 1) * y - 14) * y

      yy <- y * y
      yyyy <- yy * yy

      list(
        fn = f1 * f1 + f2 * f2,
        gr = c(2 * f1 + 2 * f2,
               12 * yyyy * y - 40 * yyyy + 8 * yy * y - 240 * yy +
                 (24 * x + 24) * y - 32 * x + 864)
      )
    },
    x0 = c(0.5, -2)
  )
}
