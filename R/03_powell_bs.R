#' Powell Badly Scaled Function
#'
#' Test function 3 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 2}, number of summand
#'   functions \code{m = 2}.
#'   \item Minima: \code{f = 0} at \code{(1.098... e-5, 9.106...) },
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Objective function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Gradient function which calculates the gradient vector
#'   given input parameter vector.
#'   \item \code{he} If available, the hessian matrix (second derivatives)
#'   of the function w.r.t. the parameters at the given values.
#'   \item \code{fg} A function which, given the parameter vector, calculates
#'   both the objective value and gradient, returning a list with members
#'   \code{fn} and \code{gr}, respectively.
#'   \item \code{x0} Standard starting point.
#'   \item \code{fmin} reported minimum
#'   \item \code{xmin} parameters at reported minimum
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \doi{doi.org/10.1145/355934.355936}
#'
#' Powell, M. J. D. (1970).
#' A hybrid method for nonlinear equations.
#' In P. Rabinowitz (Ed.),
#' \emph{Numerical Methods for Nonlinear Algebraic Equations} (pp87-114).
#' New York: Gordon & Breach.
#'
#' @examples
#' fun <- powell_bs()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
powell_bs <- function() {
  list(
    m = NA,
    fn = function(par) {
      x <- par[1]
      y <- par[2]

      f1 <- 1e4 * x * y - 1
      f2 <- exp(-x) + exp(-y) - 1.0001

      f1 * f1 + f2 * f2
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
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
       h <- matrix(NA, nrow=2, ncol=2)
       t1 <- 1.0e+04*x1*x2 - 1
       t2 <- exp( - x1 ) + exp( - x2 ) - 1.0001
       h[1,1] <- 2.0*( 1.0e+8*x2 ^ 2 + t2*exp( - x1 ) + exp( - x1 ) ^ 2 )
       h[1,2] <- 2.0*( t1*1.0e+04 + 1.0e+8*x1*x2 + exp( - x1 - x2 ) )
       h[2,2] <- 2.0*( 1.0e+8*x1 ^ 2 + t2*exp( - x2 ) + exp( - x2 ) ^ 2 )
       h[2,1] <- h[1,2]
       h
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
    x0 = c(0, 1),
    fmin = 0,
    xmin = c(1.098e-5, 9.106)
   )
}
