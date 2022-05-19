#' Powell Singular Function
#'
#' Test function 13 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 4}, number of summand
#'   functions \code{m = 4}.
#'   \item Minima: \code{f = 0} at \code{rep(0, 4)}.
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
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Powell, M. J. D. (1962).
#' An iterative method for finding stationary values of a function of several
#' variables.
#' \emph{The Computer Journal}, \emph{5}(2), 147-151.
#' \url{https://doi.org/10.1093/comjnl/5.2.147}
#'
#' @examples
#' fun <- powell_s()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
#' "L-BFGS-B")
#' @export
powell_s <- function() {

  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      x14 <- x1 - x4
      x14s <- x14 * x14
      f1 <- x1 + 10 * x2
      f2s <- 5 * (x3 - x4) ^ 2
      f3 <- (x2 - 2 * x3) ^ 2
      f4s <- 10 * x14s * x14s

      f1 * f1 + f2s + f3 * f3 + f4s
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      x12 <- x1 + 10 * x2
      x14 <- x1 - x4
      x14_3 <- x14 * x14 * x14
      x34 <- x3 - x4
      x23 <- x2 - 2 * x3
      x23_3 <- x23 * x23 * x23

      c(
        40 * x14_3 + 2 * x12,
        4 * x23_3 + 20 * x12,
        10 * x34 - 8 * x23_3,
        -10 * x34 - 40 * x14_3
      )
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      h <- matrix(0.0, ncol=4, nrow=4)
      d3 <- ( x2 - 2.0*x3 ) ^ 2
      d4 <- sqrt( 10.0 )*( x1 - x4 ) ^ 2
      t1 <- 2.0*( x2 - 2.0*x3 )
      t2 <- 2.0*sqrt( 10.0 )*( x1 - x4 )
      t3 <- 2.0*sqrt( 10.0 )
      h[1,1] <- 2.0*( 1.0 + d4*t3 + t2 ^ 2 )
      h[1,2] <- 20.0
      h[2,2] <- 2.0*( 100.0 + 2.0*d3 + t1 ^ 2 )
      h[1,3] <- 0.0
      h[2,3] <- - 4.0*( 2.0*d3 + t1 ^ 2 )
      h[3,3] <- 8.0*( 1.25 + 2.0*d3 + t1 ^ 2 )
      h[1,4] <- - 2.0*( d4*t3 + t2 ^ 2 )
      h[2,4] <- 0.0
      h[3,4] <- - 10.0
      h[4,4] <- 2.0*( 5.0 + d4*t3 + t2 ^ 2 )
      h[2,1] <- h[1,2]
      h[3,1] <- h[1,3]
      h[3,2] <- h[2,3]
      h[4,1] <- h[1,4]
      h[4,2] <- h[2,4]
      h[4,3] <- h[3,4]
      h
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      x14 <- x1 - x4
      x14_2 <- x14 * x14
      x14_3 <- x14_2 * x14
      x34 <- x3 - x4
      x23 <- x2 - 2 * x3
      x23_2 <- x23 * x23
      x23_3 <- x23_2 * x23
      x12 <- x1 + 10 * x2

      f1 <- x12
      f2s <- 5 * x34 * x34
      f3 <- x23_2
      f4s <- 10 * x14_2 * x14_2

      list(
        fn = f1 * f1 + f2s + f3 * f3 + f4s,
        gr = c(
          40 * x14_3 + 2 * x12,
          4 * x23_3 + 20 * x12,
          10 * x34 - 8 * x23_3,
          -10 * x34 - 40 * x14_3
        )
      )
    },
    x0 = c(3, -1, 0, 1)
  )
}
