#' Beale Function
#'
#' Test function 5 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 2}, number of summand
#'   functions \code{m = 3}.
#'   \item Minima: \code{f = 0} at \code{(3, 0.5) },
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
#' \doi{doi.org/10.1145/355934.355936}
#'
#' Beale, E. M. L. (1958).
#' \emph{On an iterative method for finding a local minimum of a function of more
#' than one variable} (No. 25).
#' Statistical Techniques Research Group, Section of Mathematical Statistics,
#' Department of Mathematics, Princeton University.
#'
#' @examples
#' fun <- beale()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
beale <- function() {
  list(
    m = NA,
    fn = function(par) {
      x <- par[1]
      y <- par[2]
      yy <- y * y
      yyy <- yy * y

      f1 <- 1.5 - x * (1 - y)
      f2 <- 2.25 - x * (1 - yy)
      f3 <- 2.625 - x * (1 - yyy)

      f1 * f1 + f2 * f2 + f3 * f3
    },
    gr = function(par) {
      x <- par[1]
      y <- par[2]

      yy <- y * y
      yyy <- yy * y

      f1 <- 1.5 - x * (1 - y)
      f2 <- 2.25 - x * (1 - yy)
      f3 <- 2.625 - x * (1 - yyy)

      c(
        2 * (y - 1) * f1 + 2 * (yy - 1) * f2 + 2 * (yyy - 1) * f3,
        2 * x * f1 + 4 * x * y * f2 + 6 * x * yy * f3
      )
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      h <- matrix(0.0, nrow=2, ncol=2)
      y5 <- c(1.5, 2.25, 2.625) 
      for (i in 1:3) { # fixed m at 3
         t1 <- y5[i] - x1*( 1.0 - x2 ^ i )
          h[1,1] <- h[1,1] + 2.0*( x2 ^ i - 1.0 ) ^ 2
          h[1,2] <- h[1,2] + 2.0*( t1*i*x2 ^ (i-1) + i*x1*x2 ^ (i-1)*( x2 ^ i - 1.0 ) )
          h[2,2] <- h[2,2] + 2.0*( t1*( i - 1 )*i*x1*x2 ^ ( i - 2 ) + ( i*x1*x2 ^ ( i - 1 ) ) ^ 2 )
          h[2,1] <- h[1,2]
      }
      h
    },
    fg = function(par) {
      x <- par[1]
      y <- par[2]
      yy <- y * y
      yyy <- yy * y

      f1 <- 1.5 - x * (1 - y)
      f2 <- 2.25 - x * (1 - yy)
      f3 <- 2.625 - x * (1 - yyy)

      list(
        fn = f1 * f1 + f2 * f2 + f3 * f3,
        gr = c(
          2 * (y - 1) * f1 + 2 * (yy - 1) * f2 + 2 * (yyy - 1) * f3,
          2 * x * f1 + 4 * x * y * f2 + 6 * x * yy * f3
        )
      )
    },
    x0 = c(1, 1)
  )
}
