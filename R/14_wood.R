#' Wood Function
#'
#' Test function 14 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 4}, number of summand
#'   functions \code{m = 6}.
#'   \item Minima: \code{f = 0} at \code{rep(1, 4)}.
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
#' Colville, A. R. (1968)
#' \emph{A comparative study on nonlinear programming codes} (Report No.
#' 320-2949).
#' New York, NY: IBM.
#'
#' @examples
#' fun <- wood()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
#' "L-BFGS-B")
#' @export
wood <- function() {
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      sqrt10 <- sqrt(10)

      f1 <- 10 * (x2 - x1 ^ 2)
      f2 <- 1 - x1
      f3 <- sqrt(90) * (x4 - x3 ^ 2)
      f4 <- 1 - x3
      f5 <- sqrt10 * (x2 + x4 - 2)
      f6 <- (x2 - x4) / sqrt10

      f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4 + f5 * f5 + f6 * f6
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      x1_2 <- x1 * x1
      x1_3 <- x1_2 * x1
      x3_2 <- x3 * x3
      x3_3 <- x3_2 * x3

      c(
        400 * x1_3 + (2 - 400 * x2) * x1 - 2,
        (1101 * x2 - 1000 * x1_2 + 99 * x4 - 200) * 0.2,
        360 * x3_3 + (2 - 360 * x4) * x3 - 2,
        (1001 * x4 - 900 * x3_2 + 99 * x2 - 200) * 0.2
      )
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      h <- matrix(0.0, ncol=4, nrow=4)
      h[1,1] <- 1.2e+3*x1 ^ 2 - 4.0e+2*x2 + 2.0
      h[2,2] <- 2.202e+2
      h[3,3] <- 1.08e+3*x3 ^ 2 - 3.6e+2*x4 + 2.0
      h[4,4] <- 2.002e+2
      h[1,2] <- -4.0e+2*x1
      h[1,3] <- 0.0
      h[2,3] <- 0.0
      h[1,4] <- 0.0
      h[2,4] <- 1.98e+1
      h[3,4] <- -3.6e+2*x3
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

      sqrt10 <- sqrt(10)

      x1_2 <- x1 * x1
      x1_3 <- x1_2 * x1
      x3_2 <- x3 * x3
      x3_3 <- x3_2 * x3

      f1 <- 10 * (x2 - x1_2)
      f2 <- 1 - x1
      f3 <- sqrt(90) * (x4 - x3_2)
      f4 <- 1 - x3
      f5 <- sqrt10 * (x2 + x4 - 2)
      f6 <- (x2 - x4) / sqrt10

      list(
        fn = f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4 + f5 * f5 + f6 * f6,
        gr = c(
          400 * x1_3 + (2 - 400 * x2) * x1 - 2,
          (1101 * x2 - 1000 * x1_2 + 99 * x4 - 200) * 0.2,
          360 * x3_3 + (2 - 360 * x4) * x3 - 2,
          (1001 * x4 - 900 * x3_2 + 99 * x2 - 200) * 0.2
        )
      )
    },
    x0 = c(-3, -1, -3, -1)
  )
}
