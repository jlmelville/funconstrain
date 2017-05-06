#' Wood Function
#'
#' Test function 14 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 4}, \code{m = 6}.
#'   \item Minima: \code{f = 0} at \code{rep(1, 4)}.
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
#' Colville, A. R. (1968)
#' \emph{A comparative study on nonlinear programming codes} (Report No.
#' 320-2949).
#' New York, NY: IBM.
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
