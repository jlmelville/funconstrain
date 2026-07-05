#' Rosenbrock Function
#'
#' Test function 1 from the Moré, Garbow and Hillstrom paper.
#'
#' The objective function is the sum of `m` functions, each of `n`
#' parameters.
#'
#' - Dimensions: Number of parameters `n = 2`, number of summand
#'   functions `m = 2`.
#' - Minima: `f = 0` at `(1, 1)`
#'
#' @template factory-return
#' @references
#' Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' *ACM Transactions on Mathematical Software (TOMS)*, *7*(1), 17-41.
#' <https://doi.org/10.1145/355934.355936>
#'
#' Rosenbrock, H. (1960).
#' An automatic method for finding the greatest or least value of a function.
#' *The Computer Journal*, *3*(3), 175-184.
#' <https://doi.org/10.1093/comjnl/3.3.175>
#'
#' @examples
#' fun <- rosen()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
rosen <- function() {
  list(
    m = NA,
    fn = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      100 * (x2 - x1 * x1)^2 + (1 - x1)^2
    },
    gr = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      c(
        -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
        200 * (x2 - x1 * x1)
      )
    },
    he = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      h <- matrix(NA, nrow = 2, ncol = 2)
      t1 <- 10.0 * (x2 - x1^2)
      h[1, 1] <- 2.0 * (400.0 * x1^2 - 20.0 * t1 + 1.0)
      h[1, 2] <- -400.0 * x1
      h[2, 2] <- 200.0
      h[2, 1] <- h[1, 2]
      h
    },
    fg = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x2x11 <- x2 - x1 * x1
      list(
        fn = 100 * x2x11^2 + (1 - x1)^2,
        gr = c(
          -400 * x1 * x2x11 - 2 * (1 - x1),
          200 * x2x11
        )
      )
    },
    x0 = c(-1.2, 1),
    fmin = 0,
    xmin = c(1, 1)
  )
}
