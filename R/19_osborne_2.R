#' Osborne 2 Function
#'
#' Test function 19 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 11}, \code{m = 65}.
#'   \item Minima: \code{f = 4.01377...e-2}.
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
#' Osborne, M. R. (1972).
#' Some aspects of nonlinear least squares calculations.
#' In F. A. Lootsma (Ed.),
#' \emph{Numerical methods for nonlinear optimization} (pp. 171-189).
#' New York, NY: Academic Press.
#' @export
osborne_2 <- function() {
  m <- 65

  y <- c(1.366, 1.191, 1.112, 1.013, 0.991, 0.885, 0.831, 0.847, 0.786, 0.725,
         0.746, 0.679, 0.608, 0.655, 0.616, 0.606, 0.602, 0.626, 0.651, 0.724,
         0.649, 0.649, 0.694, 0.644, 0.624, 0.661, 0.612, 0.558, 0.533, 0.495,
         0.500, 0.423, 0.395, 0.375, 0.372, 0.391, 0.396, 0.405, 0.428, 0.429,
         0.523, 0.562, 0.607, 0.653, 0.672, 0.708, 0.633, 0.668, 0.645, 0.632,
         0.591, 0.559, 0.597, 0.625, 0.739, 0.710, 0.729, 0.720, 0.636, 0.581,
         0.428, 0.292, 0.162, 0.098, 0.054)

  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]
      x7 <- par[7]
      x8 <- par[8]
      x9 <- par[9]
      x10 <- par[10]
      x11 <- par[11]

      fsum <- 0
      for (i in 1:m) {
        ti <- (i - 1) * 0.1
        f <- y[i] - (
          x1 * exp(-ti * x5) +
            x2 * exp(-(ti - x9) ^ 2 * x6) +
            x3 * exp(-(ti - x10) ^ 2 * x7) +
            x4 * exp(-(ti - x11) ^ 2 * x8))
        fsum <- fsum + f * f
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]
      x7 <- par[7]
      x8 <- par[8]
      x9 <- par[9]
      x10 <- par[10]
      x11 <- par[11]

      grad <- rep(0, 11)
      for (i in 1:m) {
        ti <- (i - 1) * 0.1
        f5 <- exp(-ti * x5)
        f15 <- x1 * f5
        f9 <- ti - x9
        f9s <- f9 * f9
        f69 <- exp(-f9s * x6)
        f269 <- x2 * f69
        f10 <- ti - x10
        f10s <- f10 * f10
        f710 <- exp(-f10s * x7)
        f3710 <- x3 * f710
        f11 <- ti - x11
        f11s <- f11 * f11
        f811 <- exp(-f11s * x8)
        f4811 <- x4 * f811
        f <- y[i] - (f15 + f269 + f3710 + f4811)

        grad[1] <- grad[1] - 2 * f5 * f
        grad[2] <- grad[2] - 2 * f69 * f
        grad[3] <- grad[3] - 2 * f710 * f
        grad[4] <- grad[4] - 2 * f811 * f
        grad[5] <- grad[5] + 2 * ti * f15 * f
        grad[6] <- grad[6] + 2 * f9s * f269 * f
        grad[7] <- grad[7] + 2 * f10s * f3710 * f
        grad[8] <- grad[8] + 2 * f11s * f4811 * f
        grad[9] <- grad[9] - 4 * x6 * x2 * f9 * f69 * f
        grad[10] <- grad[10] - 4 * x7 * x3 * f10 * f710 * f
        grad[11] <- grad[11] - 4 * x8 * x4 * f11 * f811 * f

      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]
      x7 <- par[7]
      x8 <- par[8]
      x9 <- par[9]
      x10 <- par[10]
      x11 <- par[11]

      fsum <- 0
      grad <- rep(0, 11)
      for (i in 1:m) {
        ti <- (i - 1) * 0.1
        f5 <- exp(-ti * x5)
        f15 <- x1 * f5
        f9 <- ti - x9
        f9s <- f9 * f9
        f69 <- exp(-f9s * x6)
        f269 <- x2 * f69
        f10 <- ti - x10
        f10s <- f10 * f10
        f710 <- exp(-f10s * x7)
        f3710 <- x3 * f710
        f11 <- ti - x11
        f11s <- f11 * f11
        f811 <- exp(-f11s * x8)
        f4811 <- x4 * f811
        f <- y[i] - (f15 + f269 + f3710 + f4811)

        fsum <- fsum + f * f

        grad[1] <- grad[1] - 2 * f5 * f
        grad[2] <- grad[2] - 2 * f69 * f
        grad[3] <- grad[3] - 2 * f710 * f
        grad[4] <- grad[4] - 2 * f811 * f
        grad[5] <- grad[5] + 2 * ti * f15 * f
        grad[6] <- grad[6] + 2 * f9s * f269 * f
        grad[7] <- grad[7] + 2 * f10s * f3710 * f
        grad[8] <- grad[8] + 2 * f11s * f4811 * f
        grad[9] <- grad[9] - 4 * x6 * x2 * f9 * f69 * f
        grad[10] <- grad[10] - 4 * x7 * x3 * f10 * f710 * f
        grad[11] <- grad[11] - 4 * x8 * x4 * f11 * f811 * f
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(1.3, 0.65, 0.65, 0.7, 0.6, 3, 5, 7, 2, 4.5, 5.5)
  )
}
