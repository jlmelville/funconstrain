#' Osborne 2 Function
#'
#' Test function 19 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 11}, number of summand
#'   functions \code{m = 65}.
#'   \item Minima: \code{f = 4.01377...e-2}.
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
#' Osborne, M. R. (1972).
#' Some aspects of nonlinear least squares calculations.
#' In F. A. Lootsma (Ed.),
#' \emph{Numerical methods for nonlinear optimization} (pp. 171-189).
#' New York, NY: Academic Press.
#'
#' @examples
#' fun <- osborne_2()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(1:11, fun$fn, fun$gr, method = "L-BFGS-B")
#' @export
osborne_2 <- function() {
  m <- 65
  n <- 11
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

      ti <- ((1:m) - 1) * 0.1
      fi <- y - (
          x1 * exp(-ti * x5) +
          x2 * exp(-(ti - x9) ^ 2 * x6) +
          x3 * exp(-(ti - x10) ^ 2 * x7) +
          x4 * exp(-(ti - x11) ^ 2 * x8))
      sum(fi * fi)

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

      ti <- ((1:m) - 1) * 0.1

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

      f <- y - (f15 + f269 + f3710 + f4811)
      f2 <- 2 * f

      grad <- c(
        sum(-f2 * f5),
        sum(-f2 * f69),
        sum(-f2 * f710),
        sum(-f2 * f811),
        sum(f2 * f15 * ti),
        sum(f2 * f9s * f269),
        sum(f2 * f10s * f3710),
        sum(f2 * f11s * f4811),
        sum(-f2 * 2 * x6 * x2 * f9 * f69),
        sum(-f2 * 2 * x7 * x3 * f10 * f710),
        sum(-f2 * 2 * x8 * x4 * f11 * f811)
      )

      grad
    },
    he = function(par) {
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
      h <- matrix(0.0, ncol=11, nrow=11)

      for (i in 1:m) {
          t1 <- (i - 1)/ 10.0
          d1 <- exp( - t1*x5 )
          d2 <- exp( - ( t1 - x9 ) ^ 2*x6 )
          d3 <- exp( - ( t1 - x10 ) ^ 2*x7 )
          d4 <- exp( - ( t1 - x11 ) ^ 2*x8 )
          s1 <- y[i] - ( x1*d1 + x2*d2 + x3*d3 + x4*d4 )

          w1 <- - c( d1, d2, d3, d4, - t1*x1*d1, 
                 - ( t1 - x9 ) ^ 2*x2*d2,
                 - ( t1 - x10 ) ^ 2*x3*d3,
                 - ( t1 - x11 ) ^ 2*x4*d4,
                 2.0*( t1 - x9 )*x6*x2*d2,
                 2.0*( t1 - x10 )*x7*x3*d3,
                 2.0*( t1 - x11 )*x8*x4*d4 )

          for (k in 1:n) {
             for (j in 1:k) {
                h[j,k] <- h[j,k] + 2.0*w1[j]*w1[k]
             }
          }

          h[ 1, 5] <- h[ 1, 5] + 2.0*s1*t1*d1
          h[ 5, 5] <- h[ 5, 5] - 2.0*s1*t1 ^ 2*x1*d1
          h[ 2, 6] <- h[ 2, 6] + 2.0*s1*( t1 - x9 ) ^ 2*d2
          h[ 6, 6] <- h[ 6, 6] - 2.0*s1*( t1 - x9 ) ^ 4*x2*d2
          h[ 3, 7] <- h[ 3, 7] + 2.0*s1*( t1 - x10 ) ^ 2*d3
          h[ 7, 7] <- h[ 7, 7] - 2.0*s1*( t1 - x10 ) ^ 4*x3*d3
          h[ 4, 8] <- h[ 4, 8] + 2.0*s1*( t1 - x11 ) ^ 2*d4
          h[ 8, 8] <- h[ 8, 8] - 2.0*s1*( t1 - x11 ) ^ 4*x4*d4
          h[ 2, 9] <- h[ 2, 9] - 4.0*s1*( t1 - x9 )*x6*d2
          h[ 6, 9] <- h[ 6, 9] + 4.0*s1*( t1 - x9 )*x2*d2*( ( t1 - x9 ) ^ 2*x6 - 1.0 )
          h[ 9, 9] <- h[ 9, 9] - 4.0*s1*x6*x2*d2*( 2.0*( t1 - x9 ) ^ 2*x6 - 1.0 )
          h[ 3,10] <- h[ 3,10] - 4.0*s1*( t1 - x10 )*x7*d3
          h[ 7,10] <- h[ 7,10] + 4.0*s1*x3*( t1 - x10 )*d3*( ( t1 - x10 ) ^ 2*x7 - 1.0 )
          h[10,10] <- h[10,10] - 4.0*s1*x7*x3*d3*( 2.0*( t1 - x10 ) ^ 2*x7 - 1.0 )
          h[ 4,11] <- h[ 4,11] - 4.0*s1*( t1 - x11 )*x8*d4
          h[ 8,11] <- h[ 8,11] + 4.0*s1*x4*( t1 - x11 )*d4*( ( t1 - x11 ) ^ 2*x8 - 1.0 )
          h[11,11] <- h[11,11] - 4.0*s1*x8*x4*d4*( 2.0*( t1 - x11 ) ^ 2*x8 - 1.0 )
       }
      for (j in 1:10) {
        for (k in (j+1):11) {
          h[k,j] <- h[j,k]        
        }
      }
      h
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

      ti <- ((1:m) - 1) * 0.1

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

      f <- y - (f15 + f269 + f3710 + f4811)
      fsum <- sum(f * f)

      f2 <- 2 * f
      grad <- c(
        sum(-f2 * f5),
        sum(-f2 * f69),
        sum(-f2 * f710),
        sum(-f2 * f811),
        sum(f2 * f15 * ti),
        sum(f2 * f9s * f269),
        sum(f2 * f10s * f3710),
        sum(f2 * f11s * f4811),
        sum(-f2 * 2 * x6 * x2 * f9 * f69),
        sum(-f2 * 2 * x7 * x3 * f10 * f710),
        sum(-f2 * 2 * x8 * x4 * f11 * f811)
      )

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(1.3, 0.65, 0.65, 0.7, 0.6, 3, 5, 7, 2, 4.5, 5.5)
    fmin = f = 4.013774e-2,
    xmin = c(1.309977, 0.4315538, 0.6336617, 0.5994305, 0.7541832, 0.9042886, 1.3658118, 4.823699, 2.398685, 4.568875, 5.675341)
  )
}
