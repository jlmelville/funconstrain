#' Biggs EXP6 Function
#'
#' Test function 18 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 6}, number of summand
#'   functions \code{m >= n}.
#'   \item Minima: \code{f = 5.65565...e-3} if \code{m = 13};
#'   not reported in the MGH (1981) paper is \code{(f = 0)} at
#'   \code{c(1, 10, 1, 5, 4, 3)} and \code{c(4, 10, 3, 5, 1, 1)}
#'   (and probably others) for all \code{m} (probably: I stopped testing after
#'   \code{m = 1000}).
#' }
#'
#' @param m Number of summand functions in the objective function. Should be
#'   equal to or greater than 6.
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
#' Biggs, M. C. (1971).
#' Minimization algorithms making use of non-quadratic properties of the
#' objective function.
#' \emph{IMA Journal of Applied Mathematics}, \emph{8}(3), 315-327.
#'
#' @examples
#' fun <- biggs_exp6()
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(1:6, fun$fn, fun$gr, method = "L-BFGS-B")
#'
#' # Use 20 summand functions
#' fun20 <- biggs_exp6(m = 20)
#' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
#' @export
biggs_exp6 <- function(m = 13) {
  if (m < 6) {
    stop("Biggs EXP6 function must have m >= 6")
  }

  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      ti <- 0.1 * (1:m)
      yi <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
      fi <- x3 * exp(-ti * x1) - x4 * exp(-ti * x2) + x6 * exp(-ti * x5) - yi
      sum(fi * fi)

    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      ti <- 0.1 * (1:m)
      yi <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
      f1 <- exp(-ti * x1)
      f31 <- x3 * f1
      f2 <- exp(-ti * x2)
      f42 <- x4 * f2
      f5 <- exp(-ti * x5)
      f65 <- x6 * exp(-ti * x5)
      fi <- f31 - f42 + f65 - yi
      fi2 <- 2 * fi
      tfi2 <- fi2 * ti

      c(
        -sum(tfi2 * f31),
         sum(tfi2 * f42),
         sum(fi2 * f1),
        -sum(fi2 * f2),
        -sum(tfi2 * f65),
         sum(fi2 * f5)
      )
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]
      h <- matrix(0.0, ncol=6, nrow=6)

      for (i in 1:m) {
          d1 <- i/10.0
          d2 <- exp( - d1 ) - 5.0*exp( - 10.0*d1 ) + 3.0*exp( - 4.0*d1 )
          s1 <- exp( - d1*x1 )
          s2 <- exp( - d1*x2 )
          s3 <- exp( - d1*x5 )
          t <- x3*s1 - x4*s2 + x6*s3 - d2
          d2 <- d1 ^ 2
          s1s2 <- s1 * s2
          s1s3 <- s1 * s3
          s2s3 <- s2 * s3
          h[1,1] <- h[1,1] + d2*s1*( t + x3*s1 )
          h[1,2] <- h[1,2] - d2*s1s2
          h[2,2] <- h[2,2] - d2*s2*( t - x4*s2 )
          h[1,3] <- h[1,3] - d1*s1*( t + x3*s1 )
          h[2,3] <- h[2,3] + d1*s1s2
          h[3,3] <- h[3,3] + s1 ^ 2
          h[1,4] <- h[1,4] + d1*s1s2
          h[2,4] <- h[2,4] + d1*s2*( t - x4*s2 )
          h[3,4] <- h[3,4] - s1s2
          h[4,4] <- h[4,4] + s2 ^ 2
          h[1,5] <- h[1,5] + d2*s1s3
          h[2,5] <- h[2,5] - d2*s2s3
          h[3,5] <- h[3,5] - d1*s1s3
          h[4,5] <- h[4,5] + d1*s2s3
          h[5,5] <- h[5,5] + d2*s3*( t + x6*s3 )
          h[1,6] <- h[1,6] - d1*s1s3
          h[2,6] <- h[2,6] + d1*s2s3
          h[3,6] <- h[3,6] + s1s3
          h[4,6] <- h[4,6] - s2s3
          h[5,6] <- h[5,6] - d1*s3*( t + x6*s3 )
          h[6,6] <- h[6,6] + s3 ^ 2
       }
       h[1,1] <- x3*h[1,1]
       h[2,2] <- x4*h[2,2]
       h[5,5] <- x6*h[5,5]
       h[1,2] <- x3*x4*h[1,2]
       h[2,3] <- x4*h[2,3]
       h[1,4] <- x3*h[1,4]
       h[1,5] <- x3*x6*h[1,5]
       h[2,5] <- x4*x6*h[2,5]
       h[3,5] <- x6*h[3,5]
       h[4,5] <- x6*h[4,5]
       h[1,6] <- x3*h[1,6]
       h[2,6] <- x4*h[2,6]
       h <- 2.0 * h
      h[2,1] <- h[1,2]
      h[3,1] <- h[1,3]
      h[3,2] <- h[2,3]
      h[4,1] <- h[1,4]
      h[4,2] <- h[2,4]
      h[4,3] <- h[3,4]
      h[5,1] <- h[1,5]
      h[5,2] <- h[2,5]
      h[5,3] <- h[3,5]
      h[5,4] <- h[4,5]
      h[6,1] <- h[1,6]
      h[6,2] <- h[2,6]
      h[6,3] <- h[3,6]
      h[6,4] <- h[4,6]
      h[6,5] <- h[5,6]
      h
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      x5 <- par[5]
      x6 <- par[6]

      ti <- 0.1 * (1:m)
      yi <- exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti)
      f1 <- exp(-ti * x1)
      f31 <- x3 * f1
      f2 <- exp(-ti * x2)
      f42 <- x4 * f2
      f5 <- exp(-ti * x5)
      f65 <- x6 * exp(-ti * x5)
      fi <- f31 - f42 + f65 - yi
      fi2 <- 2 * fi
      tfi2 <- fi2 * ti

      fsum <- sum(fi * fi)
      grad <- c(
        -sum(tfi2 * f31),
        sum(tfi2 * f42),
        sum(fi2 * f1),
        -sum(fi2 * f2),
        -sum(tfi2 * f65),
        sum(fi2 * f5)
      )

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(1, 2, 1, 1, 1, 1)
  )
}
