#' Gulf Research and Development Function
#'
#' Test function 11 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 3}, number of summand
#'   functions \code{3 <= m <= 100}.
#'   \item Minima: \code{f = 0} at \code{(50, 25, 1.5)}.
#' }
#'
#' Note that the equation as published by More', Garbow and Hillstrom (1981)
#' contains an error, where the symbol 'mi' should be interpreted as a minus
#' sign. The corrected version can be found in Jamil and Xang (2013), and no
#' doubt several other publications. The Jamil and Xang equation unfortunately
#' contains its own minor errors, but you can piece together the correct
#' equation from these two sources without too much trouble.
#'
#' @param m Number of summand functions in the objective function. Should be
#'   between 3 and 100, according to the MGH paper. Default value is 99, which
#'   Jamil and Xang (2013) list as the only valid value.
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
#' Jamil, M., & Yang, X. S. (2013).
#' A literature survey of benchmark functions for global optimisation problems.
#' \emph{International Journal of Mathematical Modelling and Numerical
#' Optimisation}, \emph{4}(2), 150-194.
#' \doi{doi.org/10.1504/IJMMNO.2013.055204}
#' \doi{arxiv.org/abs/1308.4008}
#'
#' @examples
#' # Use 10 summand functions
#' fun <- gulf(m = 10)
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
#'
#' # Use 20 summand functions
#' fun20 <- gulf(m = 20)
#' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
#' @export
gulf <- function(m = 99) {
  # can be between 3 and 100
  if (m < 3 || m > 100) {
    stop("Gulf research and development function: m must be between 3 and 100")
  }

  p66 <- 2 / 3

  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      ti <- 1:m * 0.01
      y <- 25 + (-50 * log(ti)) ^ p66
      fi <- exp(-(abs(x2 - y) ^ x3) / x1) - ti
      sum(fi * fi)
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      ti <- 1:m * 0.01
      y <- 25 + (-50 * log(ti)) ^ p66

      x2y <- x2 - y
      ax2y <- abs(x2y)
      x2yz <- ax2y ^ x3
      e <- exp(-x2yz / x1)
      fi <- e - ti
      efxyz2 <- 2 * e * fi * x2yz

      dx <- sum(efxyz2 / (x1 * x1))
      dy <- -sum(efxyz2 * x3 / (x1 * x2y))
      dz <- -sum(efxyz2 * log(ax2y) / x1)

      c(dx, dy, dz)
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      h <- matrix(0.0, ncol=3, nrow=3)
      d1 <- p66
      for (i in 1:m) {
          arg <- 0.01*i
          r <- ( -50.0*log( arg ) ) ^ d1 + 25.0 - x2
          t1 <- abs( r ) ^ x3/x1
          t2 <- exp( -t1 )
          t3 <- t1*t2*( t1*t2 + ( t1 - 1.0 )*( t2 - arg ) )
          t <- t1*t2*( t2 - arg )
          logr <- log( abs( r ) )
          h[1,1] <- h[1,1] + t3 - t
          h[1,2] <- h[1,2] + t3/r
          h[2,2] <- h[2,2] + ( t + x3*t3 )/r ^ 2
          h[1,3] <- h[1,3] - t3*logr
          h[2,3] <- h[2,3] + (t-x3*t3*logr)/r
          h[3,3] <- h[3,3] + t3*logr ^ 2
       }

       h[1,1] <- h[1,1] / x1 ^ 2
       h[1,2] <- h[1,2]*x3/x1
       h[2,2] <- h[2,2]*x3
       h[1,3] <- h[1,3] / x1
       h <- 2.0 * h
       h[2,1] <- h[1,2]
       h[3,1] <- h[1,3]
       h[3,2] <- h[2,3]
       h
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]

      ti <- 1:m * 0.01
      y <- 25 + (-50 * log(ti)) ^ p66

      x2y <- x2 - y
      ax2y <- abs(x2y)
      x2yz <- ax2y ^ x3
      e <- exp(-x2yz / x1)
      fi <- e - ti
      efxyz2 <- 2 * e * fi * x2yz

      dx <- sum(efxyz2 / (x1 * x1))
      dy <- -sum(efxyz2 * x3 / (x1 * x2y))
      dz <- -sum(efxyz2 * log(ax2y) / x1)

      fsum <- sum(fi * fi)
      grad <- c(dx, dy, dz)

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(5, 2.5, 0.15),
    fmin = 0,
    xmin = c(50, 25, 1.5)
  )
}
