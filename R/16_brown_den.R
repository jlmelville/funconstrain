#' Brown and Dennis Function
#'
#' Test function 16 from the More', Garbow and Hillstrom paper.
#'
#' The objective function is the sum of \code{m} functions, each of \code{n}
#' parameters.
#'
#' \itemize{
#'   \item Dimensions: Number of parameters \code{n = 4}, number of summand
#'   functions \code{m >= n}.
#'   \item Minima: \code{f = 85822.2} if \code{m = 20}.
#' }
#'
#' @param m Number of summand functions in the objective function. Should be
#'   equal to or greater than 4.
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
#' Brown, K. M., & Dennis, J. E. (1971).
#' \emph{New computational algorithms for minimizing a sum of squares of
#' nonlinear functions} (Report No. 71-6).
#' New Haven, CT: Department of Computer Science, Yale University.
#'
#' @examples
#' # Use 10 summand functions
#' fun <- brown_den(m = 10)
#' # Optimize using the standard starting point
#' x0 <- fun$x0
#' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#' "L-BFGS-B")
#' # Use your own starting point
#' res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
#' "L-BFGS-B")
#'
#' # Use 20 summand functions
#' fun20 <- brown_den(m = 20)
#' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
#' @export
brown_den <- function(m = 20) {
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      ti <- (1:m) * 0.2
      l <- x1 + ti * x2 - exp(ti)
      r <- x3 + x4 * sin(ti) - cos(ti)
      f <- l * l + r * r
      sum(f * f)
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      ti <- (1:m) * 0.2
      sinti <- sin(ti)
      l <- x1 + ti * x2 - exp(ti)
      r <- x3 + x4 * sinti - cos(ti)
      f <- l * l + r * r
      lf4 <- 4 * l * f
      rf4 <- 4 * r * f
      c(
        sum(lf4),
        sum(lf4 * ti),
        sum(rf4),
        sum(rf4 * sinti)
      )
    },
    he = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]
      h <- matrix(0.0, ncol=4, nrow=4)
      for (i in 1:m) {
        d1 <- i/5.0
        d2 <- sin( d1 )
        t1 <- x1 + d1*x2 - exp( d1 )
        t2 <- x3 + d2*x4 - cos( d1 )
        t <- 8.0 * t1 * t2
        s1 <- 12.0*t1 ^ 2 + 4.0*t2 ^ 2
        s2 <- 12.0*t2 ^ 2 + 4.0*t1 ^ 2
        h[1,1] <- h[1,1] + s1
        h[1,2] <- h[1,2] + s1*d1
        h[2,2] <- h[2,2] + s1*d1 ^ 2
        h[1,3] <- h[1,3] + t
        h[2,3] <- h[2,3] + t*d1
        h[3,3] <- h[3,3] + s2
        h[1,4] <- h[1,4] + t*d2
        h[2,4] <- h[2,4] + t*d1*d2
        h[3,4] <- h[3,4] + s2*d2
        h[4,4] <- h[4,4] + s2*d2 ^ 2
      }
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

      ti <- (1:m) * 0.2
      sinti <- sin(ti)
      l <- x1 + ti * x2 - exp(ti)
      r <- x3 + x4 * sinti - cos(ti)
      f <- l * l + r * r
      lf4 <- 4 * l * f
      rf4 <- 4 * r * f

      fsum <- sum(f * f)
      grad <- c(
        sum(lf4),
        sum(lf4 * ti),
        sum(rf4),
        sum(rf4 * sinti)
      )

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(25, 5, -5, 1)
  )
}
