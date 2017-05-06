#' Extended Powell Function
#'
#' Test function 22 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n} variable but a multiple of 4, \code{m = n}.
#'   \item Minima: \code{f = 0} at \code{rep(0, n)}
#' }
#'
#' The number of variables in the function, \code{n}, is determined by the
#' length of the vector passed to the function and gradient routines. See
#' the 'Examples' section.
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
#'   \item \code{x0} Function returning the standard starting point, given
#'   \code{n}, the number of variables desired.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Spedicato, E. (1975).
#' \emph{Computational experience with quasi-Newton algorithms for minimization
#' problems of moderately large size} (Report CISE-N-175).
#' Segrate, Milano: Computing Center, CISE.
#' @examples
#' expow <- ex_powell()
#' # 12 variable problem using the standard starting point
#' x0_12 <- expow$x0(12)
#' res_12 <- stats::optim(x0_12, expow$fn, expow$gr, method = "L-BFGS-B")
#' # Standing starting point with 8 variables
#' res_8 <- stats::optim(expow$x0(8), expow$fn, expow$gr, method = "L-BFGS-B")
#' # Create your own 4 variable starting point
#' res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), expow$fn, expow$gr,
#'                       method = "L-BFGS-B")
#' @export
ex_powell <- function() {

  sqrt10 <- sqrt(10)
  sqrt5 <- sqrt(5)

  list(
    fn = function(par) {
      n <- length(par)
      if (n %% 4 != 0) {
        stop("Extended Powell: n must be a multiple of 4")
      }
      fsum <- 0
      for (i in 1:(n / 4)) {
        w <- 4 * i - 3
        x <- 4 * i - 2
        y <- 4 * i - 1
        z <- 4 * i

        fw <- par[w] + 10 * par[x]
        fx <- sqrt5 * (par[y] - par[z])
        fy <- (par[x] - 2 * par[y]) ^ 2
        fz <- sqrt10 * (par[w] - par[z]) ^ 2

        fsum <- fsum + fw * fw + fx * fx + fy * fy + fz * fz
      }

      fsum
    },
    gr = function(par) {
      n <- length(par)
      if (n %% 4 != 0) {
        stop("Extended Powell: n must be a multiple of 4")
      }
      grad <- rep(0, n)

      for (i in 1:(n / 4)) {
        w <- 4 * i - 3
        x <- 4 * i - 2
        y <- 4 * i - 1
        z <- 4 * i

        fw <- par[w] + 10 * par[x]
        pyz <- par[y] - par[z]
        px2py <- par[x] - 2 * par[y]
        px2py2 <- px2py * px2py
        px2py3 <- px2py2 * px2py
        pwz <- par[w] - par[z]
        pwz2 <- pwz * pwz
        pwz3 <- pwz2 * pwz

        grad[w] <- grad[w] + 40 * pwz3 + 2 * fw
        grad[x] <- grad[x] + 4 * px2py3 + 20 * fw
        grad[y] <- grad[y] + 10 * pyz - 8 * px2py3
        grad[z] <- grad[z] - 10 * pyz - 40 * pwz3
      }
      grad
    },
    fg = function(par) {
      n <- length(par)
      if (n %% 4 != 0) {
        stop("Extended Powell: n must be a multiple of 4")
      }

      fsum <- 0
      grad <- rep(0, n)
      for (i in 1:(n / 4)) {
        w <- 4 * i - 3
        x <- 4 * i - 2
        y <- 4 * i - 1
        z <- 4 * i

        fw <- par[w] + 10 * par[x]
        pyz <- par[y] - par[z]
        fx <- sqrt5 * pyz
        px2py <- par[x] - 2 * par[y]
        fy <- px2py * px2py
        px2py3 <- fy * px2py
        pwz <- par[w] - par[z]
        pwz2 <- pwz * pwz
        pwz3 <- pwz2 * pwz
        fz <- sqrt10 * pwz2

        fsum <- fsum + fw * fw + fx * fx + fy * fy + fz * fz

        grad[w] <- grad[w] + 40 * pwz3 + 2 * fw
        grad[x] <- grad[x] + 4 * px2py3 + 20 * fw
        grad[y] <- grad[y] + 10 * pyz - 8 * px2py3
        grad[z] <- grad[z] - 10 * pyz - 40 * pwz3
      }
      grad

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n = 12) {
      if (n %% 4 != 0) {
        stop("Extended Powell: n must be a multiple of 4")
      }
      rep(c(3, -1, 0, 1), n / 4)
    }
  )
}
