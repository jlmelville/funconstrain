#' Brown and Dennis Function
#'
#' Test function 16 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{n = 4}, \code{m >= n}.
#'   \item Minima: \code{f = 85822.2} if \code{m = 20}.
#' }
#'
#' @param m Number of terms in the objective function. Should be equal to or
#' greater than \code{n}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Function which calculates the gradient vector given input
#'   parameter vector.
#'   \item \code{fg} A function which calculates both the value and gradient,
#'   given the parameter vector, returning a list.
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
#' @export
brown_den <- function(m = 20) {
  list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      fsum <- 0
      for (i in 1:m) {
        ti <- i * 0.2
        l <- x1 + ti * x2 - exp(ti)
        r <- x3 + x4 * sin(ti) - cos(ti)
        f <- l * l + r * r
        fsum <- fsum + f * f
      }
      fsum
    },
    gr = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      grad <- c(0, 0, 0, 0)
      for (i in 1:m) {
        ti <- i * 0.2
        sinti <- sin(ti)
        l <- x1 + ti * x2 - exp(ti)
        r <- x3 + x4 * sinti - cos(ti)
        f <- l * l + r * r
        lf4 <- 4 * l * f
        rf4 <- 4 * r * f

        grad[1] <- grad[1] + lf4
        grad[2] <- grad[2] + lf4 * ti
        grad[3] <- grad[3] + rf4
        grad[4] <- grad[4] + rf4 * sinti
      }
      grad
    },
    fg = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      x3 <- par[3]
      x4 <- par[4]

      fsum <- 0
      grad <- c(0, 0, 0, 0)
      for (i in 1:m) {
        ti <- i * 0.2
        sinti <- sin(ti)
        l <- x1 + ti * x2 - exp(ti)
        r <- x3 + x4 * sinti - cos(ti)
        f <- l * l + r * r
        lf4 <- 4 * l * f
        rf4 <- 4 * r * f

        fsum <- fsum + f * f

        grad[1] <- grad[1] + lf4
        grad[2] <- grad[2] + lf4 * ti
        grad[3] <- grad[3] + rf4
        grad[4] <- grad[4] + rf4 * sinti
      }

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = c(25, 5, -5, 1)
  )
}
