#' Watson Function
#'
#' Test function 20 from the More', Garbow and Hillstrom paper.
#'
#' \itemize{
#'   \item Dimensions: \code{2 <= n <= 31}, \code{m = 31}.
#'   \item Minima: \code{f = 2.28767...e-3} if \code{n = 6};
#'   \code{f = 1.39976...e-6} if \code{n = 9};
#'   \code{f = 4.72238...e-10} if \code{n = 12}
#' }
#'
#' The number of variables in the function, \code{n}, is determined by the
#' length of the vector passed to the function and gradient routines. See
#' the 'Examples' section.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fn} Function which calculates the value given input
#'   parameter vector.
#'   \item \code{gr} Function which calculates the gradient vector given input
#'   parameter vector.
#'   \item \code{fg} Function which calculates both the value and gradient,
#'   given the parameter vector, returning a list.
#'   \item \code{x0} Function returning the standard starting point, given
#'   \code{n}, the number of variables desired.
#' }
#' @references
#' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#' \url{https://doi.org/10.1145/355934.355936}
#'
#' Kowalik, J. S., & Osborne, M. R. (1968).
#' \emph{Methods for unconstrained optimization problems.}
#' New York, NY: Elsevier North-Holland.
#' @examples
#' wat <- watson()
#' # 6 variable problem using the standard starting point
#' x0_6 <- wat$x0(6)
#' res_6 <- stats::optim(x0_6, wat$fn, wat$gr, method = "L-BFGS-B")
#' # Standing starting point with 9 variables
#' x0_9 <- wat$x0(9)
#' res_9 <- stats::optim(x0_9, wat$fn, wat$gr, method = "L-BFGS-B")
#' # Create your own 3 variable starting point
#' res_3 <- stats::optim(c(0.1, 0.2, 0.3), wat$fn, wat$gr, method = "L-BFGS-B")
#' @export
watson <- function() {
  # m = 31
  list(
    fn = function(par) {
      n <- length(par)
      if (!(2 <= n && n <= 31)) {
        stop("Watson: n must be between 2-31")
      }
      fsum <- 0
      for (i in 1:29) {
        ti <- i / 29

        fa <- 0
        fb <- 0
        tij2 <- 1 / ti
        tij1 <- 1
        for (j in 1:n) {
          fa <- fa + (j - 1) * par[j] * tij2
          fb <- fb + par[j] * tij1
          tij2 <- tij1
          tij1 <- tij1 * ti
        }
        fi <- fa - fb * fb - 1
        fsum <- fsum + fi * fi
      }
      f30 <- par[1]
      f31 <- par[2] - par[1] ^ 2 - 1

      fsum + f30 * f30 + f31 * f31
    },
    gr = function(par) {
      n <- length(par)
      if (!(2 <= n && n <= 31)) {
        stop("Watson: n must be between 2-31")
      }
      grad <- rep(0, n)
      for (i in 1:29) {
        ti <- i / 29

        fa <- 0
        fb <- 0

        tij2 <- 1 / ti
        tij1 <- 1
        for (j in 1:n) {
          fa <- fa + (j - 1) * par[j] * tij2
          fb <- fb + par[j] * tij1

          tij2 <- tij1
          tij1 <- tij1 * ti
        }
        fi <- fa - fb * fb - 1

        tij2 <- 1 / ti
        tij1 <- 1
        fi2 <- 2 * fi
        for (j in 1:n) {
          grad[j] <- grad[j] + fi2 * ((j - 1) * tij2 - 2 * fb * tij1)
          tij2 <- tij1
          tij1 <- tij1 * ti
        }
      }
      f31 <- par[2] - par[1] ^ 2 - 1

      grad[1] <- grad[1] + 2 * par[1]
      grad[1] <- grad[1] - 4 * par[1] * f31
      grad[2] <- grad[2] + 2 * f31

      grad

    },
    fg = function(par) {
      n <- length(par)
      if (!(2 <= n && n <= 31)) {
        stop("Watson: n must be between 2-31")
      }

      fsum <- 0
      grad <- rep(0, n)
      for (i in 1:29) {
        ti <- i / 29

        fa <- 0
        fb <- 0

        tij2 <- 1 / ti
        tij1 <- 1
        for (j in 1:n) {
          fa <- fa + (j - 1) * par[j] * tij2
          fb <- fb + par[j] * tij1

          tij2 <- tij1
          tij1 <- tij1 * ti
        }
        fi <- fa - fb * fb - 1
        fsum <- fsum + fi * fi

        tij2 <- 1 / ti
        tij1 <- 1
        fi2 <- 2 * fi
        for (j in 1:n) {
          grad[j] <- grad[j] + fi2 * ((j - 1) * tij2 - 2 * fb * tij1)
          tij2 <- tij1
          tij1 <- tij1 * ti
        }
      }
      f30 <- par[1]
      f31 <- par[2] - par[1] ^ 2 - 1

      fsum <- fsum + f30 * f30 + f31 * f31

      grad[1] <- grad[1] + 2 * par[1]
      grad[1] <- grad[1] - 4 * par[1] * f31
      grad[2] <- grad[2] + 2 * f31

      list(
        fn = fsum,
        gr = grad
      )
    },
    x0 = function(n) {
      if (!(2 <= n && n <= 31)) {
        stop("Watson: n must be between 2-31")
      }
      rep(0, n)
    }
  )
}
