gfd <- function(par, fn, rel_eps = sqrt(.Machine$double.eps)) {
  g <- rep(0, length(par))
  for (i in 1:length(par)) {
    oldx <- par[i]
    if (oldx != 0) {
      eps <- oldx * rel_eps
    }
    else {
      eps <- 1e-3
    }
    par[i] <- oldx + eps
    fplus <- fn(par)

    par[i] <- oldx - eps
    fminus <- fn(par)
    par[i] <- oldx

    g[i] <- (fplus - fminus) / (2 * eps)
  }
  g
}

make_gfd <- function(fn, rel_eps = sqrt(.Machine$double.eps)) {
  function(par) {
    gfd(par, fn, rel_eps)
  }
}

# test analytical gradient equals finite difference gradient at par
expect_gfd <- function(testfun, par, tol = 1e-6) {
  fd <- make_gfd(testfun$fn)(par)
  an <- testfun$gr(par)

  expect_equal(fd, an, tol = tol)
}
