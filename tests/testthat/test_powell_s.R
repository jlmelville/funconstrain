context("Powell Singular")
testfun <- powell_s()
min_x <- rep(0, 4)
min_fx <- 0
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0, tol = 1e-4)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 4))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS", control = list(maxit = 1000))
  expect_equal(res$par, min_x, tol = 1e-3)
  expect_equal(res$value, min_fx)
})
