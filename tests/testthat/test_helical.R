context("Helical Valley")

testfun <- helical()
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0, tol = 1e-3)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(1, 0, 0))
  expect_equal(gr0, c(0, 0, 0))
})
test_that("Function is zero at stated minima", {
  expect_equal(testfun$fn(c(1, 0, 0)), 0)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS", control = list(maxit = 1000))
  expect_equal(res$par, c(1, 0, 0))
  expect_equal(res$value, 0)
})
