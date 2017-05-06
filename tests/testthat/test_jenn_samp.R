context("Jennrich and Sampson")

testfun <- jenn_samp(m = 10)
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0, tol = 1e-3)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(rep(0.2578252, 2))
  expect_equal(gr0, rep(0, 2), tol = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(rep(0.2578252, 2)), 124.362, tol = 1e-5)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", lower = -1, upper = 1)
  expect_equal(res$par, rep(0.2578, 2), tol = 1e-4)
  expect_equal(res$value, 124.362, tol = 1e-4)
})
