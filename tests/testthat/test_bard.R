context("Bard")

testfun <- bard()
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0, tol = 1e-3)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(0.08241, 1.13304, 2.34370))
  expect_equal(gr0, c(0, 0, 0), tol = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(0.08241, 1.13304, 2.34370)), 8.21487e-3)
  expect_equal(testfun$fn(c(0.8406, 1e6, 1e6)), 17.4286, tol = 1e-5)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B")
  expect_equal(res$par, c(0.08241, 1.13304, 2.34370), tol = 1e-5)
  expect_equal(res$value, 8.21487e-3)
})
