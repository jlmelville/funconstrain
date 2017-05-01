context("Freudenstein and Roth")

testfun <- freud_roth()
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(c(5, 4)), c(0, 0))
  expect_equal(testfun$gr(c(11.413, -0.8968)), c(0, 0), tol = 1e-2)
})
test_that("Function is zero at stated minima", {
  expect_equal(testfun$fn(c(5, 4)), 0, tol = 1e-3)
  expect_equal(testfun$fn(c(11.413, -0.8968)), 48.984, tol = 1e-5)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS")
  expect_equal(res$par, c(11.413, -0.8968), tol = 1e-4)
  expect_equal(res$value, 48.984, tol = 1e-5)
})
