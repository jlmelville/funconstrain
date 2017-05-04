context("Extended Powell")
testfun <- ex_powell()
min_x <- rep(0, 12)
min_fx <- 0

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(12))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(12)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x), rep(0, 12))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0(12), fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS")
  expect_equal(res$par, min_x, tol = 1e-3)
  expect_equal(res$value, min_fx)
})
