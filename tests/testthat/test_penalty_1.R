context("Penalty I")
testfun <- penalty_1()
min_x4 <- rep(0.2500075, 4)
min_fx4 <- 2.24997e-05 # actually only reach 2.249978e-05

min_x10 <- rep(0.158122, 10)
min_fx10 <- 7.08765e-05

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(4))
  expect_gfd(testfun, testfun$x0(10))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(10)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x4), rep(0, 4))
  expect_equal(testfun$gr(min_x10), rep(0, 10), tol = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x4), min_fx4)
  expect_equal(testfun$fn(min_x10), min_fx10)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0(4), fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS")
  expect_equal(res$par, min_x4, tol = 1e-3)
  expect_equal(res$value, min_fx4)

  res <- stats::optim(par = testfun$x0(10), fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS")
  expect_equal(res$par, min_x10, tol = 1e-6)
  expect_equal(res$value, min_fx10)
})
