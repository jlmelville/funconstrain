context("Discrete Boundary Value")
testfun <- disc_bv()
min_x10 <- c(-0.04316498, -0.08157716, -0.11448571, -0.14097358, -0.15990870,
             -0.16987720, -0.16908998, -0.15524954, -0.12535589, -0.07541653)
min_fx10 <- 0

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(10))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(10)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, 10), tol = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tol = 1e-6)
})

test_that("Optimizer can reach minimum", {
  res <- stats::optim(par = testfun$x0(10), fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS")
  expect_equal(res$par, min_x10, tol = 1e-6)
  expect_equal(res$value, min_fx10, tol = 1e-6)
})
