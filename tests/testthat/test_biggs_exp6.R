context("Biggs EXP6")
testfun <- biggs_exp6()
min_x <- c(1.71142, 17.6833, 1.16314, 5.18662, 1.71142, 1.16314)
min_fx <- 5.65565e-3

# Found by applying BFGS to the starting point
min_x2 <- c(1, 10, 1, 5, 4, 3)
min_fx2 <- 0

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 6), tol = 1e-4)

  gr0_2 <- testfun$gr(min_x2)
  expect_equal(gr0_2, rep(0, 6))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
  expect_equal(testfun$fn(min_x2), min_fx2)
})

# Which minimum is found seems to be dependent on m and the type of minimizer
# used, so may not be a good idea for a unit test in the long run
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B")
  expect_equal(res$par, min_x, tol = 1e-6)
  expect_equal(res$value, min_fx)
})
