testfun <- biggs_exp6()
min_x <- c(1.71142, 17.6833, 1.16314, 5.18662, 1.71142, 1.16314)
min_fx <- 5.65565e-3

# Found by applying BFGS to the starting point
min_x2 <- c(1, 10, 1, 5, 4, 3)
min_fx2 <- 0

test_that("Off-start derivatives match finite differences", {
  par <- c(0.8, 1.4, 1.1, 0.9, 0.7, 1.2)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 6), tolerance = 1e-4)

  gr0_2 <- testfun$gr(min_x2)
  expect_equal(gr0_2, rep(0, 6))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tolerance = 1e-6)
  expect_equal(testfun$fn(min_x2), min_fx2)
})

# Which minimum is found seems to be dependent on m and the type of minimizer
# used, so may not be a good idea for a unit test in the long run
