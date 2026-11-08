testfun <- gauss()
test_that("Off-start derivatives match finite differences", {
  par <- c(0.5, 1.2, 0.1)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(0.398956, 1.000019, 0.0))
  expect_equal(gr0, c(0, 0, 0), tolerance = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(0.398956, 1.000019, 0.0)), 1.12793e-8)
})
