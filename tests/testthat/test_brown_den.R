testfun <- brown_den()
min_x <- c(-11.5944399, 13.2036301, -0.4034394, 0.2367788)
min_fx <- 85822.2

test_that("Off-start derivatives match finite differences", {
  par <- c(1, 0.5, 0.2, -0.3)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 4), tolerance = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tolerance = 1.e-2)
})
