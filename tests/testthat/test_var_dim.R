testfun <- var_dim()
min_x12 <- rep(1, 12)
min_fx12 <- 0

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x12), rep(0, 12), tolerance = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x12), min_fx12, tolerance = 1e-6)
})
