testfun <- broyden_tri()
min_x10 <- c(
  -0.5707221,
  -0.6818070,
  -0.7022101,
  -0.7055106,
  -0.7049062,
  -0.7014966,
  -0.6918893,
  -0.6657965,
  -0.5960351,
  -0.4164123
)
min_fx10 <- 0

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, 10), tolerance = 1e-5)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tolerance = 1e-6)
})
