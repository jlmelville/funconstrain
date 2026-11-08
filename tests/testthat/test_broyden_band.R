testfun <- broyden_band()
min_x10 <- c(
  -0.4283029,
  -0.4765964,
  -0.5196525,
  -0.5580993,
  -0.5925062,
  -0.6245037,
  -0.6232395,
  -0.6213938,
  -0.6204536,
  -0.5864693
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
