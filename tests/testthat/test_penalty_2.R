testfun <- penalty_2()
min_x4 <- c(0.1999993, 0.1913145, 0.4800895, 0.5188700)
min_fx4 <- 9.37629e-06

min_x10 <- c(
  0.19998361,
  0.01035098,
  0.01960492,
  0.03208906,
  0.04993267,
  0.07651399,
  0.11862407,
  0.19214487,
  0.34732059,
  0.36916432
)
min_fx10 <- 0.000293660 # actually only get to  0.0002936605

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x4), rep(0, 4), tolerance = 1e-6)
  expect_equal(testfun$gr(min_x10), rep(0, 10), tolerance = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x4), min_fx4, tolerance = 1e-6)
  expect_equal(testfun$fn(min_x10), min_fx10, tolerance = 1e-5)
})
