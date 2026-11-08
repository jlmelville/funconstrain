testfun <- powell_bs()
min_x <- c(1.098177e-5, 9.106)
min_fx <- 0
test_that("Gradient matches finite differences when the first coordinate is zero", {
  expect_gfd(testfun, testfun$x0)
})
test_that("Off-start derivatives match finite differences", {
  par <- c(0.02, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, c(0, 0), tolerance = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tolerance = 1e-6)
})
