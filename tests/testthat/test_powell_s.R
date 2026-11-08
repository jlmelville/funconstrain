testfun <- powell_s()
min_x <- rep(0, 4)
min_fx <- 0
test_that("Off-start derivatives match finite differences", {
  par <- c(0.3, -0.2, 0.1, -0.4)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 4))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
})
