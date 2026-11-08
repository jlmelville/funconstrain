testfun <- jenn_samp(m = 10)
test_that("Off-start derivatives match finite differences", {
  par <- c(0.1, 0.2)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(rep(0.2578252, 2))
  expect_equal(gr0, rep(0, 2), tolerance = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(rep(0.2578252, 2)), 124.362, tolerance = 1e-5)
})
