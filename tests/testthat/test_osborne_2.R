testfun <- osborne_2()
min_x <- c(
  1.30998,
  0.43155,
  0.63366,
  0.59943,
  0.75418,
  0.90429,
  1.36581,
  4.82370,
  2.39869,
  4.56888,
  5.67534
)
min_fx <- 0.04013774

test_that("Off-start derivatives match finite differences", {
  par <- c(1, 0.5, 0.6, 0.5, 0.8, 1.1, 1.2, 2, 1.8, 2.2, 2.5)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 11), tolerance = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tolerance = 1e-6)
})
