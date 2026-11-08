testfun <- osborne_1()
min_x <- c(0.375410, 1.935836, -1.464676, 0.0128676, 0.0221227)
min_fx <- 5.46489e-05 # actually 5.464895e-05

test_that("Off-start derivatives match finite differences", {
  par <- c(0.4, 1.1, -0.8, 0.03, 0.015)
  expect_gfd(testfun, par, tolerance = 1e-5)
  # A smaller step resolves the stiff exponential curvature more accurately.
  expect_hfd(
    testfun,
    par,
    tolerance = 1e-5,
    rel_eps = 1e-5
  )
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 5), tolerance = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tolerance = 1e-5)
})
