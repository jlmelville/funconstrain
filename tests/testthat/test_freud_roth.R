testfun <- freud_roth()
test_that("Off-start derivatives match finite differences", {
  par <- c(0.7, -0.4)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(c(5, 4)), c(0, 0))
  expect_equal(testfun$gr(c(11.413, -0.8968)), c(0, 0), tolerance = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(5, 4)), 0, tolerance = 1e-3)
  expect_equal(testfun$fn(c(11.413, -0.8968)), 48.984, tolerance = 1e-5)
})
