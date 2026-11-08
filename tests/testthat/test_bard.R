testfun <- bard()
test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.8, 1.5)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Hessian singular branch raises an error", {
  expect_error(
    testfun$he(c(1, 0, 0)),
    regexp = "^Hessian is undefined at this parameter value$"
  )
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(0.08241, 1.13304, 2.34370))
  expect_equal(gr0, c(0, 0, 0), tolerance = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(
    testfun$fn(c(0.08241, 1.13304, 2.34370)),
    8.21487e-3,
    tolerance = 1e-6
  )
  expect_equal(testfun$fn(c(0.8406, 1e6, 1e6)), 17.4286, tolerance = 1e-5)
})
