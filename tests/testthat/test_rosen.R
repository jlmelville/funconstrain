testfun <- rosen()
test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.7)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Hessian is correct at x0", {
  expect_equal(testfun$he(testfun$x0), matrix(c(1330, 480, 480, 200), 2))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(c(1, 1)), c(0, 0))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(1, 1)), 0)
})
