testfun <- beale()
test_that("Beale Hessian is defined when the second coordinate is zero", {
  for (a in c(-2, 0, 3)) {
    par <- c(a, 0)
    expected <- matrix(c(6, 3 - 4 * a, 3 - 4 * a, 9 * a - 2 * a^2), 2)
    expect_equal(testfun$he(par), expected)
  }
  for (b in c(-1e-5, 1e-5)) {
    expect_hfd(testfun, c(-2, b), tolerance = 1e-6)
  }
})
test_that("Off-start derivatives match finite differences", {
  par <- c(0.8, 0.3)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(3, 0.5))
  expect_equal(gr0, c(0, 0))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(3, 0.5)), 0, tolerance = 1e-6)
})
