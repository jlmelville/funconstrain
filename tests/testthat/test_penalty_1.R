testfun <- penalty_1()
min_x4 <- rep(0.2500075, 4)
min_fx4 <- 2.24997e-05 # actually only reach 2.249978e-05

min_x10 <- rep(0.158122, 10)
min_fx10 <- 7.08765e-05

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x4), rep(0, 4))
  expect_equal(testfun$gr(min_x10), rep(0, 10), tolerance = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x4), min_fx4, tolerance = 1e-5)
  expect_equal(testfun$fn(min_x10), min_fx10, tolerance = 1e-5)
})
