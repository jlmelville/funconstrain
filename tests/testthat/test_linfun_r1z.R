n <- 10
m <- 15
testfun <- linfun_r1z(m)
min_x10 <- c(
  1,
  0.69092332,
  0.53638498,
  0.38184664,
  0.22730829,
  0.07276995,
  -0.08176839,
  -0.23630673,
  -0.39084507,
  1
)
min_fx10 <- (m * m + 3 * m - 6) / (2 * (2 * m - 3))

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Minima has expected relationship to m", {
  expect_equal(
    sum(2:(n - 1) * min_x10[2:(n - 1)]),
    3 / (2 * m - 3),
    tolerance = 1e-5
  )
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, n), tolerance = 1e-3)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tolerance = 1e-6)
})
