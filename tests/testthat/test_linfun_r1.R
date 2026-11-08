n <- 10
m <- 15
testfun <- linfun_r1(m)
min_x10 <- c(
  0.857394219,
  0.714788437,
  0.572182656,
  0.429576875,
  0.286971093,
  0.144365312,
  0.001759531,
  -0.140846251,
  -0.283452032,
  -0.426057813
)
min_fx10 <- (m * (m - 1)) / (2 * (2 * m + 1))

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Minima has expected relationship to m", {
  expect_equal(sum(1:n * min_x10), 3 / (2 * m + 1), tolerance = 1e-5)
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, n), tolerance = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tolerance = 1e-6)
})
