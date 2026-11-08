testfun <- disc_bv()
min_x10 <- c(
  -0.04316498,
  -0.08157716,
  -0.11448571,
  -0.14097358,
  -0.15990870,
  -0.16987720,
  -0.16908998,
  -0.15524954,
  -0.12535589,
  -0.07541653
)
min_fx10 <- 0

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, 10), tolerance = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tolerance = 1e-6)
})

test_that("Standard starting point matches the canonical definition", {
  n <- 4
  t <- seq_len(n) / (n + 1)
  expect_equal(testfun$x0(n), t * (t - 1))
})
