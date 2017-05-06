context("Linear Function - Rank 1 with Zero Columns and Rows")
n <- 10
m <- 15
testfun <- linfun_r1z(m)
min_x10 <- c(1, 0.69092332, 0.53638498, 0.38184664, 0.22730829, 0.07276995,
             -0.08176839, -0.23630673, -0.39084507, 1)
min_fx10 <- (m * m + 3 * m - 6) / (2 * (2 * m - 3))

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(n))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(n)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Minima has expected relationship to m", {
  expect_equal(sum(2:(n - 1) * min_x10[2:(n - 1)]),  3 / (2 * m - 3), tol = 1e-5)
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, n), tol = 1e-3)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tol = 1e-6)
})

test_that("Optimizer can reach minimum", {
  res <- stats::optim(par = testfun$x0(n), fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(factr = 0.1))
  expect_equal(res$par, min_x10, tol = 1e-6)
  expect_equal(res$value, min_fx10, tol = 1e-6)
})
