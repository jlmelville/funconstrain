context("Linear Function - Rank 1")
n <- 10
m <- 15
testfun <- linfun_r1(m)
min_x10 <- c(0.857394219, 0.714788437, 0.572182656, 0.429576875, 0.286971093,
             0.144365312, 0.001759531, -0.140846251, -0.283452032, -0.426057813)
min_fx10 <- (m * (m - 1)) / (2 * (2 * m + 1))

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(10))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(n)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Minima has expected relationship to m", {
  expect_equal(sum(1:n * min_x10),  3 / (2 * m + 1), tol = 1e-5)
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x10), rep(0, n), tol = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x10), min_fx10, tol = 1e-6)
})

test_that("Optimizer can reach minimum", {
  res <- stats::optim(par = testfun$x0(n), fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B")
  expect_equal(res$par, min_x10, tol = 1e-6)
  expect_equal(res$value, min_fx10, tol = 1e-6)
})
