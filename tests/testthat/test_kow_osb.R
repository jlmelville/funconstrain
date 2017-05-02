context("Kowalik and Osborne")
testfun <- kow_osb()
min_x <- c(0.1928, 0.1914, 0.1231, 0.1361) # only gets to 3.05706e-4
min_fx <- 3.07505e-4

# given as c(Inf, -14.07..., -Inf, -Inf) in MGH paper
min_x2 <- c(4.8905e3, -14.072, -2.21246e5, -1.380389e5)
min_fx2 <- 1.02734e-3 # Actually min_x2 only gets to 1.027346e-3, close enough!

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, rep(0, 4), tol = 1e-5)
  expect_equal(testfun$gr(min_x2), rep(0, 4))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
  expect_equal(testfun$fn(min_x2), min_fx2)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS",
                      control = list(reltol = 1e-6, abstol = 1e-6))
  expect_equal(res$par, min_x, tol = 1e-4)
  expect_equal(res$value, min_fx)
})
