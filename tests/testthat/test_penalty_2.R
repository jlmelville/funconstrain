context("Penalty II")
testfun <- penalty_2()
min_x4 <- c(0.1999993, 0.1913145, 0.4800895, 0.5188700)
min_fx4 <- 9.37629e-06

min_x10 <- c(0.19998361, 0.01035098, 0.01960492, 0.03208906, 0.04993267,
             0.07651399, 0.11862407, 0.19214487, 0.34732059, 0.36916432)
min_fx10 <- 0.000293660 # actually only get to  0.0002936605

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(4))
  expect_gfd(testfun, testfun$x0(10))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(10)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x4), rep(0, 4), tol = 1e-6)
  expect_equal(testfun$gr(min_x10), rep(0, 10), tol = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x4), min_fx4, tol = 1e-6)
  expect_equal(testfun$fn(min_x10), min_fx10, tol = 1e-6)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0(4), fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(maxit = 1000,
                                                          factr = 1e-1))
  expect_equal(res$par, min_x4, tol = 1e-4)
  expect_equal(res$value, min_fx4, tol = 1e-6) # NOTE:220519 get fn not orig

  res <- stats::optim(par = testfun$x0(10), fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(maxit = 1000,
                                                          factr = 1e-1))
  expect_equal(res$par, min_x10, tol = 1e-3)
  expect_equal(res$value, min_fx10)
})
