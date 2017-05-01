context("Powell Badly Scaled")

testfun <- powell_bs()
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(1.098e-5, 9.106))
  expect_equal(gr0[1], 0, tol = 1) # gr is ~ -29.4!
  expect_equal(gr0[2], 0, tol = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(1.098e-5, 9.106)), 0, tol = 1e-6)
})
test_that("Optimizer can reach minimum from x0", {
  # Really need to crank up the number of iterations to get evenly remotely
  # close
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS", control = list(maxit = 10000))
  expect_equal(res$par[1], 1.098e-5, tol = 1e-5) # ~ 1.542e-5
  expect_equal(res$par[2], 9.106, tol = 1) # ~6.483
  expect_equal(res$value, 0, tol = 1e-5) # ~2e-6
})
