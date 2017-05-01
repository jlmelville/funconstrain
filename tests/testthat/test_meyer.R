# This function is a pain to minimize via optim
# BFGS does best but gives up after ~500 iterations

context("Meyer")
testfun <- meyer()
min_x <- c(5.60963766837e-03, 6.18134616768e+03, 3.45223628592e+02)
min_fx <- 87.9458 # can only reach 87.94586
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
  expect_equal(gr0[1], 0, tol = 1) # -0.65!
  expect_equal(gr0[2], 0, tol = 1e-3) # 8.8e-4
  expect_equal(gr0[3], 0, tol = 1e-1) # 0.013
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tol = 1e-5)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS", control = list(maxit = 1000))
  expect_equal(res$par, min_x, tol = 1e-4)
  expect_equal(res$value, min_fx, tol = 1e-5)
})
