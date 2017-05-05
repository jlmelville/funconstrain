context("Powell Badly Scaled")

testfun <- powell_bs()
min_x <- c(1.098177e-5, 9.106)
min_fx <- 0
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
  expect_equal(gr0, c(0, 0), tol = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx, tol = 1e-6)
})
test_that("Optimizer can reach minimum from x0", {
  # bizarrely, a memory = 1 does a really good job!
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(maxit = 1000, lmm = 1,
                                                          factr = 1e-6))
  expect_equal(res$par, min_x, tol = 1e-3)
  expect_equal(res$value, 0)
})
