context("Osborne 1")
testfun <- osborne_1()
min_x <- c(0.375410, 1.935836, -1.464676, 0.0128676, 0.0221227)
min_fx <- 5.46489e-05 # actually 5.464895e-05

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
  expect_equal(gr0, rep(0, 5), tol = 1e-2)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(maxit = 1000,
                                                          factr = 0.1))
  expect_equal(res$par, min_x, tol = 1e-5)
  expect_equal(res$value, min_fx)
})
