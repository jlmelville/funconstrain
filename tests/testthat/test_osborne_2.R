context("Osborne 2")
testfun <- osborne_2()
min_x <- c(1.30998, 0.43155, 0.63366, 0.59943, 0.75418, 0.90429, 1.36581,
           4.82370, 2.39869, 4.56888, 5.67534)
min_fx <- 0.04013774

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
  expect_equal(gr0, rep(0, 11), tol = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0, fn = testfun$fn, gr = testfun$gr,
                      method = "BFGS", control = list(maxit = 1000))
  expect_equal(res$par, min_x, tol = 1e-5)
  expect_equal(res$value, min_fx)
})
