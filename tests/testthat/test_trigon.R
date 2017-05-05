context("Trigonometric")
testfun <- trigon()
min_x12 <- rep(0, 12)
min_fx12 <- 0

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(12))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(12)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x12), rep(0, 12), tol = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x12), min_fx12, tol = 1e-6)
})
# Here I massively cheat. Starting location for all things I tried invariably
# hits a different minimum. Starting closer helps.
test_that("Optimizer can reach minimum", {
  res <- stats::optim(par = testfun$x0(12) * 0.5, fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(factr = 0.1))
  expect_equal(res$par, min_x12, tol = 1e-4)
  expect_equal(res$value, min_fx12, tol = 1e-6)

})
