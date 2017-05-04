context("Watson")
testfun <- watson()
min_x <- c(-0.01573919, 1.01241881, -0.23273987, 1.25964331, -1.51285698,
           0.99265101)
min_fx <- 0.002287674

min_x9 <- c(-1.554337e-05,  9.997932e-01,  1.460705e-02,  1.478749e-01,
             9.943048e-01, -2.603397e+00,  4.087345e+00, -3.133179e+00,
             1.050054e+00)
min_fx9 <- 1.39976e-6 # Couldn't get this low 1.3999e-6

min_x12 <- c(-1.586076e-07,  1.000015e+00, -1.028239e-03,  3.457001e-01,
            -5.631347e-02,  2.263312e-01,  5.575928e-02, -2.931759e-01,
            2.301305e-01, 3.888203e-01, -5.884314e-01,  2.495999e-01)
min_fx12 <- 4.72238e-10 # Couldn't get this low: 1.0295e-8

test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(6))
  expect_gfd(testfun, testfun$x0(9))
  expect_gfd(testfun, testfun$x0(12))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(6)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x), rep(0, 6), tol = 1e-4)
  expect_equal(testfun$gr(min_x9), rep(0, 9), tol = 1e-4)
  expect_equal(testfun$gr(min_x12), rep(0, 12), tol = 1e-4)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
  expect_equal(testfun$fn(min_x9), min_fx9)
  expect_equal(testfun$fn(min_x12), min_fx12)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(par = testfun$x0(6), fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B")
  expect_equal(res$par, min_x, tol = 1e-5)
  expect_equal(res$value, min_fx)
})
