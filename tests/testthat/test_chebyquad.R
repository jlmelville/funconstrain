context("Chebyquad")

min_fx <- 0 # for m = n = 1..7, 9
min_fx8 <- 3.51687e-3 # for m = n = 8
min_fx10 <- 6.50395e-3 # for m = n = 10

min_x1 <- c(0.5)
min_x2 <- c(0.2113249, 0.7886751)
min_x3 <- c(0.1464466, 0.5000000, 0.8535534)
min_x4 <- c(0.1026728, 0.4062038, 0.5937962, 0.8973272)
min_x5 <- c(0.08375126, 0.31272930, 0.50000000, 0.68727070, 0.91624874)
min_x6 <- c(0.06687659, 0.28874067, 0.36668230, 0.63331770, 0.71125933,
            0.93312341)
min_x7 <- c(0.05806915, 0.23517161, 0.33804409, 0.50000000, 0.66195591,
            0.76482839, 0.94193085)
min_x8 <- c(0.04315276, 0.19309084, 0.26632871, 0.50000000, 0.50000000,
            0.73367129, 0.80690916, 0.95684724)
min_x9 <- c(0.04420535, 0.19949067, 0.23561911, 0.41604691, 0.50000000,
            0.58395309, 0.76438089, 0.80050933, 0.95579465)
min_x10 <- c(0.0596199, 0.1667083, 0.2391707, 0.3988843, 0.3988843, 0.6011157,
             0.6011157, 0.7608293, 0.8332917, 0.9403801)


testfun <- chebyquad()


test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0(10))
})
test_that("f, g, and fg match at x0", {
  x0 <- testfun$x0(10)
  fg <- testfun$fg(x0)
  expect_equal(fg$fn, testfun$fn(x0))
  expect_equal(fg$gr, testfun$gr(x0))
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x1), rep(0, 1), tol = 1e-3)
  expect_equal(testfun$gr(min_x2), rep(0, 2), tol = 1e-3)
  expect_equal(testfun$gr(min_x3), rep(0, 3), tol = 1e-3)
  expect_equal(testfun$gr(min_x4), rep(0, 4), tol = 1e-3)
  expect_equal(testfun$gr(min_x5), rep(0, 5), tol = 1e-3)
  expect_equal(testfun$gr(min_x6), rep(0, 6), tol = 1e-3)
  expect_equal(testfun$gr(min_x7), rep(0, 7), tol = 1e-3)
  expect_equal(testfun$gr(min_x8), rep(0, 8), tol = 1e-3)
  expect_equal(testfun$gr(min_x9), rep(0, 9), tol = 1e-3)
  expect_equal(testfun$gr(min_x10), rep(0, 10), tol = 1e-3)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x1), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x2), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x3), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x4), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x5), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x6), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x7), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x8), min_fx8, tol = 1e-6)
  expect_equal(testfun$fn(min_x9), min_fx, tol = 1e-6)
  expect_equal(testfun$fn(min_x10), min_fx10, tol = 1e-6)
})

test_that("Optimizer can reach minimum", {
  res <- stats::optim(par = testfun$x0(10), fn = testfun$fn, gr = testfun$gr,
                      method = "L-BFGS-B", control = list(factr = 0.1))
  expect_equal(res$par, min_x10, tol = 1e-6)
  expect_equal(res$value, min_fx10, tol = 1e-6)
})
