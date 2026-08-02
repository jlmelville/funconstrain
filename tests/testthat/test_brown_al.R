testfun <- brown_al()
min_x12 <- rep(1, 12)
min_fx12 <- 0
min_x12b <- c(rep(0, 11), 13)
min_fx12b <- 1

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
  expect_equal(testfun$gr(min_x12), rep(0, 12), tolerance = 1e-6)
  expect_equal(testfun$gr(min_x12b), rep(0, 12), tolerance = 1e-6)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x12), min_fx12, tolerance = 1e-6)
  expect_equal(testfun$fn(min_x12b), min_fx12b, tolerance = 1e-6)
})

test_that("Product residual gradient handles negative products and zeros", {
  cases <- list(
    negative_product = list(
      par = c(-1, 2, 3),
      expected = c(-84, 48, 30)
    ),
    one_zero = list(
      par = c(0, 2, 3),
      expected = c(-2, 14, 8)
    ),
    multiple_zeros = list(
      par = c(0, 0, 3),
      expected = c(-6, -6, -4)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    actual <- testfun$gr(case$par)
    fg <- testfun$fg(case$par)

    expect_equal(actual, case$expected)
    expect_equal(fg$gr, actual)
    expect_equal(fg$fn, testfun$fn(case$par))
    expect_true(all(is.finite(actual)), info = case_name)
    expect_true(all(is.finite(fg$gr)), info = case_name)
    expect_gfd(testfun, case$par, tolerance = 1e-6)
  }
})

test_that("Optimizer can reach minimum", {
  res <- stats::optim(
    par = testfun$x0(12),
    fn = testfun$fn,
    gr = testfun$gr,
    method = "BFGS"
  )
  expect_equal(res$par, min_x12, tolerance = 1e-6)
  expect_equal(res$value, min_fx12, tolerance = 1e-6)
})
