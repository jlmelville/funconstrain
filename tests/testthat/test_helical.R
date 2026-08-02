testfun <- helical()
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0, tolerance = 1e-3)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(c(1, 0, 0))
  expect_equal(gr0, c(0, 0, 0))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(c(1, 0, 0)), 0)
})

test_that("Gradient is finite on the smooth positive y-axis", {
  cases <- list(
    axis_z0 = list(
      par = c(0, 1, 0),
      expected = c(-2500 / pi, 0, -500)
    ),
    axis_z2 = list(
      par = c(0, 1, 2),
      expected = c(-500 / pi, 0, -96)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    actual <- testfun$gr(case$par)
    fg <- testfun$fg(case$par)

    expect_equal(actual, case$expected, tolerance = 1e-12)
    expect_true(all(is.finite(actual)), info = case_name)
    expect_equal(fg$gr, actual)
    expect_equal(fg$fn, testfun$fn(case$par))
    expect_true(all(is.finite(testfun$he(case$par))), info = case_name)
    expect_gfd(testfun, case$par, tolerance = 1e-3)
  }
})

test_that("Gradient remains correct at an ordinary off-axis point", {
  par <- c(1, 1, 0)
  angular_radial <- 200 * (sqrt(2) - 1) / sqrt(2)
  expected <- c(
    -625 / pi + angular_radial,
    625 / pi + angular_radial,
    -250
  )
  fg <- testfun$fg(par)

  expect_equal(testfun$gr(par), expected, tolerance = 1e-12)
  expect_equal(fg$gr, testfun$gr(par))
  expect_equal(fg$fn, testfun$fn(par))
  expect_gfd(testfun, par, tolerance = 1e-6)
})

test_that("Undefined Helical domains produce clear callback errors", {
  origin <- c(0, 0, 0)
  negative_y_axis <- c(0, -1, 0)

  expect_error(testfun$fn(origin), "undefined")
  expect_true(is.finite(testfun$fn(negative_y_axis)))

  for (par in list(origin, negative_y_axis)) {
    expect_error(testfun$gr(par), "undefined")
    expect_error(testfun$fg(par), "undefined")
    expect_error(testfun$he(par), "undefined")
  }
})

test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(
    par = testfun$x0,
    fn = testfun$fn,
    gr = testfun$gr,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  expect_equal(res$par, c(1, 0, 0))
  expect_equal(res$value, 0)
})
