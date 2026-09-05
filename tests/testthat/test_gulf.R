testfun <- gulf()
test_that("Gulf endpoint limits include zero and nonzero curvature", {
  for (p in c(1, 1.5, 2)) {
    x <- c(50, 25, p)
    shorter <- gulf(99)
    full <- gulf(100)
    delta_h <- matrix(0, 3, 3)
    if (p == 1) delta_h[2, 2] <- 0.0008
    expect_warning(gradient <- full$gr(x), NA)
    expect_warning(combined <- full$fg(x), NA)
    expect_warning(hessian <- full$he(x), NA)
    expect_true(all(is.finite(c(gradient, unlist(combined), hessian))))
    expect_equal(full$fn(x), shorter$fn(x))
    expect_equal(gradient, shorter$gr(x))
    expect_equal(combined$fn, full$fn(x))
    expect_equal(combined$gr, gradient)
    expect_equal(hessian - shorter$he(x), delta_h, tolerance = 1e-10)
  }
})

test_that("Gulf interior zero distances retain negative or zero curvature", {
  t <- 0.5
  y <- 25 + (-50 * log(t))^(2 / 3)
  for (p in c(2, 3)) {
    x <- c(50, y, p)
    shorter <- gulf(49)
    full <- gulf(50)
    delta_h <- matrix(0, 3, 3)
    if (p == 2) delta_h[2, 2] <- -0.04
    expect_equal(full$fn(x) - shorter$fn(x), (1 - t)^2)
    expect_equal(full$gr(x), shorter$gr(x))
    expect_equal(full$fg(x)$gr, full$gr(x))
    expect_equal(full$fg(x)$fn, full$fn(x))
    expect_equal(full$he(x) - shorter$he(x), delta_h, tolerance = 1e-10)
    # At p=3 the third derivative jumps at zero. A small step limits the
    # O(step) bias of differentiating the gradient across that point.
    expect_hfd(full, x, tolerance = 1e-5, rel_eps = 1e-7)
  }
})

test_that("Gulf deliberately distinguishes undefined gradients and Hessians", {
  y <- 25 + (-50 * log(0.5))^(2 / 3)
  for (case in list(
    list(m = 50, y = y, p = 1.5),
    list(m = 100, y = 25, p = 0.75)
  )) {
    problem <- gulf(case$m)
    x <- c(50, case$y, case$p)
    expect_true(all(is.finite(problem$gr(x))))
    expect_equal(problem$fg(x)$gr, problem$gr(x))
    expect_equal(problem$fg(x)$fn, problem$fn(x))
    expect_error(problem$he(x), "Gulf: Hessian is undefined")
  }
  for (case in list(
    list(m = 50, y = y, p = 1),
    list(m = 100, y = 25, p = 0.5)
  )) {
    problem <- gulf(case$m)
    x <- c(50, case$y, case$p)
    expect_true(is.finite(problem$fn(x)))
    expect_error(problem$gr(x), "Gulf: gradient is undefined")
    expect_error(problem$fg(x), "Gulf: gradient is undefined")
    expect_error(problem$he(x), "Gulf: Hessian is undefined")
  }
})

test_that("Gulf uses regular derivatives on both sides of exact zero distance", {
  for (m in c(50, 100)) {
    t <- 0.01 * m
    y <- 25 + (-50 * log(t))^(2 / 3)
    for (p in c(1, 1.5, 2)) {
      for (delta in c(-1e-4, -1e-7, 1e-7, 1e-4)) {
        x <- c(50, y + delta, p)
        distance <- x[2] - y
        u <- abs(distance)^p / x[1]
        e <- exp(-u)
        residual <- expm1(-u) + (1 - t)
        du <- p * sign(distance) * abs(distance)^(p - 1) / x[1]
        ddu <- p * (p - 1) * abs(distance)^(p - 2) / x[1]
        expected_g <- -2 * e * residual * du
        expected_h <- 2 * e * ((e + residual) * du^2 - residual * ddu)
        full <- gulf(m)
        shorter <- gulf(m - 1)
        actual_g <- full$gr(x)[2] - shorter$gr(x)[2]
        actual_h <- full$he(x)[2, 2] - shorter$he(x)[2, 2]
        expect_lt(abs(actual_g - expected_g), 1e-11 * max(1, abs(expected_g)))
        expect_lt(abs(actual_h - expected_h), 1e-11 * max(1, abs(expected_h)))
        expect_equal(full$fg(x)$gr, full$gr(x))
        expect_equal(full$fg(x)$fn, full$fn(x))
      }
    }
  }
})
min_x <- c(50, 25, 1.5)
min_fx <- 0
test_that("Analytical and finite difference gradients match at x0", {
  expect_gfd(testfun, testfun$x0)
})
test_that("f, g, and fg match at x0", {
  fg <- testfun$fg(testfun$x0)
  expect_equal(fg$fn, testfun$fn(testfun$x0))
  expect_equal(fg$gr, testfun$gr(testfun$x0))
})
test_that("Off-start derivatives match finite differences", {
  par <- c(10, 5, 1)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Gradient is zero at stated minima", {
  gr0 <- testfun$gr(min_x)
  expect_equal(gr0, c(0, 0, 0))
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x), min_fx)
})
test_that("Optimizer can reach minimum from x0", {
  res <- stats::optim(
    par = testfun$x0,
    fn = testfun$fn,
    gr = testfun$gr,
    method = "L-BFGS-B"
  )
  expect_equal(res$par, min_x, tolerance = 1e-6)
  expect_equal(res$value, min_fx)
})
