test_that("integral accumulations agree with the dense residual Jacobian", {
  problem <- disc_ie()
  for (n in c(1, 2, 4, 10, 20)) {
    t <- seq_len(n) / (n + 1)
    h <- 1 / (n + 1)
    kernel <- outer(t, t, function(a, b) pmin(a, b) * (1 - pmax(a, b)))
    for (x in list(rep(0, n), rep(-1, n), seq_len(n) / n - 0.5)) {
      shifted <- x + t + 1
      residual <- x + as.vector(0.5 * h * kernel %*% shifted^3)
      jacobian <- diag(n) + 1.5 * h * kernel %*% diag(shifted^2, n)
      expected_fn <- sum(residual^2)
      expected_gr <- as.vector(2 * crossprod(jacobian, residual))
      expect_equal(problem$fn(x), expected_fn, tolerance = 1e-12)
      expect_equal(problem$gr(x), expected_gr, tolerance = 1e-12)
      expect_equal(
        problem$fg(x),
        list(fn = expected_fn, gr = expected_gr),
        tolerance = 1e-12
      )
      expect_gfd(problem, x, tolerance = 1e-6)
      expect_hfd(problem, x, tolerance = 1e-5)
    }
  }
})

test_that("integral accumulations retain unnamed vector gradients", {
  problem <- disc_ie()
  x <- c(a = -0.2, b = 0.4)
  expect_identical(problem$gr(x), problem$gr(unname(x)))
  expected_fn <- setNames(problem$fn(unname(x)), "a")
  expect_identical(problem$fn(x), expected_fn)
  expect_identical(problem$fg(x), list(fn = expected_fn, gr = problem$gr(x)))
  expect_equal(problem$gr(matrix(x)), problem$gr(unname(x)))
})
