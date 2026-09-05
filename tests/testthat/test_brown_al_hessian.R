test_that("Brown Hessian matches residual derivatives including zero products", {
  problem <- brown_al()
  for (n in c(1, 2, 4, 10, 20)) {
    for (kind in c("positive", "negative", "one_zero", "two_zeros", "zeros")) {
      x <- 0.5 + seq_len(n) / n
      if (kind == "negative") x[1] <- -x[1]
      if (kind == "one_zero") x[1] <- 0
      if (kind == "two_zeros") x[seq_len(min(2, n))] <- 0
      if (kind == "zeros") x[] <- 0
      linear_jacobian <- matrix(1, n - 1, n)
      if (n > 1) {
        for (i in seq_len(n - 1)) linear_jacobian[i, i] <- 2
      }
      product_gradient <- vapply(
        seq_len(n),
        function(i) prod(x[-i]),
        numeric(1)
      )
      product_hessian <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i != j) product_hessian[i, j] <- prod(x[-c(i, j)])
        }
      }
      expected <- 2 *
        crossprod(linear_jacobian) +
        2 * tcrossprod(product_gradient) +
        2 * (prod(x) - 1) * product_hessian
      actual <- problem$he(x)
      expect_true(all(is.finite(actual)), info = paste(n, kind))
      expect_equal(actual, expected, tolerance = 1e-12)
      expect_hfd(problem, x, tolerance = 1e-5)
    }
  }
})
