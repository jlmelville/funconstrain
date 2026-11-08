test_that("Brown Hessian matches residual derivatives including zero products", {
  problem <- brown_al()
  cases <- list(
    negative_product = c(-0.5, 1, 1.5, 2),
    one_zero = c(0.5, 0, 1.5, 2),
    two_zeros = c(0, 1, 0, 2),
    more_than_two_zeros = c(0, 0, 0, 2),
    empty_product = c(0, 0)
  )
  for (case in names(cases)) {
    x <- cases[[case]]
    n <- length(x)
    linear_jacobian <- diag(n)[-n, , drop = FALSE] + 1
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
    expect_equal(problem$he(x), expected, tolerance = 1e-12, info = case)
  }
})
