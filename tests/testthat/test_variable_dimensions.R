test_that("low-dimensional Hessians match finite differences", {
  cases <- data.frame(
    factory = c(
      "broyden_tri",
      "chebyquad",
      "disc_bv",
      "disc_ie",
      "disc_ie"
    ),
    n = c(2, 1, 2, 2, 10),
    tolerance = c(1e-5, 1e-8, 1e-5, 1e-5, 1e-5),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    testfun <- get_problem_factory(case$factory)()
    x0 <- testfun$x0(case$n)
    fg <- testfun$fg(x0)
    hessian <- testfun$he(x0)

    expect_equal(length(testfun$fn(x0)), 1, info = case$factory)
    expect_equal(length(testfun$gr(x0)), case$n, info = case$factory)
    expect_equal(length(fg$gr), case$n, info = case$factory)
    expect_equal(fg$fn, testfun$fn(x0), info = paste(case$factory, "fg fn"))
    expect_equal(fg$gr, testfun$gr(x0), info = paste(case$factory, "fg gr"))
    expect_true(is.matrix(hessian), info = case$factory)
    expect_equal(dim(hessian), c(case$n, case$n), info = case$factory)
    expect_true(isSymmetric(hessian), info = case$factory)
    expect_equal(
      hessian,
      hessian_fd(x0, testfun$gr),
      tolerance = case$tolerance,
      info = paste(case$factory, "n =", case$n)
    )
  }
})
