# Domain-safe integer fixtures: large polynomial inputs expose early integer
# multiplication, while exponential/rational problems use nonsingular values.
integer_parameter_cases <- list(
  rosen = c(50000L, 1L),
  freud_roth = c(1L, 50000L),
  powell_bs = c(216L, 50000L),
  brown_bs = c(50000L, 2L),
  beale = c(1L, 50000L),
  jenn_samp = c(1L, 2L),
  helical = c(50000L, 1L, 2L),
  bard = c(216L, 1L, 2L),
  gauss = c(216L, 1L, 2L),
  meyer = c(216L, 50000L, 50000L),
  gulf = c(50L, 25L, 2L),
  box_3d = c(1L, 2L, 3L),
  powell_s = c(50000L, 1L, 2L, 3L),
  wood = c(50000L, 1L, 2L, 3L),
  kow_osb = c(1L, 2L, 3L, 4L),
  brown_den = c(50000L, 1L, 2L, 3L),
  osborne_1 = c(1L, 2L, 3L, 1L, 2L),
  biggs_exp6 = c(1L, 2L, 3L, 4L, 5L, 6L),
  osborne_2 = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 1L, 2L, 3L, 4L),
  watson = rep(c(50000L, 1L, 2L, 3L), 3),
  ex_rosen = rep(c(50000L, 1L, 2L, 3L), 2),
  ex_powell = rep(c(50000L, 1L, 2L, 3L), 2),
  penalty_1 = c(50000L, 1L, 2L, 3L),
  penalty_2 = c(216L, 1L, 2L, 3L),
  var_dim = c(50000L, 1L, 2L, 3L),
  trigon = c(50000L, 1L, 2L, 3L),
  brown_al = c(216L, 1L, 2L, 3L),
  disc_bv = c(50000L, 1L, 2L, 3L),
  disc_ie = c(50000L, 1L, 2L, 3L),
  broyden_tri = c(50000L, 1L, 2L, 3L),
  broyden_band = c(50000L, 1L, 2L, 3L),
  linfun_fr = c(216L, 1L, 2L, 3L),
  linfun_r1 = c(216L, 1L, 2L, 3L),
  linfun_r1z = c(216L, 1L, 2L, 3L),
  chebyquad = c(0L, 1L)
)

test_that("all factories promote integer arithmetic consistently in every callback", {
  expect_setequal(names(integer_parameter_cases), problem_factory_names())
  for (name in names(integer_parameter_cases)) {
    problem <- get_problem_factory(name)()
    for (named in c(FALSE, TRUE)) {
      x <- integer_parameter_cases[[name]]
      if (named) names(x) <- paste0("x", seq_along(x))
      double_x <- x + 0
      for (callback in c("fn", "gr", "fg", "he")) {
        label <- paste(name, callback, if (named) "named" else "unnamed")
        expected <- problem[[callback]](double_x)
        expect_true(all(is.finite(unlist(expected))), info = label)
        expect_warning(actual <- problem[[callback]](x), NA)
        expect_equal(actual, expected, tolerance = 0, info = label)
      }
      combined <- problem$fg(x)
      expect_equal(combined$fn, problem$fn(x))
      expect_equal(combined$gr, problem$gr(x))
    }
  }
})

test_that("integer promotion preserves attributes and does not coerce other types", {
  # The helper's identity contract protects arbitrary input classes without
  # introducing assumptions about their arithmetic methods into factory tests.
  promote <- getFromNamespace("promote_integer", "funconstrain")
  x <- structure(
    matrix(1:4, 2, dimnames = list(c("a", "b"), c("c", "d"))),
    note = "retained"
  )
  actual <- promote(x)
  expect_type(actual, "double")
  expect_identical(attributes(actual), attributes(x))
  expect_identical(unname(as.vector(actual)), as.double(x))
  for (other in list(
    c(a = 1, b = 2),
    c(1 + 1i, 2),
    c(TRUE, FALSE),
    c("1", "2"),
    factor(c("1", "2")),
    structure(1:2, class = "custom_integer")
  )) {
    expect_identical(promote(other), other)
  }
  for (callback in c("fn", "gr", "fg", "he")) {
    expect_error(rosen()[[callback]](c("1", "2")), "non-numeric")
  }
})
