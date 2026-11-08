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
  checks <- logical()
  expect_no_warning({
    for (name in names(integer_parameter_cases)) {
      problem <- get_problem_factory(name)()
      x <- integer_parameter_cases[[name]]
      names(x) <- paste0("x", seq_along(x))
      double_x <- x + 0
      for (callback in c("fn", "gr", "fg", "he")) {
        label <- paste(name, callback)
        expected <- problem[[callback]](double_x)
        actual <- problem[[callback]](x)
        checks[[paste(label, "double result is finite")]] <-
          all(is.finite(unlist(expected)))
        checks[[paste(label, "integer result matches")]] <-
          identical(actual, expected)
      }
    }
  })
  expect_all_checks(checks)
})

test_that("Rosenbrock callbacks reject nonnumeric parameters", {
  for (callback in c("fn", "gr", "fg", "he")) {
    expect_error(rosen()[[callback]](c("1", "2")), "non-numeric")
  }
})
