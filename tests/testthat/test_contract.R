test_that("problem factories expose documented core fields", {
  core_fields <- problem_factory_core_fields()
  factory_names <- problem_factory_names()
  actual <- lapply(factory_names, function(name) {
    testfun <- get_problem_factory(name)()
    x0 <- standard_x0(testfun)

    list(
      has_core_fields = all(core_fields %in% names(testfun)),
      fn_is_function = is.function(testfun$fn),
      gr_is_function = is.function(testfun$gr),
      he_is_function = is.function(testfun$he),
      fg_is_function = is.function(testfun$fg),
      fmin_is_numeric = is.numeric(testfun$fmin),
      fmin_is_scalar = length(testfun$fmin) == 1L,
      xmin_is_nonempty = length(testfun$xmin) > 0L,
      x0_is_numeric = is.numeric(x0),
      x0_is_nonempty = length(x0) > 0L
    )
  })
  names(actual) <- factory_names

  expect_all_checks(actual)
})

test_that("fixed-dimensional callbacks reject wrong-length parameters", {
  cases <- list(
    list(name = "rosen", n = 2L, problem = "Rosenbrock"),
    list(name = "freud_roth", n = 2L, problem = "Freudenstein-Roth"),
    list(name = "powell_bs", n = 2L, problem = "Powell Badly Scaled"),
    list(name = "brown_bs", n = 2L, problem = "Brown Badly Scaled"),
    list(name = "beale", n = 2L, problem = "Beale"),
    list(name = "jenn_samp", n = 2L, problem = "Jennrich-Sampson"),
    list(name = "helical", n = 3L, problem = "Helical Valley"),
    list(name = "bard", n = 3L, problem = "Bard"),
    list(name = "gauss", n = 3L, problem = "Gaussian"),
    list(name = "meyer", n = 3L, problem = "Meyer"),
    list(name = "gulf", n = 3L, problem = "Gulf"),
    list(name = "box_3d", n = 3L, problem = "Box 3D"),
    list(name = "powell_s", n = 4L, problem = "Powell Singular"),
    list(name = "wood", n = 4L, problem = "Wood"),
    list(name = "kow_osb", n = 4L, problem = "Kowalik-Osborne"),
    list(name = "brown_den", n = 4L, problem = "Brown Dennis"),
    list(name = "osborne_1", n = 5L, problem = "Osborne 1"),
    list(name = "biggs_exp6", n = 6L, problem = "Biggs EXP6"),
    list(name = "osborne_2", n = 11L, problem = "Osborne 2")
  )

  actual <- character()
  expected <- character()
  for (case in cases) {
    testfun <- get_problem_factory(case$name)()
    for (callback in c("fn", "gr", "he", "fg")) {
      for (bad_n in c(case$n - 1L, case$n + 1L)) {
        key <- paste(case$name, callback, "n =", bad_n)
        actual[[key]] <- capture_error_message(
          testfun[[callback]](rep(0, bad_n))
        )
        expected[[key]] <- paste0(
          case$problem,
          ": n is outside the allowed range"
        )
      }
    }
  }

  expect_identical(actual, expected)
})

test_that("reported xmin values evaluate to reported fmin values", {
  factory_names <- problem_factory_names()
  factories <- setNames(
    lapply(factory_names, function(name) {
      get_problem_factory(name)()
    }),
    factory_names
  )
  expect_true(anyNA(factories$box_3d$xmin))
  factories$box_3d <- NULL
  factories[["gulf m=10"]] <- gulf(m = 10)

  checks <- lapply(factories, function(testfun) {
    if (anyNA(testfun$xmin)) {
      return(c(xmin_is_available = FALSE))
    }
    objective <- testfun$fn(testfun$xmin)
    tolerance <- max(1e-8, abs(testfun$fmin) * 1e-5)
    c(
      objective_is_finite = is.finite(objective),
      objective_matches_fmin = abs(objective - testfun$fmin) <= tolerance
    )
  })

  expect_all_checks(checks)
})

test_that("linfun_r1z handles empty interiors at n=1 and n=2", {
  m_values <- c(1, 100)
  actual <- lapply(1:2, function(n) {
    testfun <- linfun_r1z(m_values[[n]])
    x <- seq_len(n)
    combined <- testfun$fg(x)

    list(
      fn = testfun$fn(x),
      gr = testfun$gr(x),
      fg_fn = combined$fn,
      fg_gr = combined$gr,
      he = testfun$he(x)
    )
  })
  names(actual) <- paste0("n=", 1:2)
  expected <- lapply(1:2, function(n) {
    m <- m_values[[n]]
    list(
      fn = m,
      gr = rep(0, n),
      fg_fn = m,
      fg_gr = rep(0, n),
      he = matrix(0, nrow = n, ncol = n)
    )
  })
  names(expected) <- names(actual)

  expect_equal(actual, expected)
})
