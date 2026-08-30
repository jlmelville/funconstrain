test_that("factory contract list covers all exported problem factories", {
  exported_factories <- setdiff(
    getNamespaceExports("funconstrain"),
    c(
      "fufn",
      "fufnrun",
      "funconstrain_catalog",
      "funconstrain_problem"
    )
  )

  actual <- list(
    names = sort(problem_factory_names()),
    count = length(problem_factory_names())
  )
  expected <- list(names = sort(exported_factories), count = 35L)

  expect_identical(actual, expected)
})

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
  checks <- logical()
  for (name in problem_factory_names()) {
    testfun <- get_problem_factory(name)()

    if (anyNA(testfun$xmin)) {
      next
    }

    objective <- testfun$fn(testfun$xmin)
    tolerance <- max(1e-8, abs(testfun$fmin) * 1e-5)

    checks[[paste(name, "finite")]] <- is.finite(objective)
    checks[[paste(name, "matches fmin")]] <-
      abs(objective - testfun$fmin) <= tolerance
  }

  expect_all_checks(checks)
})

test_that("documented stored reference configurations are internally valid", {
  cases <- list(
    list(
      name = "jenn_samp",
      args = list(m = 10),
      configuration = "m=10",
      xmin_length = 2L
    ),
    list(
      name = "gulf",
      args = list(m = 10),
      configuration = "m=10",
      xmin_length = 3L
    ),
    list(
      name = "box_3d",
      args = list(m = 20),
      configuration = "m=20",
      xmin_length = 3L
    ),
    list(
      name = "brown_den",
      args = list(m = 20),
      configuration = "m=20",
      xmin_length = 4L
    ),
    list(
      name = "biggs_exp6",
      args = list(m = 13),
      configuration = "m=13",
      xmin_length = 6L
    ),
    list(
      name = "watson",
      args = list(),
      configuration = "n=6",
      xmin_length = 6L
    ),
    list(
      name = "ex_rosen",
      args = list(),
      configuration = "n=8",
      xmin_length = 8L
    ),
    list(
      name = "ex_powell",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "penalty_1",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "penalty_2",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "var_dim",
      args = list(),
      configuration = "n=6",
      xmin_length = 6L
    ),
    list(
      name = "trigon",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "brown_al",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "disc_bv",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "disc_ie",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "broyden_tri",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "broyden_band",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "linfun_fr",
      args = list(m = 100),
      configuration = "n=4,m=100",
      xmin_length = 4L
    ),
    list(
      name = "linfun_r1",
      args = list(m = 100),
      configuration = "n=5,m=100",
      xmin_length = 5L
    ),
    list(
      name = "linfun_r1z",
      args = list(m = 100),
      configuration = "n=5,m=100",
      xmin_length = 5L
    ),
    list(
      name = "chebyquad",
      args = list(),
      configuration = "n=8",
      xmin_length = 8L
    )
  )

  actual <- lapply(cases, function(case) {
    testfun <- do.call(get_problem_factory(case$name), case$args)
    xmin_available <- !anyNA(testfun$xmin)
    objective <- if (xmin_available) testfun$fn(testfun$xmin) else NA_real_
    tolerance <- max(1e-8, abs(testfun$fmin) * 1e-5)

    list(
      xmin_length = length(testfun$xmin),
      xmin_available = xmin_available,
      objective_is_finite = if (xmin_available) is.finite(objective) else NA,
      objective_matches_fmin = if (xmin_available) {
        abs(objective - testfun$fmin) <= tolerance
      } else {
        NA
      }
    )
  })
  names(actual) <- vapply(
    cases,
    function(case) paste(case$name, case$configuration),
    character(1L)
  )
  expected <- lapply(cases, function(case) {
    xmin_available <- case$name != "box_3d"
    list(
      xmin_length = case$xmin_length,
      xmin_available = xmin_available,
      objective_is_finite = if (xmin_available) TRUE else NA,
      objective_matches_fmin = if (xmin_available) TRUE else NA
    )
  })
  names(expected) <- names(actual)

  expect_identical(actual, expected)
})

test_that("linfun_r1z handles empty interiors at n=1 and n=2", {
  testfun <- linfun_r1z(m = 100)
  actual <- lapply(1:2, function(n) {
    x <- rep(0, n)
    combined <- testfun$fg(x)

    list(
      fn = testfun$fn(x),
      gr = testfun$gr(x),
      fg_fn = combined$fn,
      fg_gr = combined$gr,
      he = testfun$he(x),
      stored_fmin_is_inapplicable = abs(testfun$fn(x) - testfun$fmin) > 1
    )
  })
  names(actual) <- paste0("n=", 1:2)
  expected <- lapply(1:2, function(n) {
    list(
      fn = 100,
      gr = rep(0, n),
      fg_fn = 100,
      fg_gr = rep(0, n),
      he = matrix(0, nrow = n, ncol = n),
      stored_fmin_is_inapplicable = TRUE
    )
  })
  names(expected) <- names(actual)

  expect_equal(actual, expected)
})

test_that("current m metadata is recorded without making it a core field", {
  factories_with_na_m <- c(
    "beale",
    "brown_bs",
    "chebyquad",
    "freud_roth",
    "jenn_samp",
    "powell_bs",
    "rosen"
  )
  factories_with_numeric_m <- c(
    linfun_fr = 100,
    linfun_r1 = 100,
    linfun_r1z = 100,
    trigon = 30
  )
  expected <- data.frame(
    factory = problem_factory_names(),
    m_present = FALSE,
    m_type = NA_character_,
    m_value = NA_character_,
    stringsAsFactors = FALSE
  )
  na_rows <- expected$factory %in% factories_with_na_m
  numeric_rows <- expected$factory %in% names(factories_with_numeric_m)

  expected$m_present[na_rows | numeric_rows] <- TRUE
  expected$m_type[na_rows] <- "logical"
  expected$m_value[na_rows] <- "NA"
  expected$m_type[numeric_rows] <- "double"
  expected$m_value[numeric_rows] <- as.character(
    factories_with_numeric_m[expected$factory[numeric_rows]]
  )

  expect_equal(factory_m_metadata(), expected)
})
