problem_manifest <- function() {
  getFromNamespace("funconstrain_problem_manifest", "funconstrain")()
}

reference_configuration_value <- function(rule, name) {
  value <- rule$source_configuration[[name]]
  if (is.null(value)) {
    return(NULL)
  }
  value
}

reference_signature <- function(rule) {
  configuration <- rule$source_configuration
  configuration_text <- if (length(configuration) == 0L) {
    "*"
  } else {
    paste0(
      names(configuration),
      "=",
      unlist(configuration, use.names = FALSE),
      collapse = ","
    )
  }
  reason_text <- if (is.null(rule$reason)) "" else paste0(":", rule$reason)

  paste0(rule$availability, ":", configuration_text, reason_text)
}

expected_problem_titles <- function() {
  c(
    rosen = "Rosenbrock Function",
    freud_roth = "Freudenstein and Roth Function",
    powell_bs = "Powell Badly Scaled Function",
    brown_bs = "Brown Badly Scaled Function",
    beale = "Beale Function",
    jenn_samp = "Jennrich and Sampson Function",
    helical = "Helical Valley Function",
    bard = "Bard Function",
    gauss = "Gaussian Function",
    meyer = "Meyer Function",
    gulf = "Gulf Research and Development Function",
    box_3d = "Box Three-Dimensional Function",
    powell_s = "Powell Singular Function",
    wood = "Wood Function",
    kow_osb = "Kowalik and Osborne Function",
    brown_den = "Brown and Dennis Function",
    osborne_1 = "Osborne 1 Function",
    biggs_exp6 = "Biggs EXP6 Function",
    osborne_2 = "Osborne 2 Function",
    watson = "Watson Function",
    ex_rosen = "Extended Rosenbrock Function",
    ex_powell = "Extended Powell Function",
    penalty_1 = "Penalty Function I",
    penalty_2 = "Penalty Function II",
    var_dim = "Variably Dimensioned Function",
    trigon = "Trigonometric Function",
    brown_al = "Brown Almost-Linear Function",
    disc_bv = "Discrete Boundary Value Function",
    disc_ie = "Discrete Integral Equation Function",
    broyden_tri = "Broyden Tridiagonal Function",
    broyden_band = "Broyden Banded Function",
    linfun_fr = "Linear Function - Full Rank",
    linfun_r1 = "Linear Function - Rank 1",
    linfun_r1z = "Linear Function - Rank 1 with Zero Columns and Rows",
    chebyquad = "Chebyquad Function"
  )
}

test_that("internal manifest is the authoritative ordered factory registry", {
  manifest <- problem_manifest()
  expected_titles <- expected_problem_titles()
  exported_factories <- setdiff(
    getNamespaceExports("funconstrain"),
    c(
      "fufn",
      "fufnrun",
      "funconstrain_catalog",
      "funconstrain_problem"
    )
  )

  expect_s3_class(manifest, "data.frame")
  expect_equal(nrow(manifest), 35L)
  expect_identical(manifest$number, seq_len(35L))
  expect_equal(length(unique(manifest$name)), 35L)
  expect_identical(manifest$name, names(expected_titles))
  expect_identical(setNames(manifest$title, manifest$name), expected_titles)
  expect_setequal(manifest$name, exported_factories)
  expect_identical(problem_factory_names(), manifest$name)
})

test_that("manifest defaults agree with current factory behavior", {
  manifest <- problem_manifest()

  for (i in seq_len(nrow(manifest))) {
    spec <- manifest[i, ]
    factory <- get_problem_factory(spec$name)
    testfun <- factory()
    info <- spec$name

    expect_identical(
      is.function(testfun$x0),
      spec$n_kind == "variable",
      info = paste(info, "n kind")
    )

    x0 <- standard_x0(testfun)
    expect_equal(length(x0), spec$n_default, info = paste(info, "default n"))

    has_m_argument <- "m" %in% names(formals(factory))
    expect_identical(
      has_m_argument,
      spec$m_kind == "configurable",
      info = paste(info, "m kind")
    )

    if (has_m_argument) {
      if (!is.na(spec$m_n_multiplier)) {
        expect_null(formals(factory)$m, info = paste(info, "derived default m"))
      } else {
        expect_equal(
          formals(factory)$m,
          spec$m_default,
          info = paste(info, "default m")
        )
      }
    }

    if (!is.na(spec$m_n_multiplier)) {
      expect_equal(
        spec$m_default,
        spec$m_n_multiplier * spec$n_default + spec$m_n_offset,
        info = paste(info, "derived default m")
      )
    }
  }
})

test_that("fixed and formula-derived m rules match the MGH definitions", {
  manifest <- problem_manifest()
  expected_fixed <- c(
    rosen = 2L,
    freud_roth = 2L,
    powell_bs = 2L,
    brown_bs = 3L,
    beale = 3L,
    helical = 3L,
    bard = 15L,
    gauss = 15L,
    meyer = 16L,
    powell_s = 4L,
    wood = 6L,
    kow_osb = 11L,
    osborne_1 = 33L,
    osborne_2 = 65L,
    watson = 31L
  )
  expected_derived <- data.frame(
    name = c(
      "ex_rosen",
      "ex_powell",
      "penalty_1",
      "penalty_2",
      "var_dim",
      "trigon",
      "brown_al",
      "disc_bv",
      "disc_ie",
      "broyden_tri",
      "broyden_band",
      "chebyquad"
    ),
    multiplier = c(1L, 1L, 1L, 2L, rep(1L, 8L)),
    offset = c(0L, 0L, 1L, 0L, 2L, rep(0L, 7L)),
    stringsAsFactors = FALSE
  )

  fixed <- manifest[manifest$m_kind == "fixed", ]
  derived <- manifest[!is.na(manifest$m_n_multiplier), ]

  expect_identical(setNames(fixed$m_default, fixed$name), expected_fixed)
  expect_identical(derived$name, expected_derived$name)
  expect_identical(derived$m_n_multiplier, expected_derived$multiplier)
  expect_identical(derived$m_n_offset, expected_derived$offset)
})

test_that("manifest n rules agree with starting-point validation", {
  manifest <- problem_manifest()
  variable <- manifest[manifest$n_kind == "variable", ]

  for (i in seq_len(nrow(variable))) {
    spec <- variable[i, ]
    factory <- get_problem_factory(spec$name)
    testfun <- if (spec$m_kind == "configurable") {
      factory(spec$m_default)
    } else {
      factory()
    }
    info <- spec$name

    expect_equal(length(testfun$x0(spec$n_min)), spec$n_min, info = info)
    expect_error(testfun$x0(spec$n_min - 1L), info = paste(info, "below min"))

    if (is.finite(spec$n_max)) {
      expect_equal(length(testfun$x0(spec$n_max)), spec$n_max, info = info)
      expect_error(
        testfun$x0(spec$n_max + 1L),
        info = paste(info, "above max")
      )
    }

    if (spec$n_multiple > 1L) {
      expect_error(
        testfun$x0(spec$n_min + 1L),
        "multiple",
        info = paste(info, "required multiple")
      )
    }
  }
})

test_that("manifest configurable m rules agree with factory validation", {
  manifest <- problem_manifest()
  configurable <- manifest[manifest$m_kind == "configurable", ]

  for (i in seq_len(nrow(configurable))) {
    spec <- configurable[i, ]
    factory <- get_problem_factory(spec$name)
    info <- spec$name

    expect_silent(factory(spec$m_min))
    expect_error(factory(spec$m_min - 1L), info = paste(info, "below min"))

    if (is.finite(spec$m_max)) {
      expect_silent(factory(spec$m_max))
      expect_error(
        factory(spec$m_max + 1L),
        info = paste(info, "above max")
      )
    }

    if (spec$m_gte_n && spec$n_kind == "variable") {
      testfun <- factory(spec$m_min)
      expect_equal(length(testfun$x0(spec$m_min)), spec$m_min, info = info)
      expect_error(
        testfun$x0(spec$m_min + 1L),
        "m must be >= n",
        info = paste(info, "m >= n")
      )
    }
  }
})

test_that("manifest reference rules have a small coherent vocabulary", {
  manifest <- problem_manifest()
  allowed_availability <- c("stored", "unavailable", "unknown")
  allowed_configuration <- c("n", "n_min", "n_max", "m", "m_min", "m_max")

  for (field in c("fmin_reference", "xmin_reference")) {
    for (i in seq_len(nrow(manifest))) {
      spec <- manifest[i, ]
      rule <- spec[[field]][[1L]]
      info <- paste(spec$name, field)

      expect_named(
        rule,
        c("availability", "source_configuration", "reason"),
        info = info
      )
      expect_true(rule$availability %in% allowed_availability, info = info)
      expect_true(is.list(rule$source_configuration), info = info)
      expect_true(
        all(names(rule$source_configuration) %in% allowed_configuration),
        info = info
      )

      if (rule$availability == "stored") {
        expect_null(rule$reason, info = info)
      } else {
        expect_true(is.character(rule$reason), info = info)
        expect_equal(length(rule$reason), 1L, info = info)
      }
    }
  }

  fmin_availability <- vapply(
    manifest$fmin_reference,
    `[[`,
    character(1L),
    "availability"
  )
  xmin_availability <- vapply(
    manifest$xmin_reference,
    `[[`,
    character(1L),
    "availability"
  )

  expect_true(all(fmin_availability == "stored"))
  expect_identical(
    manifest$name[xmin_availability == "unavailable"],
    "box_3d"
  )
  expect_false(any(xmin_availability == "unknown"))
})

test_that("reference applicability is keyed to the documented problems", {
  manifest <- problem_manifest()
  problem_names <- names(expected_problem_titles())
  expected_fmin <- setNames(
    rep("stored:*", length(problem_names)),
    problem_names
  )
  expected_xmin <- expected_fmin

  expected_fmin[c(
    "jenn_samp",
    "brown_den",
    "watson",
    "penalty_1",
    "penalty_2",
    "linfun_fr",
    "linfun_r1",
    "linfun_r1z",
    "chebyquad"
  )] <- c(
    "stored:m=10",
    "stored:m=20",
    "stored:n=6",
    "stored:n=4",
    "stored:n=4",
    "stored:n=4,m=100",
    "stored:m=100",
    "stored:n_min=3,m=100",
    "stored:n=8,m=8"
  )
  expected_xmin[c(
    "jenn_samp",
    "box_3d",
    "brown_den",
    "watson",
    "ex_rosen",
    "ex_powell",
    "penalty_1",
    "penalty_2",
    "var_dim",
    "trigon",
    "brown_al",
    "disc_bv",
    "disc_ie",
    "broyden_tri",
    "broyden_band",
    "linfun_fr",
    "linfun_r1",
    "linfun_r1z",
    "chebyquad"
  )] <- c(
    "stored:m=10",
    "unavailable:*:no_single_minimizer",
    "stored:m=20",
    "stored:n=6",
    "stored:n=8",
    "stored:n=4",
    "stored:n=4",
    "stored:n=4",
    "stored:n=6",
    "stored:n=4",
    "stored:n=4",
    "stored:n=5",
    "stored:n=5",
    "stored:n=5",
    "stored:n=5",
    "stored:n=4,m=100",
    "stored:n=5,m=100",
    "stored:n=5,m=100",
    "stored:n=8,m=8"
  )

  actual_fmin <- setNames(
    vapply(manifest$fmin_reference, reference_signature, character(1L)),
    manifest$name
  )
  actual_xmin <- setNames(
    vapply(manifest$xmin_reference, reference_signature, character(1L)),
    manifest$name
  )

  expect_identical(actual_fmin, expected_fmin)
  expect_identical(actual_xmin, expected_xmin)
})

test_that("stored reference configurations agree with factory literals", {
  manifest <- problem_manifest()

  for (i in seq_len(nrow(manifest))) {
    spec <- manifest[i, ]
    xmin_rule <- spec$xmin_reference[[1L]]
    reference_m <- reference_configuration_value(xmin_rule, "m")
    factory <- get_problem_factory(spec$name)
    testfun <- if (is.null(reference_m)) factory() else factory(reference_m)
    info <- spec$name

    if (xmin_rule$availability == "unavailable") {
      expect_true(anyNA(testfun$xmin), info = info)
      next
    }

    expect_false(anyNA(testfun$xmin), info = info)
    reference_n <- reference_configuration_value(xmin_rule, "n")
    if (!is.null(reference_n)) {
      expect_equal(length(testfun$xmin), reference_n, info = info)
    } else if (spec$n_kind == "fixed") {
      expect_equal(length(testfun$xmin), spec$n_default, info = info)
    }
  }
})
