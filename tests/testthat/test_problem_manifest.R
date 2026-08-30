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

  actual <- list(
    is_data_frame = is.data.frame(manifest),
    row_count = nrow(manifest),
    number = manifest$number,
    unique_name_count = length(unique(manifest$name)),
    name = manifest$name,
    title = setNames(manifest$title, manifest$name),
    exported_factories = sort(manifest$name),
    helper_names = problem_factory_names()
  )
  expected <- list(
    is_data_frame = TRUE,
    row_count = 35L,
    number = seq_len(35L),
    unique_name_count = 35L,
    name = names(expected_titles),
    title = expected_titles,
    exported_factories = sort(exported_factories),
    helper_names = manifest$name
  )

  expect_identical(actual, expected)
})

test_that("manifest defaults agree with current factory behavior", {
  manifest <- problem_manifest()
  checks <- lapply(seq_len(nrow(manifest)), function(i) {
    spec <- manifest[i, ]
    factory <- get_problem_factory(spec$name)
    testfun <- factory()
    x0 <- standard_x0(testfun)
    has_m_argument <- "m" %in% names(formals(factory))
    m_formal_matches <- if (has_m_argument) {
      if (!is.na(spec$m_n_multiplier)) {
        is.null(formals(factory)$m)
      } else {
        isTRUE(all.equal(formals(factory)$m, spec$m_default))
      }
    } else {
      TRUE
    }
    derived_default_matches <- if (!is.na(spec$m_n_multiplier)) {
      isTRUE(all.equal(
        spec$m_default,
        spec$m_n_multiplier * spec$n_default + spec$m_n_offset
      ))
    } else {
      TRUE
    }

    list(
      n_kind_matches = identical(
        is.function(testfun$x0),
        spec$n_kind == "variable"
      ),
      default_n_matches = identical(length(x0), spec$n_default),
      m_kind_matches = identical(
        has_m_argument,
        spec$m_kind == "configurable"
      ),
      m_formal_matches = m_formal_matches,
      derived_default_matches = derived_default_matches
    )
  })
  names(checks) <- manifest$name

  expect_all_checks(checks)
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
  checks <- logical()

  for (i in seq_len(nrow(variable))) {
    spec <- variable[i, ]
    factory <- get_problem_factory(spec$name)
    testfun <- if (spec$m_kind == "configurable") {
      factory(spec$m_default)
    } else {
      factory()
    }
    checks[[paste(spec$name, "minimum")]] <-
      length(testfun$x0(spec$n_min)) == spec$n_min
    checks[[paste(spec$name, "below minimum")]] <- !is.na(
      capture_error_message(testfun$x0(spec$n_min - 1L))
    )

    if (is.finite(spec$n_max)) {
      checks[[paste(spec$name, "maximum")]] <-
        length(testfun$x0(spec$n_max)) == spec$n_max
      checks[[paste(spec$name, "above maximum")]] <- !is.na(
        capture_error_message(testfun$x0(spec$n_max + 1L))
      )
    }

    if (spec$n_multiple > 1L) {
      checks[[paste(spec$name, "required multiple")]] <- grepl(
        "multiple",
        capture_error_message(testfun$x0(spec$n_min + 1L)),
        fixed = TRUE
      )
    }
  }

  expect_all_checks(checks)
})

test_that("manifest configurable m rules agree with factory validation", {
  manifest <- problem_manifest()
  configurable <- manifest[manifest$m_kind == "configurable", ]
  checks <- logical()

  for (i in seq_len(nrow(configurable))) {
    spec <- configurable[i, ]
    factory <- get_problem_factory(spec$name)
    checks[[paste(spec$name, "minimum")]] <- is.na(
      capture_error_message(factory(spec$m_min))
    )
    checks[[paste(spec$name, "below minimum")]] <- !is.na(
      capture_error_message(factory(spec$m_min - 1L))
    )

    if (is.finite(spec$m_max)) {
      checks[[paste(spec$name, "maximum")]] <- is.na(
        capture_error_message(factory(spec$m_max))
      )
      checks[[paste(spec$name, "above maximum")]] <- !is.na(
        capture_error_message(factory(spec$m_max + 1L))
      )
    }

    if (spec$m_gte_n && spec$n_kind == "variable") {
      testfun <- factory(spec$m_min)
      checks[[paste(spec$name, "m equals n")]] <-
        length(testfun$x0(spec$m_min)) == spec$m_min
      checks[[paste(spec$name, "m below n")]] <- grepl(
        "m must be >= n",
        capture_error_message(testfun$x0(spec$m_min + 1L)),
        fixed = TRUE
      )
    }
  }

  expect_all_checks(checks)
})

test_that("manifest reference rules have a small coherent vocabulary", {
  manifest <- problem_manifest()
  allowed_availability <- c("stored", "unavailable", "unknown")
  allowed_configuration <- c("n", "n_min", "n_max", "m", "m_min", "m_max")
  checks <- logical()

  for (field in c("fmin_reference", "xmin_reference")) {
    for (i in seq_len(nrow(manifest))) {
      spec <- manifest[i, ]
      rule <- spec[[field]][[1L]]
      key <- paste(spec$name, field)

      checks[[paste(key, "fields")]] <- identical(
        names(rule),
        c("availability", "source_configuration", "reason")
      )
      checks[[paste(key, "availability")]] <-
        rule$availability %in% allowed_availability
      checks[[paste(key, "configuration type")]] <-
        is.list(rule$source_configuration)
      checks[[paste(key, "configuration fields")]] <-
        all(names(rule$source_configuration) %in% allowed_configuration)

      if (rule$availability == "stored") {
        checks[[paste(key, "reason")]] <- is.null(rule$reason)
      } else {
        checks[[paste(key, "reason type")]] <- is.character(rule$reason)
        checks[[paste(key, "reason length")]] <- length(rule$reason) == 1L
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

  checks[["all fmin references stored"]] <- all(fmin_availability == "stored")
  checks[["box_3d xmin unavailable"]] <- identical(
    manifest$name[xmin_availability == "unavailable"],
    "box_3d"
  )
  checks[["no unknown xmin references"]] <- !any(xmin_availability == "unknown")

  expect_all_checks(checks)
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
  checks <- logical()

  for (i in seq_len(nrow(manifest))) {
    spec <- manifest[i, ]
    xmin_rule <- spec$xmin_reference[[1L]]
    reference_m <- reference_configuration_value(xmin_rule, "m")
    factory <- get_problem_factory(spec$name)
    testfun <- if (is.null(reference_m)) factory() else factory(reference_m)

    if (xmin_rule$availability == "unavailable") {
      checks[[paste(spec$name, "xmin unavailable")]] <- anyNA(testfun$xmin)
      next
    }

    checks[[paste(spec$name, "xmin available")]] <- !anyNA(testfun$xmin)
    reference_n <- reference_configuration_value(xmin_rule, "n")
    if (!is.null(reference_n)) {
      checks[[paste(spec$name, "reference n")]] <-
        length(testfun$xmin) == reference_n
    } else if (spec$n_kind == "fixed") {
      checks[[paste(spec$name, "fixed n")]] <-
        length(testfun$xmin) == spec$n_default
    }
  }

  expect_all_checks(checks)
})
