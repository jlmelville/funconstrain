test_that("catalog is a minimal ordered public projection", {
  catalog <- funconstrain_catalog()
  manifest <- getFromNamespace(
    "funconstrain_problem_manifest",
    "funconstrain"
  )()

  expected <- data.frame(
    number = seq_len(35L),
    name = manifest$name,
    title = manifest$title,
    n_kind = manifest$n_kind,
    n_default = manifest$n_default,
    m_configurable = manifest$m_kind == "configurable",
    m_default = manifest$m_default,
    stringsAsFactors = FALSE
  )
  names(expected) <- c(
    "number",
    "name",
    "title",
    "n_kind",
    "n_default",
    "m_configurable",
    "m_default"
  )

  expect_identical(catalog, expected)
})

test_that("every catalog default resolves to the public list contract", {
  catalog <- funconstrain_catalog()
  resolved <- expect_no_warning(lapply(seq_len(nrow(catalog)), function(i) {
    entry <- catalog[i, ]
    funconstrain_problem(entry$name)
  }))
  names(resolved) <- catalog$name
  valid_reference_statuses <- c(
    "applicable",
    "not_applicable",
    "unavailable",
    "unknown"
  )
  checks <- lapply(seq_len(nrow(catalog)), function(i) {
    entry <- catalog[i, ]
    case <- resolved[[i]]
    reference_checks <- lapply(case$reference, function(reference) {
      list(
        fields_match = identical(
          names(reference),
          c("value", "status", "reason", "source_configuration")
        ),
        status_is_valid = reference$status %in% valid_reference_statuses,
        reason_is_scalar_character = is.character(reference$reason) &&
          length(reference$reason) == 1L,
        source_configuration_is_list = is.list(reference$source_configuration)
      )
    })

    list(
      is_unclassed = is.null(attr(case, "class")),
      fields_match = identical(
        names(case),
        c(
          "fn",
          "gr",
          "he",
          "fg",
          "x0",
          "problem",
          "configuration",
          "reference"
        )
      ),
      problem_fields_match = identical(
        names(case$problem),
        c("number", "name", "title")
      ),
      problem_matches = identical(
        case$problem,
        list(
          number = entry$number,
          name = entry$name,
          title = entry$title
        )
      ),
      configuration_fields_match = identical(
        names(case$configuration),
        c("requested", "effective")
      ),
      requested_configuration_matches = identical(
        case$configuration$requested,
        list(n = NULL, m = NULL)
      ),
      effective_configuration_matches = identical(
        case$configuration$effective,
        list(n = entry$n_default, m = entry$m_default)
      ),
      x0_is_numeric = is.numeric(case$x0),
      x0_has_no_missing_values = !anyNA(case$x0),
      x0_is_finite = all(is.finite(case$x0)),
      x0_length_matches = length(case$x0) == entry$n_default,
      reference_fields_match = identical(
        names(case$reference),
        c("fmin", "xmin")
      ),
      references = reference_checks
    )
  })
  names(checks) <- catalog$name

  expect_all_checks(checks)
})

test_that("resolved callbacks are pinned and satisfy the core contract", {
  factory_names <- funconstrain_catalog()$name
  checks <- lapply(factory_names, function(name) {
    case <- funconstrain_problem(name)
    n <- case$configuration$effective$n
    value <- case$fn(case$x0)
    gradient <- case$gr(case$x0)
    hessian <- case$he(case$x0)
    combined <- case$fg(case$x0)

    list(
      value_is_numeric = is.numeric(value),
      value_is_scalar = length(value) == 1L,
      gradient_length_matches = length(gradient) == n,
      hessian_is_matrix = is.matrix(hessian),
      hessian_dimensions_match = identical(dim(hessian), c(n, n)),
      hessian_is_symmetric = isSymmetric(hessian),
      combined_fn_matches = isTRUE(all.equal(combined$fn, value)),
      combined_gr_matches = isTRUE(all.equal(combined$gr, gradient))
    )
  })
  names(checks) <- factory_names

  expect_all_checks(checks)
})

test_that("representative resolved callbacks enforce their pinned dimensions", {
  cases <- list(
    fixed = funconstrain_problem("rosen"),
    variable = funconstrain_problem("ex_rosen", n = 8),
    configurable = funconstrain_problem("jenn_samp", m = 20)
  )
  actual <- character()
  expected <- character()

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    n <- case$configuration$effective$n
    for (callback in c("fn", "gr", "he", "fg")) {
      key <- paste(case_name, callback)
      actual[[key]] <- capture_error_message(
        case[[callback]](rep(0, n + 1L))
      )
      expected[[key]] <- paste0(
        case$problem$title,
        ": resolved callbacks require a parameter vector of length ",
        n
      )
    }
  }

  expect_identical(actual, expected)
})

test_that("explicit valid dimensions are normalized and recorded", {
  fixed <- funconstrain_problem("rosen", n = 2, m = 2)
  variable <- funconstrain_problem("ex_rosen", n = 10)
  configurable <- funconstrain_problem("jenn_samp", m = 20)
  chebyquad <- funconstrain_problem("chebyquad", n = 4, m = 7)
  direct_chebyquad <- getExportedValue("funconstrain", "chebyquad")(m = 7)
  par <- c(0.12, 0.37, 0.64, 0.91)
  derived <- funconstrain_problem("penalty_1", n = 4, m = 5)

  actual <- list(
    fixed = fixed$configuration,
    variable = variable$configuration,
    configurable = configurable$configuration,
    chebyquad = chebyquad$configuration,
    chebyquad_objective = chebyquad$fn(par),
    derived_effective = derived$configuration$effective
  )
  expected <- list(
    fixed = list(
      requested = list(n = 2L, m = 2L),
      effective = list(n = 2L, m = 2L)
    ),
    variable = list(
      requested = list(n = 10L, m = NULL),
      effective = list(n = 10L, m = 10L)
    ),
    configurable = list(
      requested = list(n = NULL, m = 20L),
      effective = list(n = 2L, m = 20L)
    ),
    chebyquad = list(
      requested = list(n = 4L, m = 7L),
      effective = list(n = 4L, m = 7L)
    ),
    chebyquad_objective = direct_chebyquad$fn(par),
    derived_effective = list(n = 4L, m = 5L)
  )

  expect_equal(actual, expected)
})

test_that("nondefault derived m values follow independent formulas", {
  cases <- data.frame(
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
    n = c(6L, 8L, rep(6L, 10L)),
    m = c(6L, 8L, 7L, 12L, 8L, rep(6L, 7L)),
    stringsAsFactors = FALSE
  )

  actual <- vapply(
    seq_len(nrow(cases)),
    function(i) {
      case <- cases[i, ]
      funconstrain_problem(case$name, n = case$n)$configuration$effective$m
    },
    integer(1L)
  )

  expect_identical(setNames(actual, cases$name), setNames(cases$m, cases$name))
})

test_that("invalid explicit requests error without falling back", {
  expect_error(funconstrain_problem(1), "name must be")
  expect_error(funconstrain_problem(NA_character_), "name must be")
  expect_error(funconstrain_problem("rose"), "exact name")
  expect_error(funconstrain_problem("Rosenbrock Function"), "exact name")

  expect_error(funconstrain_problem("rosen", n = 3), "n is outside")
  expect_error(funconstrain_problem("rosen", m = 3), "m is outside")
  expect_error(funconstrain_problem("watson", n = 1), "n is outside")
  expect_error(funconstrain_problem("watson", n = 32), "n is outside")
  expect_error(funconstrain_problem("ex_rosen", n = 3), "multiple of 2")
  expect_error(funconstrain_problem("ex_powell", n = 6), "multiple of 4")
  expect_error(funconstrain_problem("gulf", m = 2), "m is outside")
  expect_error(funconstrain_problem("gulf", m = 101), "m is outside")
  expect_error(funconstrain_problem("penalty_1", n = 4, m = 6), "m is outside")
  expect_error(
    funconstrain_problem("linfun_fr", n = 11, m = 10),
    "m must be greater than or equal to n"
  )
  expect_error(
    funconstrain_problem("chebyquad", n = 4, m = 3),
    "m must be greater than or equal to n"
  )

  for (bad_n in list(TRUE, NA_real_, Inf, 2.5, c(2, 4))) {
    expect_error(
      funconstrain_problem("ex_rosen", n = bad_n),
      "finite whole-number scalar"
    )
  }
})

test_that("resolver enforces every encoded dimension boundary", {
  manifest <- getFromNamespace(
    "funconstrain_problem_manifest",
    "funconstrain"
  )()
  checks <- logical()

  fixed_n <- manifest[manifest$n_kind == "fixed", ]
  for (i in seq_len(nrow(fixed_n))) {
    spec <- fixed_n[i, ]
    checks[[paste(spec$name, "fixed n")]] <- is.na(
      capture_error_message(funconstrain_problem(spec$name, n = spec$n_default))
    )
    checks[[paste(spec$name, "below fixed n")]] <- !is.na(
      capture_error_message(
        funconstrain_problem(spec$name, n = spec$n_default - 1L)
      )
    )
    checks[[paste(spec$name, "above fixed n")]] <- !is.na(
      capture_error_message(
        funconstrain_problem(spec$name, n = spec$n_default + 1L)
      )
    )
  }

  variable_n <- manifest[manifest$n_kind == "variable", ]
  for (i in seq_len(nrow(variable_n))) {
    spec <- variable_n[i, ]
    valid <- funconstrain_problem(spec$name, n = spec$n_min)
    checks[[paste(spec$name, "minimum n")]] <- identical(
      valid$configuration$effective$n,
      spec$n_min
    )
    checks[[paste(spec$name, "below minimum n")]] <- !is.na(
      capture_error_message(funconstrain_problem(
        spec$name,
        n = spec$n_min - 1L
      ))
    )

    if (is.finite(spec$n_max)) {
      checks[[paste(spec$name, "maximum n")]] <- is.na(
        capture_error_message(funconstrain_problem(spec$name, n = spec$n_max))
      )
      checks[[paste(spec$name, "above maximum n")]] <- !is.na(
        capture_error_message(
          funconstrain_problem(spec$name, n = spec$n_max + 1L)
        )
      )
    }
    if (spec$n_multiple > 1L) {
      checks[[paste(spec$name, "n multiple")]] <- grepl(
        "multiple",
        capture_error_message(
          funconstrain_problem(spec$name, n = spec$n_min + 1L)
        ),
        fixed = TRUE
      )
    }
  }

  configurable_m <- manifest[manifest$m_kind == "configurable", ]
  for (i in seq_len(nrow(configurable_m))) {
    spec <- configurable_m[i, ]
    request_n <- if (spec$n_kind == "variable") spec$n_min else NULL
    valid <- funconstrain_problem(spec$name, n = request_n, m = spec$m_min)
    checks[[paste(spec$name, "minimum m")]] <- identical(
      valid$configuration$effective$m,
      spec$m_min
    )
    checks[[paste(spec$name, "below minimum m")]] <- !is.na(
      capture_error_message(
        funconstrain_problem(spec$name, n = request_n, m = spec$m_min - 1L)
      )
    )

    if (is.finite(spec$m_max)) {
      checks[[paste(spec$name, "maximum m")]] <- is.na(
        capture_error_message(
          funconstrain_problem(spec$name, n = request_n, m = spec$m_max)
        )
      )
      checks[[paste(spec$name, "above maximum m")]] <- !is.na(
        capture_error_message(
          funconstrain_problem(spec$name, n = request_n, m = spec$m_max + 1L)
        )
      )
    }
    if (spec$m_gte_n && spec$n_kind == "variable") {
      checks[[paste(spec$name, "m below n")]] <- grepl(
        "greater than or equal to n",
        capture_error_message(
          funconstrain_problem(
            spec$name,
            n = spec$m_min + 1L,
            m = spec$m_min
          )
        ),
        fixed = TRUE
      )
    }
  }

  derived_m <- manifest[manifest$m_kind == "derived", ]
  for (i in seq_len(nrow(derived_m))) {
    spec <- derived_m[i, ]
    n <- spec$n_min
    expected_m <- spec$m_n_multiplier * n + spec$m_n_offset
    valid <- funconstrain_problem(spec$name, n = n, m = expected_m)
    checks[[paste(spec$name, "derived m")]] <- identical(
      valid$configuration$effective$m,
      as.integer(expected_m)
    )
    checks[[paste(spec$name, "invalid derived m")]] <- !is.na(
      capture_error_message(
        funconstrain_problem(spec$name, n = n, m = expected_m + 1L)
      )
    )
  }

  expect_all_checks(checks)
})

test_that("reference records distinguish their public states", {
  applicable <- funconstrain_problem("rosen")$reference
  split <- funconstrain_problem("ex_rosen", n = 4)$reference
  mismatch <- funconstrain_problem("jenn_samp", m = 20)$reference
  chebyquad_square <- funconstrain_problem("chebyquad", n = 8, m = 8)
  chebyquad_overdetermined <- funconstrain_problem("chebyquad", n = 8, m = 9)
  bounded <- funconstrain_problem("linfun_r1z", n = 2, m = 100)$reference
  unavailable <- funconstrain_problem("box_3d")$reference$xmin

  actual <- list(
    applicable_fmin = applicable$fmin,
    applicable_xmin_status = applicable$xmin$status,
    split = list(
      fmin_status = split$fmin$status,
      xmin_status = split$xmin$status,
      xmin_reason = split$xmin$reason,
      xmin_source_configuration = split$xmin$source_configuration,
      xmin_length = length(split$xmin$value)
    ),
    mismatch = list(
      fmin_status = mismatch$fmin$status,
      xmin_status = mismatch$xmin$status,
      fmin_source_configuration = mismatch$fmin$source_configuration
    ),
    chebyquad_square = list(
      fmin_status = chebyquad_square$reference$fmin$status,
      xmin_status = chebyquad_square$reference$xmin$status,
      source_configuration = chebyquad_square$reference$fmin$source_configuration
    ),
    chebyquad_overdetermined = list(
      fmin_status = chebyquad_overdetermined$reference$fmin$status,
      xmin_status = chebyquad_overdetermined$reference$xmin$status
    ),
    bounded = list(
      fmin_status = bounded$fmin$status,
      source_configuration = bounded$fmin$source_configuration
    ),
    unavailable = unavailable
  )
  expected <- list(
    applicable_fmin = list(
      value = 0,
      status = "applicable",
      reason = "documented_applicability",
      source_configuration = list()
    ),
    applicable_xmin_status = "applicable",
    split = list(
      fmin_status = "applicable",
      xmin_status = "not_applicable",
      xmin_reason = "source_configuration_mismatch",
      xmin_source_configuration = list(n = 8L),
      xmin_length = 8L
    ),
    mismatch = list(
      fmin_status = "not_applicable",
      xmin_status = "not_applicable",
      fmin_source_configuration = list(m = 10L)
    ),
    chebyquad_square = list(
      fmin_status = "applicable",
      xmin_status = "applicable",
      source_configuration = list(n = 8L, m = 8L)
    ),
    chebyquad_overdetermined = list(
      fmin_status = "not_applicable",
      xmin_status = "not_applicable"
    ),
    bounded = list(
      fmin_status = "not_applicable",
      source_configuration = list(n_min = 3L, m = 100L)
    ),
    unavailable = list(
      value = NULL,
      status = "unavailable",
      reason = "no_single_minimizer",
      source_configuration = list()
    )
  )

  expect_identical(actual, expected)
})

test_that("unknown reference applicability has an explicit record", {
  resolve_reference <- getFromNamespace(
    "resolve_reference_record",
    "funconstrain"
  )
  rule <- list(
    availability = "unknown",
    source_configuration = list(),
    reason = "unencoded_rule"
  )

  expect_identical(
    resolve_reference(7, rule, n = 2L, m = 2L),
    list(
      value = 7,
      status = "unknown",
      reason = "unencoded_rule",
      source_configuration = list()
    )
  )
})

test_that("jointly applicable resolved minima agree", {
  checks <- logical()
  for (name in funconstrain_catalog()$name) {
    case <- funconstrain_problem(name)
    fmin <- case$reference$fmin
    xmin <- case$reference$xmin
    if (fmin$status != "applicable" || xmin$status != "applicable") {
      next
    }

    actual <- case$fn(xmin$value)
    tolerance <- max(1e-8, abs(fmin$value) * 1e-5)
    checks[[paste(name, "finite")]] <- is.finite(actual)
    checks[[paste(name, "matches fmin")]] <-
      abs(actual - fmin$value) <= tolerance
  }

  expect_all_checks(checks)
})
