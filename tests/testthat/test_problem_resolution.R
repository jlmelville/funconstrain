test_that("catalog is a minimal ordered public projection", {
  catalog <- funconstrain_catalog()
  manifest <- getFromNamespace(
    "funconstrain_problem_manifest",
    "funconstrain"
  )()

  expect_s3_class(catalog, "data.frame")
  expect_named(
    catalog,
    c(
      "number",
      "name",
      "title",
      "n_kind",
      "n_default",
      "m_configurable",
      "m_default"
    )
  )
  expect_equal(nrow(catalog), 35L)
  expect_identical(catalog$number, seq_len(35L))
  expect_identical(catalog$name, manifest$name)
  expect_identical(catalog$title, manifest$title)
  expect_identical(catalog$n_kind, manifest$n_kind)
  expect_identical(catalog$n_default, manifest$n_default)
  expect_identical(
    catalog$m_configurable,
    manifest$m_kind == "configurable"
  )
  expect_identical(catalog$m_default, manifest$m_default)
})

test_that("every catalog default resolves to the public list contract", {
  catalog <- funconstrain_catalog()

  for (i in seq_len(nrow(catalog))) {
    entry <- catalog[i, ]
    info <- entry$name
    expect_warning(case <- funconstrain_problem(entry$name), NA, info = info)

    expect_null(attr(case, "class"), info = info)
    expect_named(
      case,
      c("fn", "gr", "he", "fg", "x0", "problem", "configuration", "reference"),
      info = info
    )
    expect_named(case$problem, c("number", "name", "title"), info = info)
    expect_identical(
      case$problem,
      list(
        number = entry$number,
        name = entry$name,
        title = entry$title
      ),
      info = info
    )
    expect_named(case$configuration, c("requested", "effective"), info = info)
    expect_identical(
      case$configuration$requested,
      list(n = NULL, m = NULL),
      info = info
    )
    expect_identical(
      case$configuration$effective,
      list(n = entry$n_default, m = entry$m_default),
      info = info
    )
    expect_true(is.numeric(case$x0), info = info)
    expect_false(anyNA(case$x0), info = info)
    expect_true(all(is.finite(case$x0)), info = info)
    expect_equal(length(case$x0), entry$n_default, info = info)
    expect_named(case$reference, c("fmin", "xmin"), info = info)
    for (reference in case$reference) {
      expect_named(
        reference,
        c("value", "status", "reason", "source_configuration"),
        info = info
      )
      expect_true(
        reference$status %in%
          c(
            "applicable",
            "not_applicable",
            "unavailable",
            "unknown"
          ),
        info = info
      )
      expect_true(is.character(reference$reason), info = info)
      expect_length(reference$reason, 1L)
      expect_true(is.list(reference$source_configuration), info = info)
    }
  }
})

test_that("resolved callbacks are pinned and satisfy the core contract", {
  for (name in funconstrain_catalog()$name) {
    case <- funconstrain_problem(name)
    n <- case$configuration$effective$n
    info <- name

    value <- case$fn(case$x0)
    gradient <- case$gr(case$x0)
    hessian <- case$he(case$x0)
    combined <- case$fg(case$x0)

    expect_true(is.numeric(value), info = paste(info, "fn type"))
    expect_equal(length(value), 1L, info = paste(info, "fn length"))
    expect_equal(length(gradient), n, info = paste(info, "gr length"))
    expect_true(is.matrix(hessian), info = paste(info, "he matrix"))
    expect_equal(dim(hessian), c(n, n), info = paste(info, "he dimensions"))
    expect_true(isSymmetric(hessian), info = paste(info, "he symmetry"))
    expect_equal(combined$fn, value, info = paste(info, "fg fn"))
    expect_equal(combined$gr, gradient, info = paste(info, "fg gr"))

    for (callback in c("fn", "gr", "he", "fg")) {
      expect_error(
        case[[callback]](rep(0, n + 1L)),
        paste0("parameter vector of length ", n, "$"),
        info = paste(info, callback, "pinned dimension")
      )
    }
  }
})

test_that("explicit valid dimensions are normalized and recorded", {
  fixed <- funconstrain_problem("rosen", n = 2, m = 2)
  expect_identical(fixed$configuration$requested, list(n = 2L, m = 2L))
  expect_identical(fixed$configuration$effective, list(n = 2L, m = 2L))

  variable <- funconstrain_problem("ex_rosen", n = 10)
  expect_identical(variable$configuration$requested, list(n = 10L, m = NULL))
  expect_identical(variable$configuration$effective, list(n = 10L, m = 10L))

  configurable <- funconstrain_problem("jenn_samp", m = 20)
  expect_identical(
    configurable$configuration$requested,
    list(n = NULL, m = 20L)
  )
  expect_identical(
    configurable$configuration$effective,
    list(n = 2L, m = 20L)
  )

  derived <- funconstrain_problem("penalty_1", n = 4, m = 5)
  expect_identical(derived$configuration$effective, list(n = 4L, m = 5L))
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

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    resolved <- funconstrain_problem(case$name, n = case$n)
    expect_identical(
      resolved$configuration$effective$m,
      case$m,
      info = case$name
    )
  }
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

  fixed_n <- manifest[manifest$n_kind == "fixed", ]
  for (i in seq_len(nrow(fixed_n))) {
    spec <- fixed_n[i, ]
    expect_no_error(funconstrain_problem(spec$name, n = spec$n_default))
    expect_error(funconstrain_problem(spec$name, n = spec$n_default - 1L))
    expect_error(funconstrain_problem(spec$name, n = spec$n_default + 1L))
  }

  variable_n <- manifest[manifest$n_kind == "variable", ]
  for (i in seq_len(nrow(variable_n))) {
    spec <- variable_n[i, ]
    valid <- funconstrain_problem(spec$name, n = spec$n_min)
    expect_identical(valid$configuration$effective$n, spec$n_min)
    expect_error(funconstrain_problem(spec$name, n = spec$n_min - 1L))

    if (is.finite(spec$n_max)) {
      expect_no_error(funconstrain_problem(spec$name, n = spec$n_max))
      expect_error(funconstrain_problem(spec$name, n = spec$n_max + 1L))
    }
    if (spec$n_multiple > 1L) {
      expect_error(
        funconstrain_problem(spec$name, n = spec$n_min + 1L),
        "multiple"
      )
    }
  }

  configurable_m <- manifest[manifest$m_kind == "configurable", ]
  for (i in seq_len(nrow(configurable_m))) {
    spec <- configurable_m[i, ]
    request_n <- if (spec$n_kind == "variable") spec$n_min else NULL
    valid <- funconstrain_problem(spec$name, n = request_n, m = spec$m_min)
    expect_identical(valid$configuration$effective$m, spec$m_min)
    expect_error(
      funconstrain_problem(spec$name, n = request_n, m = spec$m_min - 1L)
    )

    if (is.finite(spec$m_max)) {
      expect_no_error(
        funconstrain_problem(spec$name, n = request_n, m = spec$m_max)
      )
      expect_error(
        funconstrain_problem(spec$name, n = request_n, m = spec$m_max + 1L)
      )
    }
    if (spec$m_gte_n && spec$n_kind == "variable") {
      expect_error(
        funconstrain_problem(
          spec$name,
          n = spec$m_min + 1L,
          m = spec$m_min
        ),
        "greater than or equal to n"
      )
    }
  }

  derived_m <- manifest[manifest$m_kind == "derived", ]
  for (i in seq_len(nrow(derived_m))) {
    spec <- derived_m[i, ]
    n <- spec$n_min
    expected_m <- spec$m_n_multiplier * n + spec$m_n_offset
    valid <- funconstrain_problem(spec$name, n = n, m = expected_m)
    expect_identical(valid$configuration$effective$m, as.integer(expected_m))
    expect_error(funconstrain_problem(spec$name, n = n, m = expected_m + 1L))
  }
})

test_that("reference records distinguish their public states", {
  applicable <- funconstrain_problem("rosen")$reference
  expect_identical(
    applicable$fmin,
    list(
      value = 0,
      status = "applicable",
      reason = "documented_applicability",
      source_configuration = list()
    )
  )
  expect_identical(applicable$xmin$status, "applicable")

  split <- funconstrain_problem("ex_rosen", n = 4)$reference
  expect_identical(split$fmin$status, "applicable")
  expect_identical(split$xmin$status, "not_applicable")
  expect_identical(split$xmin$reason, "source_configuration_mismatch")
  expect_identical(split$xmin$source_configuration, list(n = 8L))
  expect_equal(length(split$xmin$value), 8L)

  mismatch <- funconstrain_problem("jenn_samp", m = 20)$reference
  expect_identical(mismatch$fmin$status, "not_applicable")
  expect_identical(mismatch$xmin$status, "not_applicable")
  expect_identical(mismatch$fmin$source_configuration, list(m = 10L))

  bounded <- funconstrain_problem("linfun_r1z", n = 2, m = 100)$reference
  expect_identical(bounded$fmin$status, "not_applicable")
  expect_identical(
    bounded$fmin$source_configuration,
    list(n_min = 3L, m = 100L)
  )

  unavailable <- funconstrain_problem("box_3d")$reference$xmin
  expect_null(unavailable$value)
  expect_identical(unavailable$status, "unavailable")
  expect_identical(unavailable$reason, "no_single_minimizer")
  expect_identical(unavailable$source_configuration, list())
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
  for (name in funconstrain_catalog()$name) {
    case <- funconstrain_problem(name)
    fmin <- case$reference$fmin
    xmin <- case$reference$xmin
    if (fmin$status != "applicable" || xmin$status != "applicable") {
      next
    }

    actual <- case$fn(xmin$value)
    tolerance <- max(1e-8, abs(fmin$value) * 1e-5)
    expect_true(is.finite(actual), info = name)
    expect_true(abs(actual - fmin$value) <= tolerance, info = name)
  }
})
