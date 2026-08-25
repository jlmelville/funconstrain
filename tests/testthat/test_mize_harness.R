load_mize_harness <- function() {
  path <- system.file(
    "examples",
    "mize-harness.R",
    package = "funconstrain"
  )
  stopifnot(nzchar(path))
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment
}

test_that("mize adapter maps callbacks and counts physical calls", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.2.5.9001")
  harness <- load_mize_harness()
  problem <- funconstrain_problem("rosen")

  callbacks <- harness$as_mize_fg(problem)
  expect_named(callbacks, c("fn", "gr", "fg", "hs"))
  expect_identical(callbacks$hs, problem$he)

  counted <- harness$count_mize_callbacks(callbacks)
  counted$fg$fn(problem$x0)
  counted$fg$gr(problem$x0)
  counted$fg$fg(problem$x0)
  counted$fg$hs(problem$x0)
  expect_identical(
    counted$counts(),
    c(fn = 1L, gr = 1L, fg = 1L, hs = 1L)
  )
})

test_that("ordinary mize runs retain native and independent observations", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.2.5.9001")
  harness <- load_mize_harness()

  run <- harness$run_mize_problem(
    "rosen",
    method = "BFGS",
    controls = list(
      max_iter = 200,
      abs_tol = NULL,
      rel_tol = NULL,
      ginf_tol = 1e-6,
      step_tol = NULL
    )
  )

  expect_identical(run$problem$name, "rosen")
  expect_identical(
    run$problem$configuration$effective,
    list(n = 2L, m = 2L)
  )
  expect_identical(run$method, "BFGS")
  expect_identical(
    run$controls,
    list(
      max_iter = 200,
      abs_tol = NULL,
      rel_tol = NULL,
      ginf_tol = 1e-6,
      step_tol = NULL
    )
  )
  expect_identical(run$reference$status, "applicable")
  expect_null(run$conditions$error)
  expect_length(run$conditions$warnings, 0L)
  expect_identical(run$outcome$status, "converged")
  expect_identical(run$outcome$terminate$what, "ginf_tol")
  expect_type(run$outcome$message, "character")
  expect_type(run$outcome$converged, "logical")
  expect_named(run$outcome$terminate, c("what", "val"))
  expect_true(is.finite(run$quality$initial$objective))
  expect_true(is.finite(run$quality$final$objective))
  expect_true(is.finite(run$quality$final$gradient_inf_norm))
  expect_lt(run$quality$final$objective, run$quality$initial$objective)
  expect_lt(run$quality$final$gradient_inf_norm, 1e-6)
  expect_named(run$work$native, c("nf", "ng", "iter"))
  expect_named(run$work$physical, c("fn", "gr", "fg", "hs"))
  expect_gte(run$work$native$nf, 0)
  expect_gte(run$work$native$ng, 0)
  expect_true(all(run$work$physical >= 0L))
  expect_true(
    numeric_version(run$provenance$mize_version) >=
      numeric_version("0.2.5.9001")
  )
  expect_match(run$provenance$mize_commit, "^[0-9a-f]{40}$")
})

test_that("hard budgets preserve absent native quality without fabrication", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.2.5.9001")
  harness <- load_mize_harness()

  run <- harness$run_mize_problem(
    "rosen",
    method = "BFGS",
    controls = list(
      max_iter = 3,
      max_fg = 1,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )
  )

  expect_identical(run$outcome$status, "budget_exhausted")
  expect_identical(run$outcome$converged, FALSE)
  expect_identical(run$outcome$terminate$what, "max_fg")
  expect_null(run$native_quality$objective)
  expect_true(is.finite(run$quality$final$objective))
  expect_true(is.finite(run$quality$final$gradient_inf_norm))
  expect_lte(run$work$native$nf + run$work$native$ng, 1)
})

test_that("mize harness rejects ambiguous controls and resolves a small suite", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.2.5.9001")
  harness <- load_mize_harness()

  expect_error(
    harness$run_mize_problem("rosen", controls = list(method = "CG")),
    "harness-owned"
  )
  expect_error(
    harness$run_mize_problem(
      "rosen",
      controls = structure(list(1, 2), names = c("max_iter", "max_iter"))
    ),
    "duplicate"
  )

  failed <- harness$run_mize_problem("rosen", method = "not-a-method")
  expect_type(failed$conditions$error, "character")
  expect_gt(nchar(failed$conditions$error), 0L)
  expect_null(failed$outcome$status)
  expect_null(failed$parameters$final)

  manifest <- data.frame(
    name = c("rosen", "ex_rosen", "jenn_samp"),
    n = c(NA_integer_, 4L, NA_integer_),
    m = c(NA_integer_, NA_integer_, 20L)
  )
  suite <- harness$run_mize_suite(
    manifest,
    controls = list(max_iter = 2)
  )

  expect_identical(names(suite$runs), manifest$name)
  expect_identical(
    suite$runs$ex_rosen$problem$configuration$effective,
    list(n = 4L, m = 4L)
  )
  expect_identical(
    suite$runs$ex_rosen$problem$configuration$requested,
    list(n = 4L, m = NULL)
  )
  expect_identical(
    suite$runs$jenn_samp$problem$configuration$effective,
    list(n = 2L, m = 20L)
  )
  expect_identical(suite$runs$jenn_samp$reference$status, "not_applicable")
  expect_identical(suite$provenance, suite$runs$rosen$provenance)
})
