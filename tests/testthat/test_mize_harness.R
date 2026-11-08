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
  harness <- load_mize_harness()
  problem <- funconstrain_problem("rosen")

  callbacks <- harness$as_mize_fg(problem)
  expect_named(callbacks, c("fn", "gr", "fg", "hs"))
  expect_identical(callbacks$hs, problem$he)

  counted <- harness$count_mize_callbacks(callbacks)
  counted$fg$fn(problem$x0)
  counted$fg$fn(problem$x0)
  counted$fg$gr(problem$x0)
  counted$fg$fg(problem$x0)
  counted$fg$hs(problem$x0)
  expect_identical(
    counted$counts(),
    c(fn = 2L, gr = 1L, fg = 1L, hs = 1L)
  )
})

test_that("mize convergence agrees with independently evaluated quality", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.3.0")
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

  expect_null(run$conditions$error)
  expect_length(run$conditions$warnings, 0L)
  expect_identical(run$outcome$status, "converged")
  expect_identical(run$outcome$terminate$what, "ginf_tol")
  expect_true(run$outcome$converged)
  problem <- rosen()
  par <- run$parameters$final
  expect_equal(
    run$quality$final,
    list(
      objective = problem$fn(par),
      gradient_inf_norm = max(abs(problem$gr(par)))
    )
  )
  expect_equal(run$quality$final$objective, 0, tolerance = 1e-8)
  expect_lt(run$quality$final$gradient_inf_norm, 1e-6)
})

test_that("Hessian work retains native and physical count meanings", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.3.0")
  harness <- load_mize_harness()

  run <- harness$run_mize_problem(
    "rosen",
    method = "NEWTON",
    controls = list(
      line_search = "constant",
      step0 = 0.25,
      max_iter = 2,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )
  )

  expect_null(run$conditions$error)
  expect_identical(run$work$native$nh, 2)
  expect_identical(run$work$native$nhi, 0)
  expect_identical(run$work$physical[["hs"]], 2L)
})

test_that("hard budgets preserve absent native quality without fabrication", {
  testthat::skip_if_not_installed("mize", minimum_version = "0.3.0")
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
  testthat::skip_if_not_installed("mize", minimum_version = "0.3.0")
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
    name = c("ex_rosen", "jenn_samp"),
    n = c(4L, NA_integer_),
    m = c(NA_integer_, 20L)
  )
  suite <- harness$run_mize_suite(
    manifest,
    controls = list(max_iter = 2)
  )

  expect_identical(names(suite$runs), manifest$name)
  expect_identical(
    suite$runs$ex_rosen$problem$configuration$requested,
    list(n = 4L, m = NULL)
  )
  expect_identical(
    suite$runs$jenn_samp$problem$configuration$requested,
    list(n = NULL, m = 20L)
  )
  expect_identical(suite$runs$jenn_samp$reference$status, "not_applicable")
})
