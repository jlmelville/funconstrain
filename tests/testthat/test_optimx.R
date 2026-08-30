open_file_connections <- function(path) {
  connections <- showConnections(all = TRUE)
  matching <- connections[, "description"] == path
  as.integer(rownames(connections)[matching])
}

close_file_connections <- function(path) {
  for (id in open_file_connections(path)) {
    connection <- getConnection(id)
    if (isOpen(connection)) {
      close(connection)
    }
  }
}

test_that("optimx parser helpers expand ranges and method strings", {
  expand_ranges <- getFromNamespace("expand_ranges", "funconstrain")
  parse_test_integers <- getFromNamespace("parse_test_integers", "funconstrain")
  parse_methods <- getFromNamespace("parse_methods", "funconstrain")

  actual <- list(
    ascending_range = expand_ranges("1:3"),
    descending_range = expand_ranges("3:1"),
    scalar = expand_ranges("5"),
    mixed_integer_specification = parse_test_integers("1, 3:5, 7"),
    repeated_integer_specification = parse_test_integers("2:4, 4, 2"),
    methods = parse_methods(
      'c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs")'
    ),
    no_methods = parse_methods("no quoted methods")
  )
  expected <- list(
    ascending_range = 1:3,
    descending_range = 3:1,
    scalar = 5L,
    mixed_integer_specification = c(1L, 3L, 4L, 5L, 7L),
    repeated_integer_specification = c(2L, 3L, 4L, 4L, 2L),
    methods = c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs"),
    no_methods = character()
  )

  expect_identical(actual, expected)
})

test_that("fufn and fufnrun report missing optimx clearly", {
  testthat::local_mocked_bindings(
    optimx_available = function() FALSE,
    .package = "funconstrain"
  )

  expect_error(
    fufn(1),
    "optimx package is required, please install it",
    fixed = TRUE
  )
  expect_error(
    fufnrun("does-not-exist.txt"),
    "optimx package is required, please install it",
    fixed = TRUE
  )
})

test_that("fufn returns each supported method once", {
  testthat::skip_if_not_installed("optimx")

  methods <- fufn(1)$ameth

  expect_identical(anyDuplicated(methods), 0L)
  expect_equal(sum(methods == "L-BFGS-B"), 1)
})

test_that("fufn dispatch table returns coherent data for every problem", {
  method_fixture <- c("nlminb", "Rvmmin", "L-BFGS-B")
  testthat::local_mocked_bindings(
    require_optimx = function() NULL,
    optimx_bounded_methods = function() method_fixture,
    .package = "funconstrain"
  )

  expected <- data.frame(
    fnum = 1:35,
    fname = c(
      "rosen",
      "freud_roth",
      "powell_bs",
      "brown_bs",
      "beale",
      "jenn_samp",
      "helical",
      "bard",
      "gauss",
      "meyer",
      "gulf",
      "box_3d",
      "powell_s",
      "wood",
      "kow_osb",
      "brown_den",
      "osborne_1",
      "biggs_exp6",
      "osborne_2",
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
    ),
    npar = c(
      2,
      2,
      2,
      2,
      2,
      2,
      3,
      3,
      3,
      3,
      3,
      3,
      4,
      4,
      4,
      4,
      5,
      6,
      11,
      8,
      10,
      20,
      10,
      10,
      6,
      8,
      8,
      6,
      8,
      8,
      8,
      8,
      8,
      8,
      8
    ),
    has_l_bfgs_b = c(rep(TRUE, 16), FALSE, rep(TRUE, 18)),
    stringsAsFactors = FALSE
  )

  expected_fields <- c(
    "npar",
    "fffn",
    "ffgr",
    "ffhe",
    "x0",
    "lo",
    "up",
    "mask",
    "fname",
    "ameth"
  )

  checks <- lapply(seq_len(nrow(expected)), function(i) {
    case <- expected[i, ]
    tfun <- fufn(case$fnum)
    factory <- getExportedValue("funconstrain", case$fname)
    direct <- factory()
    direct_x0 <- if (is.function(direct$x0)) {
      direct$x0(case$npar)
    } else {
      direct$x0
    }
    hessian <- tfun$ffhe(tfun$x0)
    expected_lo <- rep(min(direct_x0) - 0.1, case$npar)
    expected_up <- rep(max(direct_x0) + 0.1, case$npar)
    if (case$fnum == 4L) {
      expected_lo[] <- -1e20
      expected_up[] <- 1e20
    }
    if (case$fnum == 17L) {
      expected_lo[4:5] <- 0
    }
    expected_methods <- if (case$fnum == 17L) {
      method_fixture[method_fixture != "L-BFGS-B"]
    } else {
      method_fixture
    }

    list(
      fields_match = identical(names(tfun), expected_fields),
      fname_matches = identical(tfun$fname, case$fname),
      npar_matches = isTRUE(all.equal(tfun$npar, case$npar)),
      npar_is_double = identical(typeof(tfun$npar), "double"),
      mask_is_integer = identical(typeof(tfun$mask), "integer"),
      x0_matches = isTRUE(all.equal(tfun$x0, direct_x0)),
      lower_matches = isTRUE(all.equal(tfun$lo, expected_lo)),
      upper_matches = isTRUE(all.equal(tfun$up, expected_up)),
      methods_match = identical(tfun$ameth, expected_methods),
      x0_length_matches = length(tfun$x0) == tfun$npar,
      lower_length_matches = length(tfun$lo) == tfun$npar,
      upper_length_matches = length(tfun$up) == tfun$npar,
      mask_length_matches = length(tfun$mask) == tfun$npar,
      bounds_contain_x0 = all(tfun$lo <= tfun$x0) &&
        all(tfun$x0 <= tfun$up),
      callbacks_are_functions = all(vapply(
        tfun[c("fffn", "ffgr", "ffhe")],
        is.function,
        logical(1L)
      )),
      fn_is_scalar = length(tfun$fffn(tfun$x0)) == 1L,
      fn_matches = isTRUE(all.equal(
        tfun$fffn(tfun$x0),
        direct$fn(direct_x0)
      )),
      gr_length_matches = length(tfun$ffgr(tfun$x0)) == tfun$npar,
      gr_matches = isTRUE(all.equal(
        tfun$ffgr(tfun$x0),
        direct$gr(direct_x0)
      )),
      hessian_is_matrix = is.matrix(hessian),
      hessian_dimensions_match = isTRUE(all.equal(
        dim(hessian),
        c(tfun$npar, tfun$npar)
      )),
      l_bfgs_b_matches = identical(
        "L-BFGS-B" %in% tfun$ameth,
        case$has_l_bfgs_b
      )
    )
  })
  names(checks) <- expected$fname

  expect_all_checks(checks)
})

test_that("fufn retains the factories' variable-dimension callbacks", {
  testthat::local_mocked_bindings(
    require_optimx = function() NULL,
    optimx_bounded_methods = function() "L-BFGS-B",
    .package = "funconstrain"
  )

  tfun <- fufn(21)
  alternate_x <- ex_rosen()$x0(8)

  expect_equal(tfun$npar, 10)
  expect_equal(tfun$fffn(alternate_x), ex_rosen()$fn(alternate_x))
  expect_equal(tfun$ffgr(alternate_x), ex_rosen()$gr(alternate_x))
  expect_equal(tfun$ffhe(alternate_x), ex_rosen()$he(alternate_x))
})

test_that("fufnrun closes file and sink resources on parser errors", {
  testthat::skip_if_not_installed("optimx")

  rfo <- tempfile("bad-rfo-")
  sink_file <- tempfile("bad-rfo-sink-")
  writeLines(c(sink_file, "36", 'c("L-BFGS-B")', "FALSE"), rfo)

  start_sinks <- sink.number()
  on.exit(
    {
      while (sink.number() > start_sinks) {
        sink()
      }
      close_file_connections(rfo)
    },
    add = TRUE
  )

  capture.output(
    err <- tryCatch(fufnrun(rfo), error = identity)
  )
  leaked_connections <- open_file_connections(rfo)
  leaked_sink_depth <- sink.number()
  close_file_connections(rfo)
  while (sink.number() > start_sinks) {
    sink()
  }

  expect_s3_class(err, "error")
  expect_match(
    conditionMessage(err),
    "Problem number out of range",
    fixed = TRUE
  )
  expect_identical(leaked_connections, integer())
  expect_equal(leaked_sink_depth, start_sinks)
})

test_that("fufnrun closes file and sink resources on normal completion", {
  testthat::skip_if_not_installed("optimx")

  rfo <- tempfile("rfo-")
  sink_file <- tempfile("rfo-sink-")
  writeLines(c(sink_file, "1", 'c("L-BFGS-B")', "FALSE"), rfo)

  start_sinks <- sink.number()
  on.exit(
    {
      while (sink.number() > start_sinks) {
        sink()
      }
      close_file_connections(rfo)
    },
    add = TRUE
  )

  expect_no_error(
    capture.output(suppressWarnings(fufnrun(rfo)))
  )

  expect_identical(open_file_connections(rfo), integer())
  expect_equal(sink.number(), start_sinks)
  expect_true(file.exists(sink_file))
  expect_true(any(grepl("END : rosen", readLines(sink_file), fixed = TRUE)))
})
