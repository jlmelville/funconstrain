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

test_that("rejected RFO configurations never create or overwrite output", {
  testthat::local_mocked_bindings(
    require_optimx = function() NULL,
    fufn = function(...) stop("execution reached"),
    .package = "funconstrain"
  )
  valid <- c("1", 'c("L-BFGS-B")', "FALSE")
  cases <- list(
    missing_bounds = valid[1:2],
    bad_bounds = c(valid[1:2], "TURUE"),
    extra_line = c(valid, "unexpected"),
    empty_methods = c("1", 'c("")', "FALSE"),
    malformed_methods = c("1", 'junk "L-BFGS-B"', "FALSE")
  )
  for (ids in c(
    "",
    "1.9",
    "1:4:6",
    "1,,2",
    "1,",
    "0",
    "36",
    "1:99999999999999999999"
  )) {
    cases[[paste0("ids=", ids)]] <- c(ids, valid[2:3])
  }
  sentinel <- charToRaw("distinctive output\r\nwith exact bytes\n")
  checks <- list()
  for (case_name in names(cases)) {
    for (existing in c(FALSE, TRUE)) {
      rfo <- tempfile("rejected-rfo-")
      output <- tempfile("preserved-output-")
      if (existing) {
        writeBin(sentinel, output)
      }
      writeLines(c(output, cases[[case_name]]), rfo)
      depth <- sink.number()
      capture.output(message <- capture_error_message(fufnrun(rfo)))
      key <- paste(case_name, if (existing) "existing output" else "new output")
      checks[[key]] <- c(
        rejected_before_execution = !is.na(message) &&
          message != "execution reached",
        output_preserved = if (existing) {
          file.exists(output) &&
            identical(readBin(output, "raw", n = 1000), sentinel)
        } else {
          !file.exists(output)
        },
        sink_restored = sink.number() == depth,
        input_closed = length(open_file_connections(rfo)) == 0L
      )
    }
  }
  expect_all_checks(checks)
})

test_that("valid RFO whitespace, descending ranges, sorting and duplicates survive", {
  testthat::skip_if_not_installed("optimx")
  seen <- integer()
  methods <- NULL
  testthat::local_mocked_bindings(
    fufn = function(fnum) {
      seen <<- c(seen, fnum)
      list(
        fname = "fixture",
        x0 = 1,
        fffn = function(x) x^2,
        ffgr = function(x) 2 * x,
        ffhe = function(x) matrix(2)
      )
    },
    .package = "funconstrain"
  )
  testthat::local_mocked_bindings(
    opm = function(..., method) {
      methods <<- method
      data.frame(value = 0)
    },
    .package = "optimx"
  )
  rfo <- tempfile("valid-rfo-")
  output <- tempfile("valid-output-")
  writeLines(
    c(output, " 3 : 2, 1, 4 : 5, 2 ", ' c( "BFGS", "CG", "BFGS" ) ', " FALSE "),
    rfo
  )
  capture.output(fufnrun(rfo))
  expect_identical(seen, 1:5)
  expect_identical(methods, c("BFGS", "CG"))
})

test_that("a solver failure after valid configuration cleans up the sink", {
  testthat::local_mocked_bindings(
    require_optimx = function() NULL,
    fufn = function(...) stop("solver setup failed"),
    .package = "funconstrain"
  )
  rfo <- tempfile("solver-error-rfo-")
  output <- tempfile("partial-log-")
  writeLines(c(output, "1", 'c("BFGS")', "FALSE"), rfo)
  depth <- sink.number()
  capture.output(expect_error(fufnrun(rfo), "solver setup failed"))
  expect_equal(sink.number(), depth)
  expect_identical(open_file_connections(rfo), integer())
  expect_true(file.exists(output))
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

test_that("fufn validates problem numbers before requiring optimx", {
  malformed <- list(
    missing = NA_integer_,
    character = "1",
    vector = c(1, 2),
    fractional = 1.5,
    infinite = Inf,
    logical = TRUE
  )
  actual <- c(
    null = capture_error_message(fufn(NULL)),
    vapply(
      malformed,
      function(fnum) {
        capture_error_message(fufn(fnum))
      },
      character(1L)
    ),
    below_range = capture_error_message(fufn(0)),
    above_range = capture_error_message(fufn(36))
  )
  expected <- c(
    null = "ffn needs a function number fnum",
    setNames(
      rep(
        "fufn: fnum must be a finite whole-number scalar",
        length(malformed)
      ),
      names(malformed)
    ),
    below_range = "fufn: fnum is outside the allowed range",
    above_range = "fufn: fnum is outside the allowed range"
  )

  expect_identical(actual, expected)
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

  catalog <- funconstrain_catalog()
  expected_n <- setNames(catalog$n_default, catalog$name)
  expected_n[c(
    "watson",
    "ex_rosen",
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
  )] <- c(8L, 10L, 10L, 10L, 6L, 8L, 8L, 6L, rep(8L, 7))

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

  checks <- lapply(seq_len(nrow(catalog)), function(i) {
    name <- catalog$name[[i]]
    n <- expected_n[[name]]
    tfun <- fufn(catalog$number[[i]])
    direct <- get_problem_factory(name)()
    direct_x0 <- if (is.function(direct$x0)) {
      direct$x0(n)
    } else {
      direct$x0
    }
    expected_lo <- rep(min(direct_x0) - 0.1, n)
    expected_up <- rep(max(direct_x0) + 0.1, n)
    if (name == "brown_bs") {
      expected_lo[] <- -1e20
      expected_up[] <- 1e20
    }
    if (name == "osborne_1") {
      expected_lo[4:5] <- 0
    }
    expected_methods <- if (name == "osborne_1") {
      method_fixture[method_fixture != "L-BFGS-B"]
    } else {
      method_fixture
    }

    list(
      fields_match = identical(names(tfun), expected_fields),
      fname_matches = identical(tfun$fname, name),
      npar_matches = identical(tfun$npar, as.double(n)),
      x0_matches = identical(tfun$x0, direct_x0),
      lower_matches = identical(tfun$lo, expected_lo),
      upper_matches = identical(tfun$up, expected_up),
      mask_matches = identical(tfun$mask, rep(1L, n)),
      methods_match = identical(tfun$ameth, expected_methods),
      fn_matches = isTRUE(all.equal(tfun$fffn(tfun$x0), direct$fn(direct_x0))),
      gr_matches = isTRUE(all.equal(tfun$ffgr(tfun$x0), direct$gr(direct_x0))),
      he_matches = isTRUE(all.equal(tfun$ffhe(tfun$x0), direct$he(direct_x0)))
    )
  })
  names(checks) <- catalog$name

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
