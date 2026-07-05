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

  expect_identical(expand_ranges("1:3"), 1:3)
  expect_identical(expand_ranges("3:1"), 3:1)
  expect_identical(expand_ranges("5"), 5L)

  expect_identical(parse_test_integers("1, 3:5, 7"), c(1L, 3L, 4L, 5L, 7L))
  expect_identical(parse_test_integers("2:4, 4, 2"), c(2L, 3L, 4L, 4L, 2L))

  expect_identical(
    parse_methods('c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs")'),
    c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs")
  )
  expect_identical(parse_methods("no quoted methods"), character())
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

test_that("fufn returns optimx-ready data for a problem", {
  testthat::skip_if_not_installed("optimx")

  tfun <- fufn(1)

  expect_named(
    tfun,
    c(
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
  )
  expect_equal(tfun$npar, 2)
  expect_identical(tfun$fname, "rosen")
  expect_true(is.function(tfun$fffn))
  expect_true(is.function(tfun$ffgr))
  expect_true(is.function(tfun$ffhe))
  expect_equal(tfun$x0, rosen()$x0)
  expect_true("L-BFGS-B" %in% tfun$ameth)
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
