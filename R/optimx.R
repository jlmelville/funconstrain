#' Test Function Data for use with Optimx
#'
#' This function provides formatted data about each of the functions in
#' this package for ease of use with the `optimx` package.
#'
#' @param fnum Function number (1-35) to return data for.
#' @seealso [fufnrun()] to run optimx using this data.
#' @export
fufn <- function(fnum) {
  if (is.null(fnum)) {
    stop("ffn needs a function number fnum")
  }
  fnum <- validate_dimension(
    fnum,
    problem = "fufn",
    min = 1L,
    max = 35L,
    label = "fnum"
  )
  require_optimx()

  manifest <- funconstrain_problem_manifest()
  index <- match(fnum, manifest$number)
  if (is.na(index)) {
    stop("fnum must be in [1, 35]")
  }
  spec <- manifest[index, , drop = FALSE]

  fname <- spec$name
  n <- as.numeric(fufn_problem_n(spec))
  factory <- get(fname, envir = environment(fufn), inherits = FALSE)
  testfun <- factory()
  x0 <- if (is.function(testfun$x0)) testfun$x0(n) else testfun$x0

  methods <- optimx_bounded_methods()
  if (spec$number == 17L) {
    methods <- methods[methods != "L-BFGS-B"]
  }

  lower <- rep(min(x0) - 0.1, n)
  upper <- rep(max(x0) + 0.1, n)
  if (spec$number == 4L) {
    lower[] <- -1e20
    upper[] <- 1e20
  }
  if (spec$number == 17L) {
    lower[4:5] <- 0
  }

  list(
    npar = n,
    fffn = testfun$fn,
    ffgr = testfun$gr,
    ffhe = testfun$he,
    x0 = x0,
    lo = lower,
    up = upper,
    mask = rep(1L, n),
    fname = fname,
    ameth = methods
  )
}

fufn_problem_n <- function(spec) {
  # Preserve the smaller dimensions historically selected by fufn.
  legacy_overrides <- stats::setNames(
    c(8L, 10L, 10L, 10L, 6L, 8L, 8L, 6L, rep(8L, 7L)),
    c(20L, 21L, 23:35)
  )
  key <- as.character(spec$number)

  if (key %in% names(legacy_overrides)) {
    return(unname(legacy_overrides[[key]]))
  }

  spec$n_default
}

#' Run Optimx Using RFO Data
#'
#' Use the information in an RFO file to run methods from optimx on the
#' functions in this package and report the results.
#'
#' @param filename path to an RFO file contain information on the functions and
#' methods from optimx to use.
#' @export
fufnrun <- function(filename = "RFO.txt") {
  # fufnrun.R -- J C Nash 2024-4-8
  ## ?? fixing kkt
  # RFO.txt is input file
  require_optimx()

  config <- readLines(filename, warn = FALSE)
  if (length(config) < 4L || any(nzchar(trimws(config[-seq_len(4L)])))) {
    stop("RFO configuration must contain four lines", call. = FALSE)
  }
  sfname <- config[1L]
  if (!nzchar(trimws(sfname))) {
    stop("RFO output path must not be empty", call. = FALSE)
  }
  probc <- sort(unique(parse_test_integers(config[2L])))
  methc <- unique(parse_methods(config[3L]))
  tbounds <- trimws(config[4L])
  if (!tbounds %in% c("TRUE", "FALSE")) {
    stop("RFO bounds flag must be TRUE or FALSE", call. = FALSE)
  }
  have_bounds <- tbounds == "TRUE"

  for (package in c("lbfgs", "lbfgsb3c")) {
    if (package %in% methc && !requireNamespace(package, quietly = TRUE)) {
      stop(
        paste(package, "package is required, please install it"),
        call. = FALSE
      )
    }
  }

  # Rejected configuration must not create or overwrite the output file.
  sink_open <- FALSE
  on.exit(
    {
      if (sink_open) {
        sink()
      }
    },
    add = TRUE
  )
  cat("opening sink file ", sfname, "\n")
  sink(sfname, split = TRUE)
  sink_open <- TRUE
  cat("sink file name=", sfname, "\n")
  cat("Final problem numbers:\n")
  print(probc)
  cat("methods in list form:")
  print(methc)
  cat("have.bounds:", have_bounds, "\n")
  for (iprob in probc) {
    # loop over problems
    tfun <- fufn(fnum = iprob)
    # print(tfun)
    cat("Problem:", tfun$fname, "\n")
    x0 <- tfun$x0
    if (have_bounds) {
      lo <- tfun$lo
      up <- tfun$up
    } else {
      lo <- -Inf
      up <- Inf
    }
    tfn <- tfun$fffn
    attr(tfn, "fname") <- tfun$fname
    tgr <- tfun$ffgr
    the <- tfun$ffhe
    nx0 <- length(x0)
    #   cat("about to call opm\n")

    if (have_bounds) {
      t21 <- optimx::opm(
        x0,
        tfn,
        tgr,
        hess = the,
        lower = lo,
        upper = up,
        method = methc,
        control = list(trace = 0)
      )
    } else {
      t21 <-
        optimx::opm(
          x0,
          tfn,
          tgr,
          hess = the,
          method = methc,
          control = list(trace = 0)
        )
    }
    print(summary(t21, order = 'value', par.select = 1:min(nx0, 5)))
    cat("END :", tfun$fname, "\n\n")
  }
  if (sink_open) {
    sink()
    sink_open <- FALSE
  }
}

optimx_available <- function() {
  requireNamespace("optimx", quietly = TRUE)
}

optimx_bounded_methods <- function() {
  methods <- optimx::ctrldefault(2)$bdmeth
  methods <- methods[methods != "lbfgsb3c"]
  unique(c(methods, "L-BFGS-B"))
}

require_optimx <- function() {
  if (!optimx_available()) {
    stop("optimx package is required, please install it", call. = FALSE)
  }
}

# converts e.g. 1:3 -> 1, 2, 3
expand_ranges <- function(x) {
  if (!grepl("^\\s*[0-9]+\\s*(:\\s*[0-9]+\\s*)?$", x)) {
    stop("Invalid problem ID or range", call. = FALSE)
  }
  boundaries <- as.double(strsplit(x, ":", fixed = TRUE)[[1L]])
  if (any(!is.finite(boundaries) | boundaries < 1 | boundaries > 35)) {
    stop("Problem number out of range. Stopping.", call. = FALSE)
  }
  boundaries <- as.integer(boundaries)
  if (length(boundaries) == 2L) {
    return(seq.int(boundaries[1L], boundaries[2L]))
  }
  boundaries
}

# parse the integers from a string like "1, 2, 3, 4:6" -> c(1, 2, 3, 4, 5, 6)
# the integers represent the test ids of the functions in this package
parse_test_integers <- function(test_fun_str) {
  if (!nzchar(trimws(test_fun_str)) || grepl(",\\s*$", test_fun_str)) {
    stop("Problem selection must contain nonempty IDs or ranges", call. = FALSE)
  }
  elements <- strsplit(test_fun_str, ",", fixed = TRUE)[[1L]]
  unlist(lapply(elements, expand_ranges))
}

# parse the method line from the RFO.txt file e.g. a string like
# 'c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs")' to an actual R character vector
# without going through eval
parse_methods <- function(input) {
  input <- trimws(input)
  if (grepl("^c\\s*\\(", input)) {
    if (!grepl("\\)$", input)) {
      stop("Invalid methods specification", call. = FALSE)
    }
    input <- sub("^c\\s*\\(", "", input)
    input <- sub("\\)$", "", input)
  }
  if (!grepl('^\\s*"[^"]+"\\s*(,\\s*"[^"]+"\\s*)*$', input)) {
    stop("Invalid or empty methods specification", call. = FALSE)
  }
  matches <- gregexpr('"[^"]+"', input)
  string_list <- regmatches(input, matches)[[1]]
  res <- sapply(string_list, function(x) substr(x, 2, nchar(x) - 1))
  names(res) <- NULL
  if (any(!nzchar(trimws(res)))) {
    stop("Invalid or empty methods specification", call. = FALSE)
  }
  res
}
