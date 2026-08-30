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
  if ((fnum < 1) || (fnum > 35)) {
    stop("fnum must be in [1, 35]")
  }
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

  mycon <- file(filename, open = "r", blocking = TRUE)
  mycon_open <- TRUE
  on.exit(
    {
      if (mycon_open) {
        close(mycon)
      }
    },
    add = TRUE
  )

  sink_open <- FALSE
  on.exit(
    {
      if (sink_open) {
        sink()
      }
    },
    add = TRUE
  )

  sfname <- readLines(mycon, n = 1)
  if (length(sfname) == 0) {
    cat("no sink file\n")
  } else {
    cat("opening sink file ", sfname, "\n")
    sink(sfname, split = TRUE)
    sink_open <- TRUE
  } # open sink file
  cat("sink file name=", sfname, "\n")

  lin2 <- readLines(mycon, n = 1)
  cat("probs =", lin2, "\n")
  if (length(lin2) == 0) {
    stop("Unexpected null probs")
  }
  probc <- parse_test_integers(lin2)

  # ?? should we check it worked?
  cat("Problem numbers:\n")
  print(probc)
  print(unique(probc))
  if (length(unique(probc)) < length(probc)) {
    cat("Duplicated problem numbers, simplifying\n")
    probc <- unique(probc)
  }
  probc <- sort(probc)
  cat("Final problem numbers:\n")
  print(probc)
  # check loop
  for (iprob in probc) {
    # loop over problems
    if ((iprob < 1) || (iprob > 35)) {
      stop("Problem number out of range. Stopping.")
    }
  } # end check loop
  meths <- readLines(mycon, n = 1)
  if (length(meths) == 0) {
    stop("Unexpected null meths")
  }
  cat("Methods:\n")
  cat(meths, "\n")

  methc <- parse_methods(meths)
  for (package in c("lbfgs", "lbfgsb3c")) {
    if (package %in% methc) {
      if (!requireNamespace(package, quietly = TRUE)) {
        stop(
          paste(package, "package is required, please install it"),
          call. = FALSE
        )
      }
    }
  }

  if (length(unique(methc)) < length(methc)) {
    cat("Duplicated methods, simplifying\n")
    methc <- unique(methc)
  }
  cat("methods in list form:")
  print(methc)
  tbounds <- readLines(mycon, n = 1)
  have_bounds <- FALSE
  if (tbounds == "TRUE") {
    have_bounds <- TRUE
  }
  cat("have.bounds:", have_bounds, "\n")
  close(mycon)
  mycon_open <- FALSE
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
  if (grepl(":", x)) {
    range_boundaries <- as.integer(unlist(strsplit(x, ":")))
    return(seq(from = range_boundaries[1], to = range_boundaries[2]))
  }
  as.integer(x)
}

# parse the integers from a string like "1, 2, 3, 4:6" -> c(1, 2, 3, 4, 5, 6)
# the integers represent the test ids of the functions in this package
parse_test_integers <- function(test_fun_str) {
  elements <- strsplit(test_fun_str, ",\\s*")[[1]]
  unlist(lapply(elements, expand_ranges))
}

# parse the method line from the RFO.txt file e.g. a string like
# 'c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs")' to an actual R character vector
# without going through eval
parse_methods <- function(input) {
  matches <- gregexpr('"(.*?)"', input, perl = TRUE)
  if (matches[[1]][1] == -1) {
    return(character())
  }
  string_list <- regmatches(input, matches)[[1]]
  res <- sapply(string_list, function(x) substr(x, 2, nchar(x) - 1))
  names(res) <- NULL
  res
}
