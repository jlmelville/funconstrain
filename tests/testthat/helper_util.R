gfd <- function(par, fn, rel_eps = .Machine$double.eps^(1 / 3)) {
  g <- rep(0, length(par))
  for (i in seq_along(par)) {
    oldx <- par[i]
    h <- rel_eps * max(1, abs(oldx))
    par[i] <- oldx + h
    fplus <- fn(par)

    par[i] <- oldx - h
    fminus <- fn(par)
    par[i] <- oldx

    g[i] <- (fplus - fminus) / (2 * h)
  }
  g
}

make_gfd <- function(fn, rel_eps = .Machine$double.eps^(1 / 3)) {
  function(par) {
    gfd(par, fn, rel_eps)
  }
}

# test analytical gradient equals finite difference gradient at par
expect_gfd <- function(testfun, par, tolerance = 1e-6, tol = NULL) {
  if (!is.null(tol)) {
    tolerance <- tol
  }

  fd <- make_gfd(testfun$fn)(par)
  an <- testfun$gr(par)

  expect_equal(fd, an, tolerance = tolerance)
}

problem_factory_names <- function() {
  funconstrain_catalog()$name
}

capture_error_message <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = conditionMessage
  )
}

expect_all_checks <- function(checks) {
  results <- unlist(checks, recursive = TRUE, use.names = TRUE)
  stopifnot(
    is.logical(results),
    !is.null(names(results)),
    all(nzchar(names(results)))
  )

  failed <- names(results)[is.na(results) | !results]
  expect_identical(failed, character())
}

problem_factory_core_fields <- function() {
  c("fn", "gr", "he", "fg", "x0", "fmin", "xmin")
}

get_problem_factory <- function(name) {
  getExportedValue("funconstrain", name)
}

standard_x0 <- function(testfun) {
  if (is.function(testfun$x0)) {
    return(testfun$x0())
  }
  testfun$x0
}

hessian_fd <- function(par, gr, rel_eps = .Machine$double.eps^(1 / 4)) {
  hessian <- matrix(0, nrow = length(par), ncol = length(par))

  for (i in seq_along(par)) {
    oldx <- par[i]
    h <- rel_eps * max(1, abs(oldx))

    par[i] <- oldx + h
    gplus <- gr(par)

    par[i] <- oldx - h
    gminus <- gr(par)
    par[i] <- oldx

    hessian[, i] <- (gplus - gminus) / (2 * h)
  }

  hessian
}

expect_hfd <- function(
  testfun,
  par,
  tolerance = 1e-6,
  rel_eps = .Machine$double.eps^(1 / 4)
) {
  analytic <- testfun$he(par)
  n <- length(par)
  expect_equal(
    list(
      value = analytic,
      is_matrix = is.matrix(analytic),
      is_numeric = is.numeric(analytic),
      dimensions = dim(analytic),
      is_symmetric = is.matrix(analytic) && isSymmetric(analytic)
    ),
    list(
      value = hessian_fd(par, testfun$gr, rel_eps),
      is_matrix = TRUE,
      is_numeric = TRUE,
      dimensions = c(n, n),
      is_symmetric = TRUE
    ),
    tolerance = tolerance
  )
}

factory_m_metadata <- function() {
  do.call(
    rbind,
    lapply(problem_factory_names(), function(name) {
      testfun <- get_problem_factory(name)()
      if ("m" %in% names(testfun)) {
        return(data.frame(
          factory = name,
          m_present = TRUE,
          m_type = typeof(testfun$m),
          m_value = paste(testfun$m, collapse = ","),
          stringsAsFactors = FALSE
        ))
      }

      data.frame(
        factory = name,
        m_present = FALSE,
        m_type = NA_character_,
        m_value = NA_character_,
        stringsAsFactors = FALSE
      )
    })
  )
}
