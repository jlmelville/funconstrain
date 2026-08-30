#' Discover funconstrain test problems
#'
#' `funconstrain_catalog()` lists the 35 Moré--Garbow--Hillstrom test problems
#' in paper order. `funconstrain_problem()` resolves one exact factory name and
#' an optional dimension request to a concrete, solver-neutral problem.
#'
#' @param name A single exact factory name from the `name` column of
#'   `funconstrain_catalog()`.
#' @param n Optional number of parameters. `NULL` uses the problem default.
#' @param m Optional number of summand functions. `NULL` uses the problem
#'   default or the value derived from `n`.
#'
#' @return `funconstrain_catalog()` returns a data frame with one row per
#'   problem and these columns:
#'
#'   - `number`, `name`, and `title` identify the problem;
#'   - `n_kind` is either `"fixed"` or `"variable"`;
#'   - `n_default` is the default parameter dimension;
#'   - `m_configurable` reports whether `m` is a factory argument; and
#'   - `m_default` is the effective `m` at `n_default`.
#'
#'   `funconstrain_problem()` returns a plain list with `fn`, `gr`, `he`, and
#'   `fg` callbacks; a concrete numeric `x0`; and `problem`, `configuration`,
#'   and `reference` metadata. The callbacks accept only vectors whose length
#'   equals `configuration$effective$n`.
#'
#'   `configuration` records normalized explicit requests under `requested`
#'   (omitted values remain `NULL`) and concrete integer values under
#'   `effective`.
#'
#'   Each `fmin` and `xmin` reference record contains `value`, `status`,
#'   `reason`, and `source_configuration`. Status is one of `"applicable"`,
#'   `"not_applicable"`, `"unavailable"`, or `"unknown"`. A stored value is
#'   returned even when it is not applicable so that its documented source
#'   configuration remains inspectable; unavailable values are `NULL`.
#'   Reasons currently distinguish documented applicability, a source
#'   configuration mismatch, no single stored minimizer, and an unencoded
#'   applicability rule. `source_configuration` is a partial rule: `n` and `m`
#'   are exact values, `_min` and `_max` suffixes are inclusive bounds, and an
#'   omitted component applies to every package-supported value.
#'
#' @examples
#' catalog <- funconstrain_catalog()
#' catalog[c(1, 20, 35), ]
#'
#' fixed <- funconstrain_problem("rosen")
#' fixed$fn(fixed$x0)
#'
#' variable_n <- funconstrain_problem("chebyquad", n = 10)
#' length(variable_n$x0)
#'
#' configurable_m <- funconstrain_problem("jenn_samp", m = 20)
#' configurable_m$configuration$effective
#'
#' different_references <- funconstrain_problem("ex_rosen", n = 4)
#' different_references$reference$fmin$status
#' different_references$reference$xmin$status
#'
#' @references
#' Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#' Testing unconstrained optimization software.
#' *ACM Transactions on Mathematical Software (TOMS)*, *7*(1), 17-41.
#' <https://doi.org/10.1145/355934.355936>
#'
#' @export
funconstrain_catalog <- function() {
  manifest <- funconstrain_problem_manifest()

  data.frame(
    number = manifest$number,
    name = manifest$name,
    title = manifest$title,
    n_kind = manifest$n_kind,
    n_default = manifest$n_default,
    m_configurable = manifest$m_kind == "configurable",
    m_default = manifest$m_default,
    stringsAsFactors = FALSE
  )
}

#' @rdname funconstrain_catalog
#' @export
funconstrain_problem <- function(name, n = NULL, m = NULL) {
  spec <- find_problem_spec(name)
  requested_n <- if (is.null(n)) NULL else resolve_n(n, spec)
  effective_n <- if (is.null(requested_n)) spec$n_default else requested_n
  requested_m <- if (is.null(m)) NULL else resolve_m(m, effective_n, spec)
  effective_m <- if (is.null(requested_m)) {
    default_m(effective_n, spec)
  } else {
    requested_m
  }

  validate_m_n_relationship(effective_m, effective_n, spec)

  factory <- get(
    spec$name,
    envir = environment(funconstrain_problem),
    inherits = FALSE
  )
  testfun <- if (spec$m_kind == "configurable") {
    factory(m = effective_m)
  } else {
    factory()
  }
  x0 <- if (is.function(testfun$x0)) testfun$x0(effective_n) else testfun$x0
  validate_resolved_x0(x0, effective_n, spec$title)

  list(
    fn = pin_problem_callback(testfun$fn, effective_n, spec$title),
    gr = pin_problem_callback(testfun$gr, effective_n, spec$title),
    he = pin_problem_callback(testfun$he, effective_n, spec$title),
    fg = pin_problem_callback(testfun$fg, effective_n, spec$title),
    x0 = x0,
    problem = list(
      number = spec$number,
      name = spec$name,
      title = spec$title
    ),
    configuration = list(
      requested = list(n = requested_n, m = requested_m),
      effective = list(n = effective_n, m = effective_m)
    ),
    reference = list(
      fmin = resolve_reference_record(
        testfun$fmin,
        spec$fmin_reference[[1L]],
        effective_n,
        effective_m
      ),
      xmin = resolve_reference_record(
        testfun$xmin,
        spec$xmin_reference[[1L]],
        effective_n,
        effective_m
      )
    )
  )
}

find_problem_spec <- function(name) {
  valid_name <- is.character(name) &&
    length(name) == 1L &&
    !is.na(name) &&
    nzchar(name)
  if (!valid_name) {
    stop("name must be a single non-missing character string")
  }

  manifest <- funconstrain_problem_manifest()
  index <- match(name, manifest$name)
  if (is.na(index)) {
    stop(
      "Unknown funconstrain problem '",
      name,
      "'; use an exact name from funconstrain_catalog()"
    )
  }

  manifest[index, , drop = FALSE]
}

resolve_n <- function(n, spec) {
  validate_resolver_dimension(
    n,
    problem = spec$title,
    min = spec$n_min,
    max = spec$n_max,
    multiple = spec$n_multiple,
    label = "n"
  )
}

resolve_m <- function(m, n, spec) {
  if (spec$m_kind == "configurable") {
    return(validate_resolver_dimension(
      m,
      problem = spec$title,
      min = spec$m_min,
      max = spec$m_max,
      multiple = spec$m_multiple,
      label = "m"
    ))
  }

  expected <- default_m(n, spec)
  validate_resolver_dimension(
    m,
    problem = spec$title,
    min = expected,
    max = expected,
    label = "m"
  )
}

default_m <- function(n, spec) {
  if (!is.na(spec$m_n_multiplier)) {
    value <- spec$m_n_multiplier * n + spec$m_n_offset
    return(validate_resolver_dimension(
      value,
      problem = spec$title,
      min = value,
      max = value,
      label = "m"
    ))
  }

  as.integer(spec$m_default)
}

validate_resolver_dimension <- function(
  value,
  problem,
  min,
  max,
  multiple = 1L,
  label
) {
  value <- validate_dimension(
    value,
    problem = problem,
    min = min,
    max = max,
    multiple = multiple,
    label = label
  )
  if (value > .Machine$integer.max) {
    stop(problem, ": ", label, " exceeds the maximum supported value")
  }

  as.integer(value)
}

validate_m_n_relationship <- function(m, n, spec) {
  if (spec$m_gte_n && m < n) {
    stop(spec$title, ": m must be greater than or equal to n")
  }
}

validate_resolved_x0 <- function(x0, n, problem) {
  valid <- is.numeric(x0) &&
    !is.complex(x0) &&
    length(x0) == n &&
    length(x0) > 0L &&
    all(is.finite(x0))
  if (!valid) {
    stop(
      "The resolved starting point for ",
      problem,
      " must be a finite numeric vector of length ",
      n
    )
  }

  invisible(x0)
}

pin_problem_callback <- function(callback, n, problem) {
  force(callback)
  force(n)
  force(problem)

  function(par) {
    if (length(par) != n) {
      stop(
        problem,
        ": resolved callbacks require a parameter vector of length ",
        n
      )
    }
    callback(par)
  }
}

resolve_reference_record <- function(value, rule, n, m) {
  source_configuration <- rule$source_configuration
  if (rule$availability == "unavailable") {
    return(list(
      value = NULL,
      status = "unavailable",
      reason = rule$reason,
      source_configuration = source_configuration
    ))
  }
  if (rule$availability == "unknown") {
    return(list(
      value = value,
      status = "unknown",
      reason = rule$reason,
      source_configuration = source_configuration
    ))
  }

  applicable <- reference_config_applies(source_configuration, n, m)
  list(
    value = value,
    status = if (applicable) "applicable" else "not_applicable",
    reason = if (applicable) {
      "documented_applicability"
    } else {
      "source_configuration_mismatch"
    },
    source_configuration = source_configuration
  )
}

reference_config_applies <- function(configuration, n, m) {
  if (length(configuration) == 0L) {
    return(TRUE)
  }

  current <- list(n = n, m = m)
  checks <- vapply(
    names(configuration),
    function(rule_name) {
      value <- configuration[[rule_name]]
      if (rule_name %in% c("n", "m")) {
        return(current[[rule_name]] == value)
      }

      component <- substr(rule_name, 1L, 1L)
      if (endsWith(rule_name, "_min")) {
        return(current[[component]] >= value)
      }
      if (endsWith(rule_name, "_max")) {
        return(current[[component]] <= value)
      }

      stop("Unknown reference configuration rule: ", rule_name)
    },
    logical(1L)
  )

  all(checks)
}
