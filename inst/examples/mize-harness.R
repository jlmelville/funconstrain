# Sourceable helpers for small funconstrain and mize integration checks.

as_mize_fg <- function(problem) {
  callbacks <- list(
    fn = problem$fn,
    gr = problem$gr,
    fg = problem$fg,
    hs = problem$he
  )
  valid <- vapply(callbacks, is.function, logical(1L))
  if (!all(valid)) {
    stop(
      "problem must provide function callbacks named fn, gr, fg, and he",
      call. = FALSE
    )
  }

  callbacks
}

count_mize_callbacks <- function(fg) {
  callback_names <- c("fn", "gr", "fg", "hs")
  counts <- new.env(parent = emptyenv())
  for (callback_name in callback_names) {
    counts[[callback_name]] <- 0L
  }

  wrap <- function(callback_name) {
    callback <- fg[[callback_name]]
    if (!is.function(callback)) {
      return(callback)
    }

    force(callback_name)
    force(callback)
    function(...) {
      counts[[callback_name]] <- counts[[callback_name]] + 1L
      callback(...)
    }
  }

  counted <- fg
  for (callback_name in callback_names) {
    counted[[callback_name]] <- wrap(callback_name)
  }

  list(
    fg = counted,
    counts = function() {
      vapply(callback_names, function(name) counts[[name]], integer(1L))
    }
  )
}

run_mize_problem <- function(
  name,
  n = NULL,
  m = NULL,
  method = "BFGS",
  controls = list()
) {
  minimum_mize_version <- "0.2.5.9002"
  installed_mize <- requireNamespace("mize", quietly = TRUE)
  if (
    !installed_mize ||
      utils::packageVersion("mize") < minimum_mize_version
  ) {
    stop(
      "run_mize_problem() requires GitHub mize >= ",
      minimum_mize_version,
      "; install it with pak::pak(\"jlmelville/mize\")",
      call. = FALSE
    )
  }
  if (!is.list(controls)) {
    stop("controls must be a named list", call. = FALSE)
  }
  if (length(controls) > 0L) {
    control_names <- names(controls)
    valid_names <- !is.null(control_names) &&
      !anyNA(control_names) &&
      all(nzchar(control_names))
    if (!valid_names) {
      stop("controls must be a named list", call. = FALSE)
    }
    if (anyDuplicated(control_names)) {
      stop("controls must not contain duplicate names", call. = FALSE)
    }

    owned <- intersect(control_names, c("par", "fg", "method"))
    if (length(owned) > 0L) {
      stop(
        "controls must not supply harness-owned argument(s): ",
        paste(owned, collapse = ", "),
        call. = FALSE
      )
    }

    unknown <- setdiff(control_names, names(formals(mize::mize)))
    if (length(unknown) > 0L) {
      stop(
        "unknown mize control(s): ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
  }

  problem <- funconstrain::funconstrain_problem(name, n = n, m = m)
  callbacks <- as_mize_fg(problem)
  counted <- count_mize_callbacks(callbacks)
  initial_quality <- list(
    objective = problem$fn(problem$x0),
    gradient_inf_norm = max(abs(problem$gr(problem$x0)))
  )

  warnings <- character()
  error_message <- NULL
  result <- tryCatch(
    withCallingHandlers(
      do.call(
        mize::mize,
        c(
          list(par = problem$x0, fg = counted$fg, method = method),
          controls
        )
      ),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      error_message <<- conditionMessage(condition)
      NULL
    }
  )

  if (!is.null(result)) {
    required_result_fields <- c(
      "par",
      "status",
      "message",
      "converged",
      "terminate",
      "nf",
      "ng",
      "nh",
      "nhi"
    )
    missing_fields <- setdiff(required_result_fields, names(result))
    if (length(missing_fields) > 0L) {
      stop(
        "mize result is missing required field(s): ",
        paste(missing_fields, collapse = ", "),
        call. = FALSE
      )
    }
  }

  final_parameters <- if (is.null(result)) NULL else result$par
  final_quality <- if (is.null(final_parameters)) {
    list(objective = NULL, gradient_inf_norm = NULL)
  } else {
    list(
      objective = problem$fn(final_parameters),
      gradient_inf_norm = max(abs(problem$gr(final_parameters)))
    )
  }

  mize_description <- utils::packageDescription("mize")
  mize_commit <- mize_description[["RemoteSha"]]
  if (is.null(mize_commit)) {
    mize_commit <- mize_description[["GithubSHA1"]]
  }
  if (is.null(mize_commit)) {
    mize_commit <- NA_character_
  }

  list(
    problem = list(
      number = problem$problem$number,
      name = problem$problem$name,
      title = problem$problem$title,
      configuration = problem$configuration
    ),
    method = method,
    controls = controls,
    parameters = list(initial = problem$x0, final = final_parameters),
    quality = list(initial = initial_quality, final = final_quality),
    reference = problem$reference$fmin,
    outcome = list(
      status = if (is.null(result)) NULL else result$status,
      message = if (is.null(result)) NULL else result$message,
      converged = if (is.null(result)) NULL else result$converged,
      terminate = if (is.null(result)) NULL else result$terminate
    ),
    native_quality = list(
      objective = if (is.null(result)) NULL else result[["f"]],
      gradient_l2_norm = if (is.null(result)) NULL else result[["g2n"]],
      gradient_inf_norm = if (is.null(result)) NULL else result[["ginfn"]]
    ),
    conditions = list(warnings = warnings, error = error_message),
    work = list(
      native = list(
        nf = if (is.null(result)) NULL else result$nf,
        ng = if (is.null(result)) NULL else result$ng,
        nh = if (is.null(result)) NULL else result$nh,
        nhi = if (is.null(result)) NULL else result$nhi,
        iter = if (is.null(result)) NULL else result[["iter"]]
      ),
      physical = counted$counts()
    ),
    provenance = list(
      r_version = R.version.string,
      platform = R.version$platform,
      funconstrain_version = as.character(
        utils::packageVersion("funconstrain")
      ),
      mize_version = as.character(utils::packageVersion("mize")),
      mize_commit = mize_commit
    )
  )
}

run_mize_suite <- function(
  manifest,
  method = "BFGS",
  controls = list()
) {
  if (!is.data.frame(manifest)) {
    stop("manifest must be a data frame", call. = FALSE)
  }
  required_columns <- c("name", "n", "m")
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "manifest is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  optional_dimension <- function(value) {
    if (length(value) == 1L && is.na(value)) NULL else value
  }
  runs <- lapply(seq_len(nrow(manifest)), function(row) {
    run_mize_problem(
      name = as.character(manifest$name[[row]]),
      n = optional_dimension(manifest$n[[row]]),
      m = optional_dimension(manifest$m[[row]]),
      method = method,
      controls = controls
    )
  })
  names(runs) <- as.character(manifest$name)

  provenance <- if (length(runs) == 0L) NULL else runs[[1L]]$provenance
  list(
    manifest = manifest,
    method = method,
    controls = controls,
    provenance = provenance,
    runs = runs
  )
}
