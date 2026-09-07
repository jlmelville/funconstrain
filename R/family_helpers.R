validate_family_whole <- function(
  value,
  problem,
  label,
  min,
  max = .Machine$integer.max
) {
  if (!is.null(dim(value))) {
    stop(
      problem,
      ": ",
      label,
      " must be a finite whole-number scalar",
      call. = FALSE
    )
  }
  value <- validate_dimension(
    value,
    problem = problem,
    min = min,
    max = max,
    label = label
  )
  as.integer(value)
}

validate_family_scalar <- function(value, problem, label, min) {
  valid <- is.numeric(value) &&
    !is.logical(value) &&
    !is.complex(value) &&
    length(value) == 1L &&
    is.null(dim(value)) &&
    is.finite(value) &&
    value >= min
  if (!valid) {
    stop(
      problem,
      ": ",
      label,
      " must be a finite numeric scalar greater than or equal to ",
      min,
      call. = FALSE
    )
  }
  as.double(value)
}

validate_family_product <- function(left, right, problem, label) {
  if (left > .Machine$integer.max %/% right) {
    stop(
      problem,
      ": ",
      label,
      " exceeds the maximum supported size",
      call. = FALSE
    )
  }
  as.integer(left * right)
}

validate_family_hessian_size <- function(n, problem) {
  if (n > .Machine$integer.max %/% n) {
    stop(
      problem,
      ": the dense Hessian exceeds the maximum supported size",
      call. = FALSE
    )
  }
  invisible(n)
}

validate_family_parameter <- function(par, n, problem) {
  par <- promote_integer(par)
  valid <- is.numeric(par) &&
    !is.logical(par) &&
    !is.complex(par) &&
    is.null(dim(par)) &&
    length(par) == n &&
    all(is.finite(par))
  if (!valid) {
    stop(
      problem,
      ": par must be a finite real numeric vector of length ",
      n,
      call. = FALSE
    )
  }
  par
}

with_family_rng <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit(
    {
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) {
        seed_env <- globalenv()
        seed_env[[".Random.seed"]] <- old_seed
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(code)
}

family_basis <- function(x) {
  decomposition <- qr(x, LAPACK = FALSE)
  basis <- qr.Q(decomposition, complete = FALSE)
  signs <- sign(diag(qr.R(decomposition)))
  signs[signs == 0] <- 1
  sweep(basis, 2L, signs, `*`)
}
