psd_dense_oracle <- function(problem, par) {
  size <- problem$configuration$matrix_size
  rank <- problem$configuration$rank
  basis <- problem$data$basis
  values <- problem$data$eigenvalues
  factor <- matrix(par, nrow = size, ncol = rank)
  target <- tcrossprod(sweep(basis, 2L, sqrt(values), `*`))
  residual <- tcrossprod(factor) - target
  n <- length(par)

  hessian <- vapply(
    seq_len(n),
    function(column) {
      direction <- matrix(0, nrow = size, ncol = rank)
      direction[column] <- 1
      as.vector(
        direction %*%
          crossprod(factor) +
          factor %*%
            (crossprod(direction, factor) + crossprod(factor, direction)) -
          target %*% direction
      )
    },
    numeric(n)
  )

  list(
    fn = sum(residual * residual) / 4,
    gr = as.vector(residual %*% factor),
    he = hessian
  )
}

psd_rel_norm <- function(actual, expected) {
  sqrt(sum((actual - expected)^2)) / max(1, sqrt(sum(expected^2)))
}

test_that("PSD factorization exposes its standalone contract", {
  problem <- psd_factorization()

  expect_named(
    problem,
    c(
      "fn",
      "gr",
      "he",
      "fg",
      "x0",
      "fmin",
      "xmin",
      "quality",
      "configuration",
      "data"
    ),
    ignore.order = FALSE
  )
  expect_identical(
    problem$configuration,
    list(matrix_size = 64L, rank = 4L, kappa = 10, seed = 101L, n = 256L)
  )
  expect_identical(dim(problem$data$basis), c(64L, 4L))
  expect_identical(length(problem$data$eigenvalues), 4L)
  expect_identical(length(problem$x0), 256L)
  expect_identical(length(problem$xmin), 256L)
  expect_identical(problem$fmin, 0)

  minimum <- psd_factorization(matrix_size = 2, rank = 2)
  next_size <- psd_factorization(matrix_size = 3, rank = 2)
  expect_identical(minimum$configuration$n, 4L)
  expect_identical(next_size$configuration$n, 6L)

  for (current in list(minimum, next_size, problem)) {
    n <- current$configuration$n
    expect_true(is.finite(current$fn(current$x0)))
    expect_length(current$gr(current$x0), n)
    expect_identical(dim(current$he(current$x0)), c(n, n))
    expect_identical(
      current$fg(current$x0),
      list(fn = current$fn(current$x0), gr = current$gr(current$x0))
    )
    expect_equal(current$fn(current$xmin), current$fmin, tolerance = 1e-28)
    expect_named(
      current$quality(current$xmin),
      c("relative_product_error", "max_relative_mode_error")
    )
  }
})

test_that("PSD value, gradient, and Hessian match a dense oracle", {
  problem <- psd_factorization(matrix_size = 5, rank = 3, kappa = 25, seed = 17)
  par <- problem$x0 + sin(seq_along(problem$x0)) / 13
  expected <- psd_dense_oracle(problem, par)

  expect_equal(problem$fn(par), expected$fn, tolerance = 1e-14)
  expect_equal(problem$gr(par), expected$gr, tolerance = 1e-13)
  expect_equal(problem$he(par), expected$he, tolerance = 1e-13)
  expect_true(isSymmetric(problem$he(par)))
  expect_identical(
    problem$fg(par),
    list(fn = problem$fn(par), gr = problem$gr(par))
  )
})

test_that("PSD derivatives remain accurate over multiple scales", {
  problem <- psd_factorization(matrix_size = 4, rank = 2, kappa = 10, seed = 3)
  direction <- cos(seq_along(problem$x0))
  direction <- direction / sqrt(sum(direction^2))

  for (scale in c(0.01, 1, 10)) {
    par <- scale * problem$x0 + sin(seq_along(problem$x0)) / 20
    steps <- 10^(-4:-6) * max(1, sqrt(sum(par^2)))
    value_errors <- numeric(length(steps))
    gradient_errors <- numeric(length(steps))
    for (i in seq_along(steps)) {
      step <- steps[i]
      fd_value <- (problem$fn(par + step * direction) -
        problem$fn(par - step * direction)) /
        (2 * step)
      fd_gradient <- (problem$gr(par + step * direction) -
        problem$gr(par - step * direction)) /
        (2 * step)
      value_errors[i] <- abs(fd_value - sum(problem$gr(par) * direction)) /
        max(1, abs(fd_value))
      gradient_errors[i] <- psd_rel_norm(
        fd_gradient,
        as.vector(problem$he(par) %*% direction)
      )
    }

    expect_lt(max(value_errors), 2e-6)
    expect_lt(max(gradient_errors), 2e-6)
  }
})

test_that("PSD invariances and weak-mode diagnostics are visible", {
  problem <- psd_factorization(matrix_size = 6, rank = 2, kappa = 100, seed = 9)
  factor <- matrix(problem$xmin, nrow = 6L, ncol = 2L)
  angle <- 0.37
  rotation <- matrix(
    c(cos(angle), sin(angle), -sin(angle), cos(angle)),
    nrow = 2L
  )
  rotated <- as.vector(factor %*% rotation)

  expect_equal(problem$fn(rotated), problem$fn(problem$xmin), tolerance = 1e-28)
  expect_equal(
    problem$quality(rotated),
    problem$quality(problem$xmin),
    tolerance = 1e-13
  )
  expect_lt(problem$fn(problem$xmin), 1e-28)
  expect_lt(max(abs(problem$gr(problem$xmin))), 1e-14)

  skew <- matrix(c(0, 1, -1, 0), nrow = 2L)
  tangent <- as.vector(factor %*% skew)
  expect_lt(max(abs(problem$he(problem$xmin) %*% tangent)), 1e-13)

  omitted <- factor
  omitted[, 2L] <- 0
  omitted_quality <- problem$quality(as.vector(omitted))
  expected_product <- problem$data$eigenvalues[2L] /
    sqrt(sum(problem$data$eigenvalues^2))
  expect_equal(
    omitted_quality$relative_product_error,
    expected_product,
    tolerance = 1e-13
  )
  expect_equal(omitted_quality$max_relative_mode_error, 1, tolerance = 1e-13)
})

test_that("PSD generation is seeded, local, and configuration-sensitive", {
  first <- psd_factorization(matrix_size = 5, rank = 2, seed = 41)
  repeated <- psd_factorization(matrix_size = 5, rank = 2, seed = 41)
  changed <- psd_factorization(matrix_size = 5, rank = 2, seed = 42)
  conditioned <- psd_factorization(
    matrix_size = 5,
    rank = 2,
    kappa = 100,
    seed = 41
  )

  expect_identical(first$x0, repeated$x0)
  expect_identical(first$xmin, repeated$xmin)
  expect_false(identical(first$x0, changed$x0))
  expect_false(identical(first$xmin, changed$xmin))
  expect_identical(first$data$basis, conditioned$data$basis)
  expect_false(identical(first$data$eigenvalues, conditioned$data$eigenvalues))

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed)
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
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

  RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
  set.seed(811)
  saved_kind <- RNGkind()
  saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  psd_factorization(matrix_size = 3, rank = 2, seed = 7)
  expect_identical(RNGkind(), saved_kind)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), saved_seed)

  rm(".Random.seed", envir = .GlobalEnv)
  psd_factorization(matrix_size = 3, rank = 2, seed = 7)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(812)
  stats::rnorm(1)
  stats::runif(1)
  expected_normal <- stats::rnorm(1)
  expected_uniform <- stats::runif(1)
  set.seed(812)
  stats::rnorm(1)
  stats::runif(1)
  psd_factorization(matrix_size = 3, rank = 2, seed = 7)
  expect_identical(stats::rnorm(1), expected_normal)
  expect_identical(stats::runif(1), expected_uniform)
})

test_that("PSD validates controls and callback parameters", {
  bad_size <- list(1, 2.5, NA_real_, Inf, TRUE, 1 + 1i, "4", matrix(4))
  for (value in bad_size) {
    expect_error(psd_factorization(matrix_size = value), "matrix_size")
  }
  expect_error(psd_factorization(matrix_size = 4, rank = 1), "rank")
  expect_error(psd_factorization(matrix_size = 4, rank = 5), "rank")
  expect_error(psd_factorization(matrix_size = 4, rank = matrix(2)), "rank")
  expect_error(psd_factorization(kappa = 0.9), "kappa")
  expect_error(psd_factorization(kappa = Inf), "kappa")
  expect_error(psd_factorization(kappa = matrix(2)), "kappa")
  expect_error(psd_factorization(seed = -1), "seed")
  expect_error(psd_factorization(seed = 1.5), "seed")
  expect_error(psd_factorization(seed = matrix(1)), "seed")
  expect_error(
    psd_factorization(matrix_size = 50000, rank = 2),
    "target matrix dimension"
  )

  problem <- psd_factorization(matrix_size = 3, rank = 2)
  integer_par <- seq_len(problem$configuration$n)
  numeric_par <- as.numeric(integer_par)
  for (callback in c("fn", "gr", "he", "fg", "quality")) {
    expect_equal(
      problem[[callback]](integer_par),
      problem[[callback]](numeric_par)
    )
  }

  malformed <- list(
    numeric(5),
    numeric(7),
    rep(NA_real_, 6),
    rep(Inf, 6),
    rep(TRUE, 6),
    rep(1 + 0i, 6),
    matrix(0, nrow = 2, ncol = 3),
    letters[1:6]
  )
  for (callback in c("fn", "gr", "he", "fg", "quality")) {
    for (par in malformed) {
      expect_error(
        problem[[callback]](par),
        "finite real numeric vector of length 6"
      )
    }
  }
})

test_that("PSD generated problems survive serialization", {
  problem <- psd_factorization(matrix_size = 4, rank = 2, seed = 57)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(problem, path)
  restored <- readRDS(path)

  expect_identical(restored$x0, problem$x0)
  expect_identical(restored$xmin, problem$xmin)
  expect_identical(restored$data, problem$data)
  expect_identical(restored$fg(restored$x0), problem$fg(problem$x0))
  expect_identical(
    restored$quality(restored$xmin),
    problem$quality(problem$xmin)
  )
})
