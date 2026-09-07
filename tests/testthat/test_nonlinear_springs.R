spring_matrix <- function(problem) {
  n <- problem$configuration$n
  from <- problem$data$edge_from
  to <- problem$data$edge_to
  edge_count <- length(from)
  operator <- matrix(0, nrow = edge_count, ncol = n)

  for (edge in seq_len(edge_count)) {
    if (from[edge] > 0L) operator[edge, from[edge]] <- -1 / problem$data$h
    if (to[edge] > 0L) operator[edge, to[edge]] <- 1 / problem$data$h
  }
  operator
}

spring_dense_oracle <- function(problem, par) {
  operator <- spring_matrix(problem)
  gamma <- problem$configuration$gamma
  h <- problem$data$h
  current <- as.vector(operator %*% par)
  target <- as.vector(operator %*% problem$xmin)
  difference <- current - target
  potential <- function(edge_value) {
    edge_value^2 / 2 + gamma * edge_value^4 / 4
  }
  edge_gradient <- difference *
    (1 + gamma * (current^2 + current * target + target^2))
  weights <- 1 + 3 * gamma * current^2

  list(
    fn = h^2 *
      sum(
        potential(current) -
          potential(target) -
          (target + gamma * target^3) * difference
      ),
    gr = as.vector(h^2 * crossprod(operator, edge_gradient)),
    he = h^2 * crossprod(operator, weights * operator)
  )
}

spring_rel_norm <- function(actual, expected) {
  sqrt(sum((actual - expected)^2)) / max(1, sqrt(sum(expected^2)))
}

test_that("nonlinear springs exposes its standalone contract", {
  problem <- nonlinear_springs()

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
    list(grid_size = 8L, gamma = 10, seed = 101L, n = 64L)
  )
  expect_identical(problem$data$h, 1 / 9)
  expect_identical(length(problem$data$edge_from), 144L)
  expect_identical(length(problem$data$edge_to), 144L)
  expect_true(any(problem$data$edge_from == 0L))
  expect_true(any(problem$data$edge_to == 0L))
  expect_identical(length(problem$x0), 64L)
  expect_identical(length(problem$xmin), 64L)
  expect_identical(problem$fmin, 0)

  minimum <- nonlinear_springs(grid_size = 3)
  next_size <- nonlinear_springs(grid_size = 4)
  expect_identical(minimum$configuration$n, 9L)
  expect_identical(next_size$configuration$n, 16L)

  for (current in list(minimum, next_size, problem)) {
    n <- current$configuration$n
    expect_true(is.finite(current$fn(current$x0)))
    expect_length(current$gr(current$x0), n)
    expect_identical(dim(current$he(current$x0)), c(n, n))
    expect_identical(
      current$fg(current$x0),
      list(fn = current$fn(current$x0), gr = current$gr(current$x0))
    )
    expect_identical(current$fn(current$xmin), current$fmin)
    expect_named(
      current$quality(current$xmin),
      c("relative_field_error", "relative_edge_error")
    )
  }
})

test_that("spring value, gradient, and Hessian match an incidence oracle", {
  problem <- nonlinear_springs(grid_size = 4, gamma = 7, seed = 17)
  par <- problem$x0 + sin(seq_along(problem$x0)) / 13
  expected <- spring_dense_oracle(problem, par)

  expect_equal(problem$fn(par), expected$fn, tolerance = 1e-13)
  expect_equal(problem$gr(par), expected$gr, tolerance = 1e-13)
  expect_equal(problem$he(par), expected$he, tolerance = 1e-13)
  expect_true(isSymmetric(problem$he(par)))
  expect_identical(
    problem$fg(par),
    list(fn = problem$fn(par), gr = problem$gr(par))
  )
})

test_that("spring derivatives remain accurate over multiple scales", {
  problem <- nonlinear_springs(grid_size = 3, gamma = 10, seed = 3)
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
      gradient_errors[i] <- spring_rel_norm(
        fd_gradient,
        as.vector(problem$he(par) %*% direction)
      )
    }

    expect_lt(max(value_errors), 2e-6)
    expect_lt(max(gradient_errors), 2e-6)
  }
})

test_that("spring target, start, and quality diagnostics preserve structure", {
  harmonic <- nonlinear_springs(grid_size = 6, gamma = 0, seed = 29)
  nonlinear <- nonlinear_springs(grid_size = 6, gamma = 10, seed = 29)
  operator <- spring_matrix(harmonic)

  expect_identical(harmonic$x0, nonlinear$x0)
  expect_identical(harmonic$xmin, nonlinear$xmin)
  expect_identical(harmonic$data, nonlinear$data)
  expect_equal(max(abs(operator %*% harmonic$x0)), 1, tolerance = 1e-15)
  expect_equal(max(abs(operator %*% harmonic$xmin)), 1, tolerance = 1e-15)
  expect_identical(
    harmonic$quality(harmonic$xmin),
    list(
      relative_field_error = 0,
      relative_edge_error = 0
    )
  )
  expect_identical(harmonic$fn(harmonic$xmin), 0)
  expect_equal(harmonic$gr(harmonic$xmin), rep(0, 36L))

  grid <- seq_len(6L)
  h <- harmonic$data$h
  low <- as.vector(outer(grid, grid, function(i, j) {
    sin(pi * i * h) * sin(pi * j * h)
  }))
  high <- as.vector(outer(grid, grid, function(i, j) {
    sin(6 * pi * i * h) * sin(6 * pi * j * h)
  }))
  low <- low / sqrt(sum(low^2))
  high <- high / sqrt(sum(high^2))
  low_quality <- harmonic$quality(harmonic$xmin + 0.1 * low)
  high_quality <- harmonic$quality(harmonic$xmin + 0.1 * high)

  expect_equal(
    low_quality$relative_field_error,
    high_quality$relative_field_error,
    tolerance = 1e-15
  )
  expect_gt(
    high_quality$relative_edge_error,
    3 * low_quality$relative_edge_error
  )

  broad <- nonlinear_springs(grid_size = 8, gamma = 0, seed = 101)
  broad_grid <- seq_len(8L)
  broad_basis <- outer(broad_grid, broad_grid, function(i, mode) {
    sin(pi * i * mode / 9)
  })
  broad_error <- matrix(broad$x0 - broad$xmin, nrow = 8L, ncol = 8L)
  coefficients <- (2 / 9) * crossprod(broad_basis, broad_error) %*% broad_basis
  values_1d <- 4 * sin(broad_grid * pi / 18)^2
  mode_values <- outer(values_1d, values_1d, `+`)
  energy <- 0.5 * mode_values * coefficients^2
  upper_half <- mode_values >= median(as.vector(mode_values))
  expect_gt(sum(energy[upper_half]) / sum(energy), 0.1)
})

test_that("harmonic springs have the discrete sine spectrum", {
  size <- 5L
  problem <- nonlinear_springs(grid_size = size, gamma = 0)
  hessian <- problem$he(problem$x0)
  grid <- seq_len(size)
  mode_p <- 2L
  mode_q <- 4L
  mode <- as.vector(outer(grid, grid, function(i, j) {
    sin(mode_p * pi * i / (size + 1)) * sin(mode_q * pi * j / (size + 1))
  }))
  value <- 4 *
    sin(mode_p * pi / (2 * (size + 1)))^2 +
    4 * sin(mode_q * pi / (2 * (size + 1)))^2

  expect_identical(diag(hessian), rep(4, size^2))
  expect_equal(as.vector(hessian %*% mode), value * mode, tolerance = 1e-13)
  expect_identical(problem$he(problem$xmin), hessian)
})

test_that("spring generation is seeded and restores caller RNG state", {
  first <- nonlinear_springs(grid_size = 4, seed = 41)
  repeated <- nonlinear_springs(grid_size = 4, seed = 41)
  changed <- nonlinear_springs(grid_size = 4, seed = 42)

  expect_identical(first$x0, repeated$x0)
  expect_identical(first$xmin, repeated$xmin)
  expect_false(identical(first$x0, changed$x0))
  expect_identical(first$xmin, changed$xmin)

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
  nonlinear_springs(grid_size = 3, seed = 7)
  expect_identical(RNGkind(), saved_kind)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), saved_seed)

  rm(".Random.seed", envir = .GlobalEnv)
  nonlinear_springs(grid_size = 3, seed = 7)
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
  nonlinear_springs(grid_size = 3, seed = 7)
  expect_identical(stats::rnorm(1), expected_normal)
  expect_identical(stats::runif(1), expected_uniform)
})

test_that("springs validates controls and callback parameters", {
  bad_size <- list(2, 3.5, NA_real_, Inf, TRUE, 1 + 1i, "4", matrix(4))
  for (value in bad_size) {
    expect_error(nonlinear_springs(grid_size = value), "grid_size")
  }
  expect_error(nonlinear_springs(gamma = -1), "gamma")
  expect_error(nonlinear_springs(gamma = Inf), "gamma")
  expect_error(nonlinear_springs(gamma = matrix(2)), "gamma")
  expect_error(nonlinear_springs(seed = -1), "seed")
  expect_error(nonlinear_springs(seed = 1.5), "seed")
  expect_error(nonlinear_springs(seed = matrix(1)), "seed")
  expect_error(nonlinear_springs(grid_size = 30000), "incident dimension")
  expect_error(nonlinear_springs(grid_size = 40000), "edge dimension")

  problem <- nonlinear_springs(grid_size = 3)
  integer_par <- seq_len(problem$configuration$n)
  numeric_par <- as.numeric(integer_par)
  for (callback in c("fn", "gr", "he", "fg", "quality")) {
    expect_equal(
      problem[[callback]](integer_par),
      problem[[callback]](numeric_par)
    )
  }

  malformed <- list(
    numeric(8),
    numeric(10),
    rep(NA_real_, 9),
    rep(Inf, 9),
    rep(TRUE, 9),
    rep(1 + 0i, 9),
    matrix(0, nrow = 3, ncol = 3),
    letters[1:9]
  )
  for (callback in c("fn", "gr", "he", "fg", "quality")) {
    for (par in malformed) {
      expect_error(
        problem[[callback]](par),
        "finite real numeric vector of length 9"
      )
    }
  }
})

test_that("spring generated problems survive serialization", {
  problem <- nonlinear_springs(grid_size = 4, gamma = 3, seed = 57)
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
