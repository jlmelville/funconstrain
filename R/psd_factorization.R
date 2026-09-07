#' Symmetric Positive Semidefinite Factorization Problem
#'
#' Constructs the nonconvex factorization problem
#' \deqn{F(U) = \lVert UU^T - M \rVert_F^2 / 4,}
#' where \eqn{r} is `rank` and
#' \deqn{M = Q \mathop{\rm diag}(\lambda) Q^T, \qquad
#' \lambda_j = \kappa^{-(j-1)/(r-1)}.}
#' The columns of `Q` are the sign-fixed thin QR basis of a seeded Gaussian
#' matrix. The returned minimizer is
#' \deqn{U_* = Q \mathop{\rm diag}(\sqrt{\lambda}).}
#' The starting factor is an independent seeded Gaussian matrix scaled to have
#' squared Frobenius norm `sum(lambda)`. Parameters contain `U` in column-major
#' order. Multiplying `U` on the right by an orthogonal matrix leaves the
#' objective unchanged.
#'
#' @param matrix_size Number of rows in `U`. Must be a whole number.
#' @param rank Number of columns in `U`. Must be a whole number between 2 and
#'   `matrix_size`.
#' @param kappa Ratio of the largest to smallest positive target eigenvalue.
#'   Must be finite and at least one. This is not the condition number of the
#'   full Hessian, which is singular at a minimizer.
#' @param seed Whole-number random seed in `[0, .Machine$integer.max]`.
#'
#' @return A list with the core problem fields `fn`, `gr`, `he`, `fg`, `x0`,
#'   `fmin`, and `xmin`. The Hessian callback returns a dense matrix and
#'   therefore uses quadratic storage in the parameter dimension. The
#'   generated `xmin` is one of many equivalent minimizing factors.
#'
#'   The additional `quality(par)` callback returns:
#'
#'   - `relative_product_error`,
#'     \eqn{\lVert UU^T-M\rVert_F/\lVert\lambda\rVert_2}; and
#'   - `max_relative_mode_error`, the maximum over columns `q_j` of
#'     \eqn{\lVert(UU^T-M)q_j\rVert_2/\lambda_j}.
#'
#'   Product error is equivalent to objective value, while relative mode error
#'   separately exposes weak-mode recovery. `configuration` records the
#'   concrete constructor controls and parameter dimension `n`. `data`
#'   contains snapshots of the target `basis` and positive `eigenvalues`.
#'   Objective and quality evaluation form the dense `matrix_size`-squared
#'   residual, while gradient evaluation uses factored algebra without that
#'   dense residual.
#'   Construction uses a local Mersenne-Twister/Inversion/Rejection stream and
#'   restores the caller's `RNGkind()` settings and `.Random.seed`, including an
#'   absent seed. This cannot restore hidden cached state used by Box-Muller or
#'   custom user generators; use the default Inversion normal generator when an
#'   uninterrupted ambient normal stream is required. Exact generated values
#'   can depend on the R and numerical-library environment, particularly the QR
#'   factorization.
#'
#' @examples
#' problem <- psd_factorization(matrix_size = 8, rank = 2)
#' problem$fn(problem$x0)
#' problem$quality(problem$xmin)
#'
#' @export
psd_factorization <- function(
  matrix_size = 64,
  rank = 4,
  kappa = 10,
  seed = 101
) {
  problem <- "PSD Factorization"
  matrix_size <- validate_family_whole(
    matrix_size,
    problem,
    "matrix_size",
    min = 2L
  )
  rank <- validate_family_whole(
    rank,
    problem,
    "rank",
    min = 2L,
    max = matrix_size
  )
  kappa <- validate_family_scalar(kappa, problem, "kappa", min = 1)
  seed <- validate_family_whole(seed, problem, "seed", min = 0L)
  n <- validate_family_product(
    matrix_size,
    rank,
    problem,
    "parameter dimension"
  )
  validate_family_product(
    matrix_size,
    matrix_size,
    problem,
    "target matrix dimension"
  )

  generated <- with_family_rng(seed, {
    list(
      target = matrix(stats::rnorm(n), nrow = matrix_size, ncol = rank),
      start = matrix(stats::rnorm(n), nrow = matrix_size, ncol = rank)
    )
  })
  basis <- family_basis(generated$target)
  eigenvalues <- kappa^(-(seq_len(rank) - 1) / (rank - 1))
  target <- tcrossprod(sweep(basis, 2L, sqrt(eigenvalues), `*`))
  start_scale <- sqrt(sum(eigenvalues)) / sqrt(sum(generated$start^2))
  x0 <- as.vector(start_scale * generated$start)
  xmin <- as.vector(sweep(basis, 2L, sqrt(eigenvalues), `*`))

  as_factor <- function(par) {
    matrix(
      validate_family_parameter(par, n, problem),
      nrow = matrix_size,
      ncol = rank
    )
  }
  fn <- function(par) {
    factor <- as_factor(par)
    residual <- tcrossprod(factor) - target
    sum(residual * residual) / 4
  }
  gr <- function(par) {
    factor <- as_factor(par)
    target_factor <- basis %*%
      sweep(
        crossprod(basis, factor),
        1L,
        eigenvalues,
        `*`
      )
    as.vector(factor %*% crossprod(factor) - target_factor)
  }
  fg <- function(par) {
    factor <- as_factor(par)
    residual <- tcrossprod(factor) - target
    target_factor <- basis %*%
      sweep(
        crossprod(basis, factor),
        1L,
        eigenvalues,
        `*`
      )
    list(
      fn = sum(residual * residual) / 4,
      gr = as.vector(factor %*% crossprod(factor) - target_factor)
    )
  }
  he <- function(par) {
    factor <- as_factor(par)
    validate_family_hessian_size(n, problem)
    factor_gram <- crossprod(factor)
    residual <- tcrossprod(factor) - target
    hessian <- matrix(0, nrow = n, ncol = n)
    for (column_block in seq_len(rank)) {
      columns <- seq_len(matrix_size) + (column_block - 1L) * matrix_size
      for (row_block in seq_len(rank)) {
        rows <- seq_len(matrix_size) + (row_block - 1L) * matrix_size
        block <- tcrossprod(factor[, column_block], factor[, row_block])
        diag(block) <- diag(block) + factor_gram[column_block, row_block]
        if (row_block == column_block) {
          block <- block + residual
        }
        hessian[rows, columns] <- block
      }
    }
    hessian
  }
  quality <- function(par) {
    factor <- as_factor(par)
    residual <- tcrossprod(factor) - target
    mode_residual <- residual %*% basis
    list(
      relative_product_error = sqrt(sum(residual^2)) / sqrt(sum(eigenvalues^2)),
      max_relative_mode_error = max(
        sqrt(colSums(mode_residual^2)) / eigenvalues
      )
    )
  }

  list(
    fn = fn,
    gr = gr,
    he = he,
    fg = fg,
    x0 = x0,
    fmin = 0,
    xmin = xmin,
    quality = quality,
    configuration = list(
      matrix_size = matrix_size,
      rank = rank,
      kappa = kappa,
      seed = seed,
      n = n
    ),
    data = list(
      basis = basis,
      eigenvalues = eigenvalues
    )
  )
}
