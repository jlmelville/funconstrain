#' Manufactured Nonlinear Springs Problem
#'
#' Constructs a spring network on an interior `grid_size` by `grid_size` grid
#' with eliminated zero boundary values. If `D` is the oriented edge-difference
#' operator divided by the grid spacing, `z = D x`, and `s = D x_*`, the
#' grid spacing is `h = 1 / (grid_size + 1)` on the unit square, the objective is
#' \deqn{h^2 \sum_e (z_e-s_e)^2\left\{\frac{1}{2}+
#' \frac{\gamma}{4}\left[(z_e+s_e)^2+2s_e^2\right]\right\}.}
#' Thus the generated target `x_*` is a known global minimizer, and `gamma = 0`
#' gives the harmonic member of the family.
#'
#' Before normalization, the target at interior grid point `(i, j)` is
#' \deqn{\sin(\pi i h)\sin(\pi j h) +
#' 0.1\sin(3\pi i h)\sin(2\pi j h).}
#' It is scaled so that `max(abs(D x_*))` is one. The start is a seeded Gaussian
#' mixture of the discrete sine basis. Its coefficient for modes `(p, q)` is
#' divided by \eqn{(p^2+q^2)^{1/2}}, then the resulting field is also scaled so
#' that its maximum absolute edge difference is one.
#'
#' @param grid_size Number of interior points along each grid direction. Must
#'   be a whole number of at least three.
#' @param gamma Nonlinearity strength. Must be finite and nonnegative.
#' @param seed Whole-number random seed in `[0, .Machine$integer.max]`.
#'
#' @return A list with the core problem fields `fn`, `gr`, `he`, `fg`, `x0`,
#'   `fmin`, and `xmin`. Parameters use column-major interior-grid order. The
#'   Hessian callback returns a dense matrix and therefore uses quadratic
#'   storage in the parameter dimension.
#'
#'   The additional `quality(par)` callback returns `relative_field_error`,
#'   \eqn{\lVert x-x_*\rVert_2/\lVert x_*\rVert_2}, and
#'   `relative_edge_error`,
#'   \eqn{\lVert D(x-x_*)\rVert_2/\lVert D x_*\rVert_2}. Check both errors:
#'   small field error need not imply small edge error. `configuration` records
#'   the concrete
#'   constructor controls and parameter dimension `n`. `data` contains the grid
#'   spacing `h` and oriented `edge_from` and `edge_to` indices; zero denotes an
#'   eliminated boundary value. These values are inspectable snapshots rather
#'   than setters for the callbacks.
#'
#'   Construction uses a local Mersenne-Twister/Inversion/Rejection stream and
#'   restores the caller's `RNGkind()` settings and `.Random.seed`, including an
#'   absent seed. This cannot restore hidden cached state used by Box-Muller or
#'   custom user generators; use the default Inversion normal generator when an
#'   uninterrupted ambient normal stream is required. Exact generated values
#'   can depend on the R and numerical-library environment.
#'
#' @examples
#' problem <- nonlinear_springs(grid_size = 4, gamma = 10)
#' problem$fn(problem$x0)
#' problem$quality(problem$xmin)
#'
#' @export
nonlinear_springs <- function(grid_size = 8, gamma = 10, seed = 101) {
  problem <- "Nonlinear Springs"
  grid_size <- validate_family_whole(
    grid_size,
    problem,
    "grid_size",
    min = 3L
  )
  gamma <- validate_family_scalar(gamma, problem, "gamma", min = 0)
  seed <- validate_family_whole(seed, problem, "seed", min = 0L)
  n <- validate_family_product(
    grid_size,
    grid_size,
    problem,
    "parameter dimension"
  )
  edge_count <- validate_family_product(
    2L * grid_size,
    grid_size + 1L,
    problem,
    "edge dimension"
  )
  validate_family_product(n, 4L, problem, "incident dimension")
  h <- 1 / (grid_size + 1)

  index <- function(i, j) i + (j - 1L) * grid_size
  edge_from <- integer(edge_count)
  edge_to <- integer(edge_count)
  edge <- 0L
  for (j in seq_len(grid_size)) {
    for (i in 0L:grid_size) {
      edge <- edge + 1L
      edge_from[edge] <- if (i == 0L) 0L else index(i, j)
      edge_to[edge] <- if (i == grid_size) 0L else index(i + 1L, j)
    }
  }
  for (i in seq_len(grid_size)) {
    for (j in 0L:grid_size) {
      edge <- edge + 1L
      edge_from[edge] <- if (j == 0L) 0L else index(i, j)
      edge_to[edge] <- if (j == grid_size) 0L else index(i, j + 1L)
    }
  }

  incident_edge <- matrix(0L, nrow = n, ncol = 4L)
  incident_sign <- matrix(0, nrow = n, ncol = 4L)
  incident_count <- integer(n)
  for (current_edge in seq_len(edge_count)) {
    if (edge_from[current_edge] > 0L) {
      node <- edge_from[current_edge]
      incident_count[node] <- incident_count[node] + 1L
      incident_edge[node, incident_count[node]] <- current_edge
      incident_sign[node, incident_count[node]] <- -1
    }
    if (edge_to[current_edge] > 0L) {
      node <- edge_to[current_edge]
      incident_count[node] <- incident_count[node] + 1L
      incident_edge[node, incident_count[node]] <- current_edge
      incident_sign[node, incident_count[node]] <- 1
    }
  }

  apply_d <- function(par) {
    from_value <- numeric(edge_count)
    to_value <- numeric(edge_count)
    from_interior <- edge_from > 0L
    to_interior <- edge_to > 0L
    from_value[from_interior] <- par[edge_from[from_interior]]
    to_value[to_interior] <- par[edge_to[to_interior]]
    (to_value - from_value) / h
  }
  apply_dt <- function(edge_value) {
    incident_value <- matrix(edge_value[incident_edge], nrow = n, ncol = 4L)
    rowSums(incident_sign * incident_value) / h
  }

  grid <- seq_len(grid_size)
  sine_basis <- outer(grid, grid, function(i, mode) sin(pi * i * mode * h))
  target_field <- outer(
    grid,
    grid,
    function(i, j) {
      sin(pi * i * h) *
        sin(pi * j * h) +
        0.1 * sin(3 * pi * i * h) * sin(2 * pi * j * h)
    }
  )
  target_field <- as.vector(target_field)
  target_scale <- 1 / max(abs(apply_d(target_field)))
  xmin <- target_scale * target_field
  target_edges <- apply_d(xmin)

  gaussian_coefficients <- with_family_rng(
    seed,
    matrix(stats::rnorm(n), nrow = grid_size, ncol = grid_size)
  )
  filtered_coefficients <- gaussian_coefficients /
    sqrt(outer(grid^2, grid^2, `+`))
  start_field <- sine_basis %*% filtered_coefficients %*% t(sine_basis)
  start_field <- as.vector(start_field)
  start_scale <- 1 / max(abs(apply_d(start_field)))
  x0 <- start_scale * start_field

  as_parameter <- function(par) {
    validate_family_parameter(par, n, problem)
  }
  fn <- function(par) {
    par <- as_parameter(par)
    current_edges <- apply_d(par)
    difference <- current_edges - target_edges
    h^2 *
      sum(
        difference^2 *
          (0.5 +
            gamma / 4 * ((current_edges + target_edges)^2 + 2 * target_edges^2))
      )
  }
  gr <- function(par) {
    par <- as_parameter(par)
    current_edges <- apply_d(par)
    difference <- current_edges - target_edges
    h^2 *
      apply_dt(
        difference *
          (1 +
            gamma *
              (current_edges^2 + current_edges * target_edges + target_edges^2))
      )
  }
  fg <- function(par) {
    par <- as_parameter(par)
    current_edges <- apply_d(par)
    difference <- current_edges - target_edges
    list(
      fn = h^2 *
        sum(
          difference^2 *
            (0.5 +
              gamma /
                4 *
                ((current_edges + target_edges)^2 + 2 * target_edges^2))
        ),
      gr = h^2 *
        apply_dt(
          difference *
            (1 +
              gamma *
                (current_edges^2 +
                  current_edges * target_edges +
                  target_edges^2))
        )
    )
  }
  he <- function(par) {
    par <- as_parameter(par)
    validate_family_hessian_size(n, problem)
    current_edges <- apply_d(par)
    weights <- 1 + 3 * gamma * current_edges^2
    hessian <- matrix(0, nrow = n, ncol = n)
    for (current_edge in seq_len(edge_count)) {
      from <- edge_from[current_edge]
      to <- edge_to[current_edge]
      weight <- weights[current_edge]
      if (from > 0L) hessian[from, from] <- hessian[from, from] + weight
      if (to > 0L) hessian[to, to] <- hessian[to, to] + weight
      if (from > 0L && to > 0L) {
        hessian[from, to] <- hessian[from, to] - weight
        hessian[to, from] <- hessian[to, from] - weight
      }
    }
    hessian
  }
  target_field_norm <- sqrt(sum(xmin^2))
  target_edge_norm <- sqrt(sum(target_edges^2))
  quality <- function(par) {
    par <- as_parameter(par)
    error <- par - xmin
    list(
      relative_field_error = sqrt(sum(error^2)) / target_field_norm,
      relative_edge_error = sqrt(sum(apply_d(error)^2)) / target_edge_norm
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
      grid_size = grid_size,
      gamma = gamma,
      seed = seed,
      n = n
    ),
    data = list(
      h = h,
      edge_from = edge_from,
      edge_to = edge_to
    )
  )
}
