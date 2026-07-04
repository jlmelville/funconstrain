gfd <- function(par, fn, rel_eps = sqrt(.Machine$double.eps)) {
  g <- rep(0, length(par))
  for (i in 1:length(par)) {
    oldx <- par[i]
    if (oldx != 0) {
      eps <- oldx * rel_eps
    }
    else {
      eps <- 1e-3
    }
    par[i] <- oldx + eps
    fplus <- fn(par)

    par[i] <- oldx - eps
    fminus <- fn(par)
    par[i] <- oldx

    g[i] <- (fplus - fminus) / (2 * eps)
  }
  g
}

make_gfd <- function(fn, rel_eps = sqrt(.Machine$double.eps)) {
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
  c(
    "bard",
    "beale",
    "biggs_exp6",
    "box_3d",
    "brown_al",
    "brown_bs",
    "brown_den",
    "broyden_band",
    "broyden_tri",
    "chebyquad",
    "disc_bv",
    "disc_ie",
    "ex_powell",
    "ex_rosen",
    "freud_roth",
    "gauss",
    "gulf",
    "helical",
    "jenn_samp",
    "kow_osb",
    "linfun_fr",
    "linfun_r1",
    "linfun_r1z",
    "meyer",
    "osborne_1",
    "osborne_2",
    "penalty_1",
    "penalty_2",
    "powell_bs",
    "powell_s",
    "rosen",
    "trigon",
    "var_dim",
    "watson",
    "wood"
  )
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
    eps <- rel_eps
    if (oldx != 0) {
      eps <- abs(oldx) * rel_eps
    }

    par[i] <- oldx + eps
    gplus <- gr(par)

    par[i] <- oldx - eps
    gminus <- gr(par)
    par[i] <- oldx

    hessian[, i] <- (gplus - gminus) / (2 * eps)
  }

  hessian
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
