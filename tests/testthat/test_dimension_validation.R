test_that("supported one-dimensional factories have coherent callbacks", {
  factories <- list(
    penalty_1 = penalty_1(),
    penalty_2 = penalty_2(),
    var_dim = var_dim(),
    trigon = trigon(),
    brown_al = brown_al(),
    disc_bv = disc_bv(),
    disc_ie = disc_ie(),
    broyden_tri = broyden_tri(),
    broyden_band = broyden_band(),
    linfun_fr = linfun_fr(m = 1),
    linfun_r1 = linfun_r1(m = 1),
    linfun_r1z = linfun_r1z(m = 1),
    chebyquad = chebyquad()
  )

  actual <- lapply(factories, function(testfun) {
    x <- if (is.function(testfun$x0)) testfun$x0(1) else testfun$x0
    fn <- testfun$fn(x)
    gr <- testfun$gr(x)
    he <- testfun$he(x)
    fg <- testfun$fg(x)

    list(
      fn_is_scalar = length(fn) == 1L,
      gr_is_scalar = length(gr) == 1L,
      he_is_matrix = is.matrix(he),
      he_dimensions_match = identical(dim(he), c(1L, 1L)),
      he_is_finite = all(is.finite(he)),
      fg_fn_matches = isTRUE(all.equal(fg$fn, fn)),
      fg_gr_matches = isTRUE(all.equal(fg$gr, gr))
    )
  })

  expect_all_checks(actual)
})

test_that("one-dimensional Brown Almost-Linear callbacks use the scalar formula", {
  testfun <- brown_al()
  x <- 0.3

  expect_equal(testfun$fn(x), (x - 1)^2)
  expect_equal(testfun$gr(x), 2 * (x - 1))
  expect_equal(testfun$he(x), matrix(2, 1, 1))
})

test_that("one-dimensional Discrete Boundary Value callbacks use one residual", {
  testfun <- disc_bv()
  x <- 0.3
  h <- 1 / 2
  t <- 1 / 2
  r <- 2 * x + (h^2 / 2) * (x + t + 1)^3
  rp <- 2 + 1.5 * h^2 * (x + t + 1)^2
  rpp <- 3 * h^2 * (x + t + 1)

  expect_equal(testfun$fn(x), r^2)
  expect_equal(testfun$gr(x), 2 * r * rp)
  expect_equal(testfun$he(x), matrix(2 * (rp^2 + r * rpp), 1, 1))
})

test_that("one-dimensional Broyden Tridiagonal callbacks use one residual", {
  testfun <- broyden_tri()
  x <- 0.3
  r <- (3 - 2 * x) * x + 1
  rp <- 3 - 4 * x

  expect_equal(testfun$fn(x), r^2)
  expect_equal(testfun$gr(x), 2 * r * rp)
  expect_equal(testfun$he(x), matrix(2 * (rp^2 - 4 * r), 1, 1))
})

test_that("one-dimensional Linear Rank 1 with Zero Columns and Rows is constant", {
  m_values <- c(1L, 2L, 3L)
  actual <- lapply(m_values, function(m) {
    testfun <- linfun_r1z(m)
    x <- testfun$x0(1)
    combined <- testfun$fg(x)

    list(
      fn = testfun$fn(x),
      gr = testfun$gr(x),
      he = testfun$he(x),
      fg_fn = combined$fn,
      fg_gr = combined$gr
    )
  })
  expected <- lapply(m_values, function(m) {
    list(fn = m, gr = 0, he = matrix(0, 1, 1), fg_fn = m, fg_gr = 0)
  })

  expect_equal(actual, expected)
})

test_that("structural variable-dimension rules are deliberate", {
  wat <- watson()
  expect_length(wat$x0(2), 2)
  expect_length(wat$x0(31), 31)
  expect_error(wat$x0(1), "outside")
  expect_error(wat$x0(32), "outside")

  exros <- ex_rosen()
  expect_length(exros$x0(2), 2)
  expect_length(exros$x0(4), 4)
  expect_error(exros$fn(numeric()), "outside")
  expect_error(exros$x0(3), "multiple")

  expow <- ex_powell()
  expect_length(expow$x0(4), 4)
  expect_length(expow$x0(8), 8)
  expect_error(expow$fn(numeric()), "outside")
  expect_error(expow$x0(2), "outside")
  expect_error(expow$x0(6), "multiple")
})

test_that("variable callbacks cover boundary dimensions", {
  cases <- list(
    list(name = "watson", n = c(2L, 6L, 3L)),
    list(name = "ex_rosen", n = c(2L, 8L, 4L)),
    list(name = "ex_powell", n = c(4L, 20L, 8L)),
    list(name = "penalty_1", n = c(1L, 25L, 2L)),
    list(name = "penalty_2", n = c(1L, 25L, 2L)),
    list(name = "var_dim", n = c(1L, 30L, 2L)),
    list(name = "trigon", n = c(1L, 30L, 2L)),
    list(name = "brown_al", n = c(1L, 30L, 2L)),
    list(name = "disc_bv", n = c(1L, 35L, 2L)),
    list(name = "disc_ie", n = c(1L, 35L, 2L)),
    list(name = "broyden_tri", n = c(1L, 40L, 2L)),
    list(name = "broyden_band", n = c(1L, 40L, 2L)),
    list(name = "linfun_fr", n = c(1L, 45L, 2L)),
    list(name = "linfun_r1", n = c(1L, 45L, 2L)),
    list(name = "linfun_r1z", n = c(1L, 45L, 2L)),
    list(name = "chebyquad", n = c(1L, 50L, 2L))
  )

  checks <- list()
  for (case in cases) {
    factory <- get_problem_factory(case$name)
    testfun <- if (grepl("^linfun_", case$name)) {
      factory(m = 100)
    } else {
      factory()
    }

    for (n in case$n) {
      key <- paste(case$name, "n =", n)
      x <- testfun$x0(n)
      fn <- testfun$fn(x)
      gr <- testfun$gr(x)
      he <- testfun$he(x)
      fg <- testfun$fg(x)

      checks[[key]] <- list(
        fn_is_finite_scalar = length(fn) == 1L && is.finite(fn),
        gr_length_matches = length(gr) == n,
        he_is_finite_matrix = is.matrix(he) && all(is.finite(he)),
        he_dimensions_match = identical(dim(he), c(n, n)),
        fg_fn_matches = isTRUE(all.equal(fg$fn, fn)),
        fg_gr_matches = isTRUE(all.equal(fg$gr, gr))
      )
    }
  }

  expect_all_checks(checks)
})

test_that("linear factories reject m less than n in every callback", {
  factories <- list(
    linfun_fr = linfun_fr,
    linfun_r1 = linfun_r1,
    linfun_r1z = linfun_r1z
  )
  actual <- character()

  for (name in names(factories)) {
    factory <- factories[[name]]
    testfun <- factory(m = 1)
    x <- c(0, 0)

    for (callback in c("fn", "gr", "he", "fg")) {
      key <- paste(name, callback)
      actual[[key]] <- capture_error_message(testfun[[callback]](x))
    }
    key <- paste(name, "x0")
    actual[[key]] <- capture_error_message(testfun$x0(2))
  }

  matches <- setNames(
    grepl("m must be >= n", actual, fixed = TRUE),
    names(actual)
  )
  expect_all_checks(matches)
})

test_that("factory m bounds are validated", {
  expect_silent(jenn_samp(2))
  expect_error(jenn_samp(1), "outside")

  expect_silent(gulf(3))
  expect_silent(gulf(100))
  expect_error(gulf(2), "outside")
  expect_error(gulf(101), "outside")

  expect_silent(box_3d(3))
  expect_error(box_3d(2), "outside")

  expect_silent(brown_den(4))
  expect_error(brown_den(3), "outside")

  expect_silent(biggs_exp6(6))
  expect_error(biggs_exp6(5), "outside")
})

test_that("accepted linear factory m values retain storage type", {
  factories <- list(
    linfun_fr = linfun_fr,
    linfun_r1 = linfun_r1,
    linfun_r1z = linfun_r1z
  )
  actual <- lapply(factories, function(factory) {
    integer_m <- factory(3L)
    double_m <- factory(3)

    list(integer = integer_m$m, double = double_m$m)
  })
  expected <- rep(
    list(list(integer = 3L, double = 3)),
    length(factories)
  )
  names(expected) <- names(factories)

  expect_identical(actual, expected)
})
