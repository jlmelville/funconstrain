test_that("validate_dimension preserves valid numeric inputs", {
  validate <- getFromNamespace("validate_dimension", "funconstrain")

  expect_identical(
    validate(4L, "test", min = 2L, max = 4L, multiple = 2L),
    4L
  )
  expect_identical(
    validate(4, "test", min = 2L, max = 4L, multiple = 2L),
    4
  )
  expect_error(
    validate(3, "test", min = 2L, max = 4L, multiple = 2L),
    "multiple"
  )
  expect_error(validate(1, "test", min = 2L), "outside")
  expect_error(validate(5, "test", max = 4L), "outside")
})

test_that("validate_dimension rejects malformed dimension values", {
  validate <- getFromNamespace("validate_dimension", "funconstrain")
  invalid <- list(
    1.5,
    Inf,
    NaN,
    TRUE,
    "1",
    numeric(),
    c(1, 2)
  )

  for (value in invalid) {
    expect_error(
      validate(value, "test"),
      "finite whole-number scalar"
    )
  }
  expect_error(validate(0, "test"), "outside")
  expect_error(validate(-1, "test"), "outside")
})

test_that("supported one-dimensional factories have coherent callbacks", {
  factories <- list(
    penalty_1(),
    penalty_2(),
    var_dim(),
    trigon(),
    brown_al(),
    disc_bv(),
    disc_ie(),
    broyden_tri(),
    broyden_band(),
    linfun_fr(m = 1),
    linfun_r1(m = 1),
    linfun_r1z(m = 1),
    chebyquad()
  )

  for (testfun in factories) {
    x <- if (is.function(testfun$x0)) testfun$x0(1) else testfun$x0
    fn <- testfun$fn(x)
    gr <- testfun$gr(x)
    he <- testfun$he(x)
    fg <- testfun$fg(x)

    expect_length(fn, 1)
    expect_length(gr, 1)
    expect_true(is.matrix(he))
    expect_equal(dim(he), c(1, 1))
    expect_true(all(is.finite(he)))
    expect_equal(fg$fn, fn)
    expect_equal(fg$gr, gr)
  }
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
  for (m in c(1L, 2L, 3L)) {
    testfun <- linfun_r1z(m)
    x <- testfun$x0(1)

    expect_equal(testfun$fn(x), m)
    expect_equal(testfun$gr(x), 0)
    expect_equal(testfun$he(x), matrix(0, 1, 1))
    expect_equal(testfun$fg(x)$fn, m)
    expect_equal(testfun$fg(x)$gr, 0)
  }
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

  for (case in cases) {
    factory <- get_problem_factory(case$name)
    testfun <- if (grepl("^linfun_", case$name)) {
      factory(m = 100)
    } else {
      factory()
    }

    for (n in case$n) {
      label <- paste(case$name, "n =", n)
      x <- testfun$x0(n)
      fn <- testfun$fn(x)
      gr <- testfun$gr(x)
      he <- testfun$he(x)
      fg <- testfun$fg(x)

      expect_true(length(fn) == 1 && is.finite(fn), info = label)
      expect_true(length(gr) == n, info = label)
      expect_true(
        is.matrix(he) && identical(dim(he), c(n, n)) && all(is.finite(he)),
        info = label
      )
      expect_equal(fg$fn, fn, info = label)
      expect_equal(fg$gr, gr, info = label)
    }
  }
})

test_that("linear factories reject m less than n in every callback", {
  for (factory in list(linfun_fr, linfun_r1, linfun_r1z)) {
    testfun <- factory(m = 1)
    x <- c(0, 0)

    for (callback in c("fn", "gr", "he", "fg")) {
      expect_error(testfun[[callback]](x), "m must be >= n")
    }
    expect_error(testfun$x0(2), "m must be >= n")
  }
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
  for (factory in list(linfun_fr, linfun_r1, linfun_r1z)) {
    integer_m <- factory(3L)
    double_m <- factory(3)

    expect_identical(integer_m$m, 3L)
    expect_identical(double_m$m, 3)
  }
})
