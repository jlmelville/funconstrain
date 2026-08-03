test_that("factory contract list covers all exported problem factories", {
  exported_factories <- setdiff(
    getNamespaceExports("funconstrain"),
    c("fufn", "fufnrun")
  )

  expect_equal(sort(problem_factory_names()), sort(exported_factories))
  expect_equal(length(problem_factory_names()), 35)
})

test_that("problem factories expose documented core fields", {
  core_fields <- problem_factory_core_fields()

  for (name in problem_factory_names()) {
    testfun <- get_problem_factory(name)()

    expect_true(all(core_fields %in% names(testfun)), info = name)
    expect_true(is.function(testfun$fn), info = paste(name, "fn"))
    expect_true(is.function(testfun$gr), info = paste(name, "gr"))
    expect_true(is.function(testfun$he), info = paste(name, "he"))
    expect_true(is.function(testfun$fg), info = paste(name, "fg"))
    expect_true(is.numeric(testfun$fmin), info = paste(name, "fmin"))
    expect_equal(length(testfun$fmin), 1, info = paste(name, "fmin"))
    expect_true(length(testfun$xmin) > 0, info = paste(name, "xmin"))

    x0 <- standard_x0(testfun)
    expect_true(is.numeric(x0), info = paste(name, "x0"))
    expect_true(length(x0) > 0, info = paste(name, "x0"))
  }
})

test_that("fixed-dimensional callbacks reject wrong-length parameters", {
  cases <- list(
    list(name = "rosen", n = 2L, problem = "Rosenbrock"),
    list(name = "freud_roth", n = 2L, problem = "Freudenstein-Roth"),
    list(name = "powell_bs", n = 2L, problem = "Powell Badly Scaled"),
    list(name = "brown_bs", n = 2L, problem = "Brown Badly Scaled"),
    list(name = "beale", n = 2L, problem = "Beale"),
    list(name = "jenn_samp", n = 2L, problem = "Jennrich-Sampson"),
    list(name = "helical", n = 3L, problem = "Helical Valley"),
    list(name = "bard", n = 3L, problem = "Bard"),
    list(name = "gauss", n = 3L, problem = "Gaussian"),
    list(name = "meyer", n = 3L, problem = "Meyer"),
    list(name = "gulf", n = 3L, problem = "Gulf"),
    list(name = "box_3d", n = 3L, problem = "Box 3D"),
    list(name = "powell_s", n = 4L, problem = "Powell Singular"),
    list(name = "wood", n = 4L, problem = "Wood"),
    list(name = "kow_osb", n = 4L, problem = "Kowalik-Osborne"),
    list(name = "brown_den", n = 4L, problem = "Brown Dennis"),
    list(name = "osborne_1", n = 5L, problem = "Osborne 1"),
    list(name = "biggs_exp6", n = 6L, problem = "Biggs EXP6"),
    list(name = "osborne_2", n = 11L, problem = "Osborne 2")
  )

  for (case in cases) {
    testfun <- get_problem_factory(case$name)()
    for (callback in c("fn", "gr", "he", "fg")) {
      for (bad_n in c(case$n - 1L, case$n + 1L)) {
        expect_error(
          testfun[[callback]](rep(0, bad_n)),
          paste0(
            "^",
            case$problem,
            ": n is outside the allowed range$"
          ),
          info = paste(case$name, callback, "n =", bad_n)
        )
      }
    }
  }
})

test_that("fg matches fn and gr at each factory standard x0", {
  for (name in problem_factory_names()) {
    testfun <- get_problem_factory(name)()
    x0 <- standard_x0(testfun)
    fg <- testfun$fg(x0)

    expect_true(all(c("fn", "gr") %in% names(fg)), info = name)
    expect_equal(fg$fn, testfun$fn(x0), info = paste(name, "fn"))
    expect_equal(fg$gr, testfun$gr(x0), info = paste(name, "gr"))
  }
})

test_that("Hessians are numeric symmetric matrices at each factory standard x0", {
  for (name in problem_factory_names()) {
    testfun <- get_problem_factory(name)()
    x0 <- standard_x0(testfun)
    hessian <- testfun$he(x0)

    expect_true(is.matrix(hessian), info = name)
    expect_true(is.numeric(hessian), info = name)
    expect_equal(dim(hessian), c(length(x0), length(x0)), info = name)
    expect_true(isSymmetric(hessian), info = name)
  }
})

test_that("representative Hessians match finite differences of gradients", {
  for (name in c("bard", "beale", "rosen", "wood")) {
    testfun <- get_problem_factory(name)()
    x0 <- standard_x0(testfun)

    expect_equal(
      testfun$he(x0),
      hessian_fd(x0, testfun$gr),
      tolerance = 1e-4,
      info = name
    )
  }
})

test_that("reported xmin values evaluate to reported fmin values", {
  for (name in problem_factory_names()) {
    testfun <- get_problem_factory(name)()

    if (anyNA(testfun$xmin)) {
      next
    }

    actual <- testfun$fn(testfun$xmin)
    tolerance <- max(1e-8, abs(testfun$fmin) * 1e-5)

    expect_true(is.finite(actual), info = name)
    expect_true(
      abs(actual - testfun$fmin) <= tolerance,
      info = paste(
        name,
        "fn(xmin) =",
        format(actual, digits = 16),
        "fmin =",
        format(testfun$fmin, digits = 16)
      )
    )
  }
})

test_that("documented stored reference configurations are internally valid", {
  cases <- list(
    list(
      name = "jenn_samp",
      args = list(m = 10),
      configuration = "m=10",
      xmin_length = 2L
    ),
    list(
      name = "gulf",
      args = list(m = 10),
      configuration = "m=10",
      xmin_length = 3L
    ),
    list(
      name = "box_3d",
      args = list(m = 20),
      configuration = "m=20",
      xmin_length = 3L
    ),
    list(
      name = "brown_den",
      args = list(m = 20),
      configuration = "m=20",
      xmin_length = 4L
    ),
    list(
      name = "biggs_exp6",
      args = list(m = 13),
      configuration = "m=13",
      xmin_length = 6L
    ),
    list(
      name = "watson",
      args = list(),
      configuration = "n=6",
      xmin_length = 6L
    ),
    list(
      name = "ex_rosen",
      args = list(),
      configuration = "n=8",
      xmin_length = 8L
    ),
    list(
      name = "ex_powell",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "penalty_1",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "penalty_2",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "var_dim",
      args = list(),
      configuration = "n=6",
      xmin_length = 6L
    ),
    list(
      name = "trigon",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "brown_al",
      args = list(),
      configuration = "n=4",
      xmin_length = 4L
    ),
    list(
      name = "disc_bv",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "disc_ie",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "broyden_tri",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "broyden_band",
      args = list(),
      configuration = "n=5",
      xmin_length = 5L
    ),
    list(
      name = "linfun_fr",
      args = list(m = 100),
      configuration = "n=4,m=100",
      xmin_length = 4L
    ),
    list(
      name = "linfun_r1",
      args = list(m = 100),
      configuration = "n=5,m=100",
      xmin_length = 5L
    ),
    list(
      name = "linfun_r1z",
      args = list(m = 100),
      configuration = "n=5,m=100",
      xmin_length = 5L
    ),
    list(
      name = "chebyquad",
      args = list(),
      configuration = "n=8",
      xmin_length = 8L
    )
  )

  for (case in cases) {
    test_that(
      paste("stored reference is valid for", case$name, case$configuration),
      {
        testfun <- do.call(get_problem_factory(case$name), case$args)
        info <- paste(case$name, case$configuration)

        expect_equal(
          length(testfun$xmin),
          case$xmin_length,
          info = paste(info, "xmin length")
        )

        if (!anyNA(testfun$xmin)) {
          actual <- testfun$fn(testfun$xmin)
          tolerance <- max(1e-8, abs(testfun$fmin) * 1e-5)

          expect_true(is.finite(actual), info = paste(info, "fn(xmin)"))
          expect_true(
            abs(actual - testfun$fmin) <= tolerance,
            info = paste(
              info,
              "fn(xmin) =",
              format(actual, digits = 16),
              "fmin =",
              format(testfun$fmin, digits = 16)
            )
          )
        }
      }
    )
  }
})

test_that("linfun_r1z handles empty interiors at n=1 and n=2", {
  testfun <- linfun_r1z(m = 100)

  for (n in 1:2) {
    x <- rep(0, n)
    info <- paste("linfun_r1z", "m=100", paste0("n=", n))

    expect_equal(testfun$fn(x), 100, info = paste(info, "fn"))
    expect_equal(testfun$gr(x), rep(0, n), info = paste(info, "gr"))
    expect_equal(testfun$fg(x)$fn, 100, info = paste(info, "fg fn"))
    expect_equal(testfun$fg(x)$gr, rep(0, n), info = paste(info, "fg gr"))
    expect_equal(
      testfun$he(x),
      matrix(0, nrow = n, ncol = n),
      info = paste(info, "he")
    )
    expect_true(
      abs(testfun$fn(x) - testfun$fmin) > 1,
      info = paste(info, "stored fmin is not applicable")
    )
  }
})

test_that("current m metadata is recorded without making it a core field", {
  expected <- data.frame(
    factory = problem_factory_names(),
    m_present = c(
      FALSE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      FALSE,
      TRUE,
      TRUE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      FALSE,
      TRUE,
      TRUE,
      FALSE,
      FALSE,
      FALSE
    ),
    m_type = c(
      NA,
      "logical",
      NA,
      NA,
      NA,
      "logical",
      NA,
      NA,
      NA,
      "logical",
      NA,
      NA,
      NA,
      NA,
      "logical",
      NA,
      NA,
      NA,
      "logical",
      NA,
      "double",
      "double",
      "double",
      NA,
      NA,
      NA,
      NA,
      NA,
      "logical",
      NA,
      "logical",
      "double",
      NA,
      NA,
      NA
    ),
    m_value = c(
      NA,
      "NA",
      NA,
      NA,
      NA,
      "NA",
      NA,
      NA,
      NA,
      "NA",
      NA,
      NA,
      NA,
      NA,
      "NA",
      NA,
      NA,
      NA,
      "NA",
      NA,
      "100",
      "100",
      "100",
      NA,
      NA,
      NA,
      NA,
      NA,
      "NA",
      NA,
      "NA",
      "30",
      NA,
      NA,
      NA
    ),
    stringsAsFactors = FALSE
  )

  expect_equal(factory_m_metadata(), expected)
})
