min_fx <- 0
min_fx8 <- 3.51687e-3 # for m = n = 8

min_x2 <- c(0.2113249, 0.7886751)
min_x3 <- c(0.1464466, 0.5000000, 0.8535534)
min_x8 <- c(
  0.04315276,
  0.19309084,
  0.26632871,
  0.50000000,
  0.50000000,
  0.73367129,
  0.80690916,
  0.95684724
)

chebyquad_objective <- function(par, m) {
  y <- 2 * par - 1
  residuals <- vapply(
    seq_len(m),
    function(i) {
      value <- mean(cos(i * acos(y)))
      if (i %% 2L == 0L) {
        value <- value + 1 / (i^2 - 1)
      }
      value
    },
    numeric(1L)
  )

  sum(residuals^2)
}

testfun <- chebyquad()
test_that("large Chebyquad residual counts do not overflow integer indices", {
  m <- 50000L
  problem <- chebyquad(m)
  x <- 0.5
  i <- as.double(seq_len(m))
  # At the shifted coordinate zero, T_i(0) cycles exactly through 0,-1,0,1.
  residual <- rep(c(0, -1, 0, 1), length.out = m)
  even <- i %% 2 == 0
  residual[even] <- residual[even] + 1 / (i[even]^2 - 1)
  expected_fn <- sum(residual^2)
  first <- rep(c(1, 0, -1, 0), length.out = m) * i
  second <- -i^2 * rep(c(0, -1, 0, 1), length.out = m)
  expected_he <- matrix(8 * sum(first^2 + residual * second))
  expect_warning(value <- problem$fn(x), NA)
  expect_warning(gradient <- problem$gr(x), NA)
  expect_warning(combined <- problem$fg(x), NA)
  expect_warning(hessian <- problem$he(x), NA)
  expect_true(all(is.finite(c(value, gradient, unlist(combined), hessian))))
  expect_equal(value, expected_fn)
  expect_equal(gradient, 0)
  expect_equal(combined$fn, value)
  expect_equal(combined$gr, gradient)
  expect_equal(hessian, expected_he)
})

test_that("Off-start derivatives match finite differences", {
  par <- c(0.2, 0.4, 0.6, 0.8)
  expect_gfd(testfun, par, tolerance = 1e-6)
  expect_hfd(testfun, par, tolerance = 1e-5)
  fg <- testfun$fg(par)
  expect_equal(fg$fn, testfun$fn(par))
  expect_equal(fg$gr, testfun$gr(par))
})
test_that("Explicit square m matches the dynamic default", {
  par <- c(0.12, 0.37, 0.64, 0.91)
  explicit <- chebyquad(m = length(par))

  expect_equal(explicit$fn(par), testfun$fn(par))
  expect_equal(explicit$gr(par), testfun$gr(par))
  expect_equal(explicit$he(par), testfun$he(par))
  expect_equal(explicit$fg(par), testfun$fg(par))
})
test_that("Overdetermined objectives and derivatives match independent oracles", {
  cases <- list(
    list(par = 0.31, m = 4L),
    list(par = c(0.12, 0.37, 0.64, 0.91), m = 7L),
    list(par = seq(0.08, 0.92, length.out = 8L), m = 12L)
  )

  for (case in cases) {
    overdetermined <- chebyquad(m = case$m)
    info <- paste0("n=", length(case$par), ", m=", case$m)

    expect_equal(
      overdetermined$fn(case$par),
      chebyquad_objective(case$par, case$m),
      tolerance = 1e-12,
      info = info
    )
    expect_gfd(overdetermined, case$par, tolerance = 1e-6)
    expect_hfd(overdetermined, case$par, tolerance = 1e-5)
    expect_true(isSymmetric(overdetermined$he(case$par)), info = info)

    fg <- overdetermined$fg(case$par)
    expect_equal(fg$fn, overdetermined$fn(case$par), info = info)
    expect_equal(fg$gr, overdetermined$gr(case$par), info = info)
  }
})
test_that("Explicit m is validated at factory and callback boundaries", {
  for (bad_m in list(0, -1, 1.5, Inf, NA_real_, TRUE, c(4, 5))) {
    expect_error(chebyquad(m = bad_m), "finite whole-number scalar|outside")
  }

  too_few_residuals <- chebyquad(m = 3)
  for (callback in c("fn", "gr", "he", "fg")) {
    expect_error(
      too_few_residuals[[callback]](rep(0.2, 4)),
      "m must be >= n",
      info = callback
    )
  }
  expect_error(too_few_residuals$x0(4), "m must be >= n")
})
test_that("Gradient is zero at stated minima", {
  expect_equal(testfun$gr(min_x2), rep(0, 2), tolerance = 1e-3)
  expect_equal(testfun$gr(min_x3), rep(0, 3), tolerance = 1e-3)
  expect_equal(testfun$gr(min_x8), rep(0, 8), tolerance = 1e-3)
})
test_that("Function value is correct at stated minima", {
  expect_equal(testfun$fn(min_x2), min_fx, tolerance = 1e-6)
  expect_equal(testfun$fn(min_x3), min_fx, tolerance = 1e-6)
  expect_equal(testfun$fn(min_x8), min_fx8, tolerance = 1e-5)
})
