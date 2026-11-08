test_that("low-dimensional Hessians match finite differences", {
  expect_hfd(broyden_tri(), broyden_tri()$x0(2), tolerance = 1e-5)
  expect_hfd(chebyquad(), chebyquad()$x0(1), tolerance = 1e-8)
  expect_hfd(disc_bv(), disc_bv()$x0(2), tolerance = 1e-5)
})
