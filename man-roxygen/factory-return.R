#' @return A list containing the core problem contract:
#'
#' - `fn`: Objective function. It takes a numeric parameter vector and returns
#'   the scalar objective value.
#' - `gr`: Gradient function. It takes a numeric parameter vector and returns
#'   the gradient vector.
#' - `he`: Hessian function. It takes a numeric parameter vector and returns the
#'   Hessian matrix of second derivatives.
#' - `fg`: Combined objective and gradient function. It takes a numeric
#'   parameter vector and returns a list with members `fn` and `gr`.
#' - `x0`: Suggested starting point. For fixed-dimension problems this is a
#'   numeric vector; for variable-dimension problems this is a function that
#'   returns a numeric vector for a requested `n`.
#' - `fmin`: A reported minimum objective value.
#' - `xmin`: A corresponding reported parameter vector, or an `NA` vector when
#'   no single minimizer is stored.
#'
#' For problems with variable `n` or configurable `m`, these stored references
#' may apply only to the configuration described in the factory's Minima
#' section; they are not recalculated for other choices of `n` or `m`. In some
#' cases, `fmin` applies more broadly than the stored `xmin`.
#'
#' Some factories also include `m`, a metadata field for the number of summand
#' functions. It is absent for most factories, `NA` for several legacy fixed
#' problems, and numeric for `linfun_fr()`, `linfun_r1()`, `linfun_r1z()`, and
#' `trigon()`. It is not part of the core return contract needed by optimizers.
