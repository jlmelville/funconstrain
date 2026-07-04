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
#' - `fmin`: Reported minimum objective value.
#' - `xmin`: Numeric vector at a reported minimum.
#'
#' Some factories also include `m`, a metadata field for the number of summand
#' functions. It is absent for most factories, `NA` for several legacy fixed
#' problems, and numeric for `linfun_fr()`, `linfun_r1()`, `linfun_r1z()`, and
#' `trigon()`. It is not part of the core return contract needed by optimizers.
