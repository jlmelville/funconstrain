validate_dimension <- function(
  value,
  problem,
  min = 1L,
  max = Inf,
  multiple = 1L,
  label = "n"
) {
  valid_scalar <- is.numeric(value) &&
    !is.logical(value) &&
    !is.complex(value) &&
    length(value) == 1L &&
    is.finite(value) &&
    value == floor(value)
  if (!valid_scalar) {
    stop(problem, ": ", label, " must be a finite whole-number scalar")
  }

  valid_multiple <- is.numeric(multiple) &&
    !is.logical(multiple) &&
    !is.complex(multiple) &&
    length(multiple) == 1L &&
    is.finite(multiple) &&
    multiple > 0 &&
    multiple == floor(multiple)
  if (!valid_multiple) {
    stop(problem, ": multiple must be a positive whole number")
  }

  if (value < min || value > max) {
    stop(problem, ": ", label, " is outside the allowed range")
  }
  if (value %% multiple != 0) {
    stop(problem, ": ", label, " must be a multiple of ", multiple)
  }

  value
}
