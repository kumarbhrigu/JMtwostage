#' td: Identity Function (Internal)
#'
#' A simple identity function that returns the input as is. This function is typically used as a placeholder or for specific purposes like indicating time-dependent variables in model formulas.
#'
#' @param x A vector or object that is returned unchanged.
#'
#' @return The same value as the input.
#'
#' @details
#' The `td` function does not modify the input value and simply returns it. It can be used in situations where a time-dependent covariate is needed in model formulas, typically in survival analysis models.
#' @keywords internal
#' @export
td<-function(x)x


