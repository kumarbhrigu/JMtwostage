#' Summary method for jmLOCF objects
#'
#' @param object An object of class `jmLOCF`, returned by the `jmLOCF` function.
#' @param ... Additional arguments (not used).
#'
#' @return A summary of the `jmLOCF` object.
#' @rdname summary.jmLOCF
#' @method summary jmLOCF
#' @export
summary.jmLOCF <- function(object, ...) {
  if(!inherits(object,'jmLOCF'))
    stop("\n Not a 'JMtwostate' object.\n")
  result <- list()

  # General information about the method and variables
  result$method <- object$method
  result$tname <- object$tname
  result$ename <- object$ename

  # Coefficients and variance
  result$coefs <- object$coefs
  result$vars <- object$vars

  # Fitted model summary
  if (!is.null(object$fitted_model)) {
    result$fitted_model_summary <- summary(object$fitted_model[[1]])
  } else {
    result$fitted_model_summary <- NULL
  }

  # Return the result as a summary list
  return(result)
}



#' Summary method for jmMI objects
#'
#' @param object An object of class `jmMI`, returned by the `jmMI` function.
#' @param ... Additional arguments (not used).
#'
#' @return A summary of the `jmMI` object.
#' @rdname summary.jmMI
#' @method summary jmMI
#' @export
summary.jmMI <- function(object, ...) {
  if(!inherits(object,'jmMI'))
    stop("\n Not a 'JMtwostate' object.\n")
  # Prepare a summary
  result <- list()

  # General information about the method and variables
  result$method <- object$method
  result$tname <- object$tname
  result$ename <- object$ename

  # Coefficients and variance
  result$coefs <- object$coefs
  result$vars <- object$vars

  # Fitted model summary
  if (!is.null(object$fitted_model)) {
    result$fitted_model_summary <- summary(object$fitted_model[[1]])
  } else {
    result$fitted_model_summary <- NULL
  }

  # Return the result as a summary list
  return(result)
}



#' Summary method for jmwMI objects
#'
#' @param object An object of class `jmwMIIPW`, returned by the `jmwMI` function.
#' @param ... Additional arguments (not used).
#'
#' @return A summary of the `jmwMI` object.
#' @rdname summary.jmwMI
#' @method summary jmwMI
#' @export
summary.jmwMI <- function(object, ...) {
  if(!inherits(object,'jmwMI'))
    stop("\n Not a 'JMtwostate' object.\n")
  # Prepare a summary
  result <- list()

  # General information about the method and variables
  result$method <- object$method
  result$tname <- object$tname
  result$ename <- object$ename

  # Coefficients and variance
  result$coefs <- object$coefs
  result$vars <- object$vars

  # Fitted model summary
  if (!is.null(object$fitted_model)) {
    result$fitted_model_summary <- summary(object$fitted_model[[1]])
  } else {
    result$fitted_model_summary <- NULL
  }

  # Return the result as a summary list
  return(result)
}


#' Summary method for jmCC objects
#'
#' @param object An object of class `jmCC`, returned by the `jmCC` function.
#' @param ... Additional arguments (not used).
#'
#' @return A summary of the `jmCC` object.
#' @rdname summary.jmCC
#' @method summary jmCC
#' @export
summary.jmCC <- function(object, ...) {
  if(!inherits(object,'jmCC'))
    stop("\n Not a 'JMtwostate' object.\n")
  # Prepare a summary
  result <- list()

  # General information about the method and variables
  result$method <- object$method
  result$tname <- object$tname
  result$ename <- object$ename

  # Coefficients and variance
  result$coefs <- object$coefs
  result$vars <- object$vars

  # Fitted model summary
  if (!is.null(object$fitted_model)) {
    result$fitted_model_summary <- summary(object$fitted_model[[1]])
  } else {
    result$fitted_model_summary <- NULL
  }

  # Return the result as a summary list
  return(result)
}


#' Summary method for jm2s objects
#'
#' @param object An object of class `jm2s`, returned by the `jm2s` function.
#' @param ... Additional arguments (not used).
#'
#' @return A summary of the `jm2s` object, including method, coefficients, variance, and fitted model.
#' @rdname summary.jm2s
#' @method summary jm2s
#' @export
summary.jm2s <- function(object, ...) {
  if(!inherits(object,'jm2s'))
    stop("\n Not a 'JMtwostate' object.\n")
  # Prepare a summary object
  result <- list()

  # General information about the method and variables
  result$method <- object$method
  result$tname <- object$tname
  result$ename <- object$ename

  # Coefficients and variance of the model
  result$coefs <- object$coefs
  result$vars <- object$vars

  # Fitted model summary if available
  if (!is.null(object$fitted_model)) {
    result$fitted_model_summary <- summary(object$fitted_model[[1]])
  } else {
    result$fitted_model_summary <- NULL
  }

  # Return the detailed summary
  return(result)
}
