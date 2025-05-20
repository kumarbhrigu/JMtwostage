#' Print method for jmLOCF objects
#'
#' @param x An object of class `jmLOCF`, returned by the `jmLOCF` function.
#' @param ... Additional arguments (not used).
#'
#' @return Printed summary of the `jmLOCF` object.
#' @rdname print.jmLOCF
#' @method print jmLOCF
#' @export
print.jmLOCF <- function(x, ...) {

  if(!inherits(x,'jmLOCF'))
    stop("\n Not a 'JMtwostate' object.\n")
  cat("Fitting two stage joint model :\n")

  # Print a brief description of the result
  cat("Method: ", x$method, "\n")
  cat("Time variable: ", x$tname, "\n")
  cat("Event variable: ", x$ename, "\n")


  cat("\nFitted model summary:\n")
  print(round(x$res,3))

  invisible(x)  # Return the object invisibly
}





#' Print method for jmMI objects
#'
#' @param x An object of class `jmMI`, returned by the `jmMI` function.
#' @param ... Additional arguments (not used).
#'
#' @return Printed summary of the `jmMI` object.
#' @rdname print.jmMI
#' @method print jmMI
#' @export
print.jmMI <- function(x, ...) {
  if(!inherits(x,'jmMI'))
    stop("\n Not a 'JMtwostate' object.\n")
  cat("Fitting two stage joint model :\n")

  # Print a brief description of the result
  cat("Method: ", x$method, "\n")
  cat("Time variable: ", x$tname, "\n")
  cat("Event variable: ", x$ename, "\n")



  # If the fitted model is available, print a summary
  cat("\nFitted model summary:\n")
  print(round(x$res,3))

  invisible(x)  # Return the object invisibly
}




#' Print method for jmMIIPW objects
#'
#' @param x An object of class `jmwMI`, returned by the `jmwMI` function.
#' @param ... Additional arguments (not used).
#'
#' @return Printed summary of the `jmwMI` object.
#' @rdname print.jmwMI
#' @method print jmwMI
#' @export
print.jmwMI <- function(x, ...) {
  if(!inherits(x,'jmwMI'))
    stop("\n Not a 'JMtwostate' object.\n")
  cat("Fitting two stage joint model :\n")

  # Print a brief description of the result
  cat("Method: ", x$method, "\n")
  cat("Time variable: ", x$tname, "\n")
  cat("Event variable: ", x$ename, "\n")

  # If the fitted model is available, print a summary
  cat("\nFitted model summary:\n")
  print(round(x$res,3))

  invisible(x)  # Return the object invisibly
}


#' Print method for jmCC objects
#'
#' @param x An object of class `jmCC`, returned by the `jmCC` function.
#' @param ... Additional arguments (not used).
#'
#' @return Printed summary of the `jmCC` object.
#' @rdname print.jmCC
#' @method print jmCC
#' @export
print.jmCC <- function(x, ...) {
  if(!inherits(x,'jmCC'))
    stop("\n Not a 'JMtwostate' object.\n")
  cat("Fitting two stage joint model :\n")

  # Print general information about the method
  cat("Method: ", x$method, "\n")
  cat("Time variable: ", x$tname, "\n")
  cat("Event variable: ", x$ename, "\n")

  # If the fitted model is available, print a summary
  cat("\nFitted model summary:\n")
  print(round(x$res,3))

  invisible(x)  # Return the object invisibly
}



#' Print method for jm2s objects
#'
#' @param x An object of class `jm2s`, returned by the `jm2s` function.
#' @param ... Additional arguments (not used).
#'
#' @return Printed summary of the `jm2s` object.
#' @rdname print.jm2s
#' @method print jm2s
#' @export
print.jm2s <- function(x, ...) {
  if(!inherits(x,'jm2s'))
    stop("\n Not a 'JMtwostate' object.\n")
  cat("Fitting two stage joint model :\n")

  # Print the method used
  cat("Method: ", x$method, "\n")

  # Print the time and event variables
  cat("Time variable: ", x$tname, "\n")
  cat("Event variable: ", x$ename, "\n")


  # Print the fitted model if available

    cat("\nFitted model summary:\n")
    print(round(x$res,3))


  invisible(x)  # Return the object invisibly
}
