#' Summarize Results from Multiple Imputation for Survival Analysis
#'
#' This function combines results from multiple imputation for survival analysis, calculating the combined estimates,
#' standard errors, confidence intervals, and p-values for each coefficient using Rubin's rules. It applies the `combMI` function
#' for each column of coefficients and variances, and provides a summary output based on the specified model type (e.g., additive or hazard ratio).
#'
#' @param coefs A matrix or data frame where each column represents coefficients from an individual imputation.
#' @param vars A matrix or data frame where each column represents variances corresponding to the coefficients from each imputation.
#' @param model A string indicating the type of model. The options are `"Add"` for additive models, or any other string (e.g., `"Cox"`) for hazard ratio models.
#' @param nam A vector of names corresponding to each coefficient in the output.
#' @param M The number of imputations.
#'
#' @return A data frame with the following columns:
#' \item{HD}{The combined estimate (or hazard ratio) for each coefficient.}
#' \item{SE}{The combined standard error for each coefficient.}
#' \item{CIlow}{The lower bound of the 95\% confidence interval for each coefficient.}
#' \item{CIupp}{The upper bound of the 95\% confidence interval for each coefficient.}
#' \item{p-value}{The p-value for testing the null hypothesis that each coefficient is zero.}
#'
#'
#' @keywords internal
#' @export
summ_survtd<-function(coefs,vars,model,nam,M)
{
  out<-data.frame()
  for(j in 1:ncol(coefs))
    out<-rbind(out,comb_MI(coefs[,j],vars[,j],M))
  if(model=="Add")
    names(out)<-c("HD","SE","CIlow","CIupp","p-value")
  else
  {
    names(out)<-c("logHR","SE","CIlow","CIupp","p-value")
  }
  row.names(out)<-nam
  return(out)
}
