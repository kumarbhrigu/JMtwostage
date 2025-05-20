#' Summarize Results from Aalen's Additive Regression Model
#'
#' This function extracts and summarizes the results from an Aalen's additive regression model fit. It calculates the estimates,
#' standard errors, 95% confidence intervals, and p-values for each coefficient, based on the model fit object.
#'
#' @param fit A fitted model object from Aalen's additive regression, typically obtained from the `aalen()` function.
#'
#' @return A data frame containing the following columns:
#' \item{HD}{The estimated hazard coefficient for each variable.}
#' \item{SE}{The standard error of the hazard coefficient.}
#' \item{CIlow}{The lower bound of the 95\% confidence interval for the hazard coefficient.}
#' \item{CIupp}{The upper bound of the 95\% confidence interval for the hazard coefficient.}
#' \item{p-value}{The p-value for testing the null hypothesis that the hazard coefficient is zero.}
#'
#' @keywords internal
#' @export
summ_aalen<-function(fit)
{
  Estimate<-t(fit$gamma)
  SE<-sqrt(diag(fit$var.gamma))
  CILow<-Estimate-1.96*SE
  CIUpp<-Estimate+1.96*SE
  Pvalue<-2*(1-pnorm(abs(Estimate/SE)))
  out<-data.frame(t(Estimate),SE,t(CILow),t(CIUpp),t(Pvalue))
  names(out)<-c("HD","SE","CIlow","CIupp","p-value")
  nam<-row.names(out)
  nam<-str_replace(nam,"const\\(","")
  nam<-str_replace(nam,"\\)","")
  row.names(out)<-nam
  return(out)
}
