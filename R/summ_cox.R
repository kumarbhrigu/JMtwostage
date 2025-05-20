#' Summarize Cox Proportional Hazards Model
#'
#' This function takes the output of a Cox proportional hazards model (`coxph`) and returns a summary data frame containing
#' the log hazard ratio (logHR), standard error (SE), confidence intervals (CI), and p-values for the model coefficients.
#'
#' @param fit A fitted Cox proportional hazards model object (of class `coxph`), typically returned by the `coxph` function.
#'
#' @return A data frame with the following columns:
#' \item{logHR}{The log hazard ratio for each coefficient.}
#' \item{SE}{The standard error of each coefficient.}
#' \item{CIlow}{The lower bound of the 95\% confidence interval for each coefficient.}
#' \item{CIupp}{The upper bound of the 95\% confidence interval for each coefficient.}
#' \item{p-value}{The p-value for testing the null hypothesis that each coefficient is zero.}
#'
#' @keywords internal
#' @export
summ_cox<-function(fit)
{
  if(length(summary(fit)$coef[,c("coef")])!=1)
    out<-data.frame(summary(fit)$coef[,c("coef")],summary(fit)$coef[,c("se(coef)")],
                    log(summary(fit)$conf.int[,3:4]),summary(fit)$coef[,5]) else
                      out<-data.frame(t(c(summary(fit)$coef[,c("coef")],summary(fit)$coef[,c("se(coef)")],
                                          log(summary(fit)$conf.int[,3:4]),summary(fit)$coef[,5])))

                    names(out)<-c("logHR","SE","CIlow","CIupp","p-value")
                    return(out)
}

