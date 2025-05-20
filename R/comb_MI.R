#' Combine Results from Multiple Imputation
#'
#' This function combines estimates and variances from multiple imputation using Rubin's rules. It calculates the combined estimate,
#' standard error, confidence intervals, and p-value based on the imputed coefficient estimates and variances.
#'
#' @param coefmi A numeric vector of coefficients from each imputation.
#' @param varmi A numeric vector of variances from each imputation.
#' @param M The number of imputations.
#'
#' @return A numeric vector containing:
#' \item{est}{The combined estimate (mean of the imputed coefficients).}
#' \item{se}{The combined standard error.}
#' \item{lowerCI}{The lower bound of the 95\% confidence interval.}
#' \item{upperCI}{The upper bound of the 95\% confidence interval.}
#' \item{pval}{The p-value for the hypothesis test of the combined estimate being zero.}
#' @importFrom stats qt pt
#' @keywords internal
#' @export
comb_MI<-function(coefmi,varmi,M)
{est<-mean(coefmi)
varinter<-(1/(M-1))*sum((coefmi-est)^2)
varintra<-mean(varmi)
se<-sqrt(varintra+(1+1/M)*varinter)
r<-(1+1/M)*varinter/varintra
vvv<-(M-1)*(1+1/r)^2
tal<-qt(0.025,df=vvv,lower.tail=F)
pval<-2*(1-pt(q=abs(est/se),df=vvv))
return(c(est,se,est-tal*se,est+tal*se,pval))}

