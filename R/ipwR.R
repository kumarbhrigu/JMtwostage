#' Calculate Inverse Probability Weights (IPW) for Longitudinal Response Variable
#'
#' This function calculates inverse probability weights (IPW) for a longitudinal response variable based on a given imputation model and
#' IPW model. It uses a logistic regression model to estimate the probability of observing a non-missing value for the response variable,
#' and calculates the weights accordingly.
#'
#' @param data1 A data frame containing the data for the longitudinal model.
#' @param ipwmodel A formula specifying the model for the inverse probability weights, typically a logistic regression model (e.g., `~ x1 + x2`).
#' @param impmodel A formula specifying the imputation model for the response variable (e.g., `y ~ x1 + x2 + (1|id)`).
#' @param idvar A string specifying the name of the variable identifying the subjects (ID variable).
#'
#' @return A list containing:
#' \item{ipw.y}{A vector of the inverse probability weights for each observation.}
#' \item{fittedvalue.id}{A data frame with the fitted values for each subject ID.}
#'
#' @details The function first creates a binary variable indicating whether the response variable is missing or observed. It then fits a
#' logistic regression model for the probability of observing a non-missing response using the specified IPW model. The function calculates
#' the inverse probability weights as the inverse of the fitted values. Finally, it aggregates the weights at the subject level and returns
#' the calculated IPW and fitted values.
#'
#' @keywords internal
#' @export
ipwR<-function(data1,ipwmodel,impmodel,idvar){
  id<-data1$id
  ipwmodel<-deparse(ipwmodel)
  resp.y<-all.vars(nobars(impmodel))[1]
  ipwmodel<-as.formula(paste0("r.ipw.y",ipwmodel,collapse=""))
  ipw.var<-all.vars(nobars(ipwmodel))
  r.ipw.y<-ifelse(is.na(data1[,resp.y])==TRUE,0,1)
  data.r<-cbind(data1,r.ipw.y=r.ipw.y)
  lme.ipw.y<-glmer(ipwmodel,family=binomial(link='logit'),data=data.r)
  fitvalue<-fitted(lme.ipw.y)
  fitvalue.ipw<-aggregate(fitvalue~id,FUN=prod,data=data.frame(id,fitvalue))$fitvalue
  result<-list()
  result$ipw.y<-fitvalue.ipw
  result$fittedvalue.id<-data.frame(id=id,fittedvalue=fitvalue)
  return(result)
}
