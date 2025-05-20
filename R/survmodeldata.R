#' Prepare Data for Survival Model
#'
#' This function prepares the data for fitting a survival model by extracting the relevant variables from a formula and returning a list
#' containing the variables for the survival analysis. It processes the survival formula, extracts the status and time variables, and
#' constructs the necessary data structure for further modeling.
#'
#' @param data A data frame containing the data for the survival model.
#' @param survformula A formula specifying the survival model, typically of the form `Surv(time, status) ~ predictors`.
#' @param idvar A string specifying the name of the variable identifying the subjects (ID variable). Default is `NULL`.
#' @param timeVar A string specifying the time variable to be used in the survival analysis. Default is `NULL`.
#'
#' @return A list containing:
#' \item{status}{A vector representing the event status variable from the survival data.}
#' \item{time}{A vector representing the time variable from the survival data.}
#' \item{time.var.long}{A vector representing the time variable for longitudinal data, if applicable.}
#' \item{x}{A matrix of covariates for the survival model (design matrix).}
#' \item{idvar}{The ID variable used to identify the subjects in the data.}
#' \item{idlist}{A vector of subject IDs from the data.}
#'
#' @details The function takes a survival formula (e.g., `Surv(time, status) ~ age + trt`) and extracts the status and time variables from
#' the `Surv()` function. It then extracts the predictor variables from the formula and prepares the necessary data for the survival model.
#'
#' @keywords internal
#' @export
survmodeldata<-function(data,survformula,idvar=NULL,timeVar=NULL){
  data<-data;formula<-survformula
  formula1<-deparse(survformula)
  survmatch<-gregexpr("Surv\\([^()]+\\)", formula1)
  surv.part<-regmatches(formula1,survmatch)
  surv.statustime<-strsplit(gsub("Surv\\(([^()]+)\\)", "\\1",surv.part ),",")[[1]]
  surv.status<-surv.statustime[1]
  surv.time<-surv.statustime[2]
  surv.x<-strsplit(formula1,"\\~")[[1]][2]
  framelist<-list()
  framelist$status<-data[,gsub("\\s+","",surv.status)]
  framelist$time<-data[,gsub("\\s+","",surv.time)]
  framelist$time.var.long<-data[,gsub("\\s+","",surv.time)]
  framelist$x<-model.matrix(as.formula(paste0("~",surv.x)),data)
  framelist$idvar<-idvar
  framelist$idlist<-data[idvar]
  return(framelist)
}
