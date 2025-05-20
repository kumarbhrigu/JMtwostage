#' Two stage Joint Model with CC(complete cases analysis method)
#'
#' This function performs joint modeling of survival and longitudinal data with time-dependent covariates using either a Cox proportional hazards model or an Aalen additive hazards model. The function handles the missing data in longitudinal measurements through complete case analysis(CC).
#'
#' @param ldata A data frame containing the longitudinal data (e.g., biomarkers over time).
#' @param sdata A data frame containing the survival data with time-to-event and event indicators.
#' @param id A character string indicating the column name for the subject identifier in both `ldata` and `sdata`.
#' @param visitTime A character string indicating the column name for the visit time in the longitudinal data.
#' @param timeDep A vector of time-dependent covariates that vary over time in the longitudinal data.
#' @param coxModel A `coxph` model object (optional). The formula specifies the time-to-event outcome, including time-dependent covariates.
#' @param model A character string specifying the model to use. Options are "Cox" for Cox proportional hazards or "Aalen" for Aalen additive hazards.
#'
#' @return A list containing the following elements:
#' \item{res}{A summary of the Cox model fit.}
#' \item{coefs}{A data frame of model coefficients for the fitted model.}
#' \item{vars}{A data frame of the variances of the coefficients.}
#' \item{nam}{The names of the variables in the model.}
#' \item{fitted_model}{A list containing the fitted model object.}
#' \item{imp.data}{The imputed data used for the analysis.}
#' \item{miss.data}{The missing values in the time-dependent covariates.}
#' \item{complete_data}{The complete data with imputed values.}
#' \item{tname}{The name of the time variable in the survival data.}
#' \item{ename}{The name of the event variable in the survival data.}
#' \item{method}{The method used for joint modeling (always "CC" for this function).}
#'
#' @details
#' This function requires the survival data (`sdata`) to contain a time-to-event variable and a binary event indicator. The longitudinal data (`ldata`) should contain time-varying biomarkers or covariates. The time-dependent covariates are specified via the `timeDep` argument, which must be present in `ldata`.
#'
#' The model formula for the Cox model is generated automatically based on the provided `coxModel`. If no `coxModel` is provided, the function will use a default Cox model.
#'
#' The time-dependent covariates are handled through imputation and the appropriate survival model (Cox or Aalen) is fit to the data.
#' @export
#' @references
#' Goodrich, B., et al. "rstanarm: Bayesian applied regression modeling via Stan. R package version 2.17. 4." Online< http://mc-stan. org (2018).
#'
#' Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
#' @examples
#' library(survival)
#' model_jmCC <-jmCC(ldata=long_data,sdata=surv_data,
#'               timeDep=c("marker_1","marker_2","marker_3"),
#'               visitTime="Time",
#'               coxModel=Surv(survival_time,survival_status)~Z_1+Z_2+
#'               td(marker_1)+td(marker_2)+td(marker_3),
#'               model="Cox",
#'               id="ID")
#' model_jmCC
#' @seealso \code{\link[survival]{coxph}}, \code{\link[timereg]{aalen}}
#' @importFrom utils tail
#' @importFrom survival coxph
#' @importFrom dplyr group_by mutate
#' @importFrom stringr str_replace
#' @importFrom stats model.frame
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
jmCC <- function(ldata,
                 sdata,
                 id,
                 visitTime,
                 timeDep,
                 coxModel=NULL,
                 model="Cox"){
  call<-match.call()
  m<-match.call(expand.dots=FALSE)
  #preparing data
  tdlist<-timeDep

  if(sum(tdlist%in%names(ldata))!=length(timeDep)){
    stop("some biomarker is not included in the longitudinal data")
  }

  sVar<-all.vars(coxModel)
  special<-c("td")
  Terms<-terms(coxModel,special,data=sdata)
  m$formula<-Terms
  m[[1]]<-as.name("model.frame")

  #sdata, m = 5,id,visitTime,timeDep,impModel = NULL, coxModel = NULL
  m$ldata<-m$visitTime<-m$timeDep<-m$coxModel<-m$model<-NULL
  m$id<-NULL
  m$sdata<-NULL
  m$na.action<-NULL
  m$data<-sdata
  m<-eval(m,parent.frame())
  mt<-attr(m,"terms")
  intercept<-attr(mt,"intercept")
  S<-model.extract(m,"response")
  if(!inherits(S,"Surv"))
    stop("Response must be a survival object")

  # Recover names of all variables in model
  tname <- as.character(Terms[[2]][2])
  ename <- as.character(Terms[[2]][3])
  fnames <- attr(Terms, "term.labels")[-(attr(Terms, "specials")$td -1)]
  tdnames <- attr(Terms, "term.labels")[(attr(Terms, "specials")$td -1)]
  tdnames <- str_replace(tdnames, "td\\(", "")
  tdnames <- str_replace(tdnames, "\\)", "")
  ##############################################################################
  ##############################################################################
  event_times<-sdata[,tname]
  new_sdata<-sdata[,c(id,tname,visitTime,ename)]
  #new_sdata_group<-group_by(new_sdata,by=id)
  #new_sdata_group<-mutate(new_sdata_group,tstart=head(visitTime),tstop=tail)
  st_function<-function(ldata,sdata,id,visitTime,stime,event){
    unique_id<-unique(ldata[,id])
    unique_stime<-sdata[,stime]
    unique_event<-sdata[,event]
    ldata_list<-list()
    for(i in 1:length(unique_id)){
      ldata_list[[i]]<-ldata[ldata[,id]==unique_id[i],]
      ldata_list[[i]]$tstart<-ldata_list[[i]][,visitTime]

      if(isTRUE(unique_stime[i]!=tail(ldata_list[[i]][,visitTime],1))){
        tstop<-c(tail(ldata_list[[i]][,visitTime],-1),unique_stime[i])}else{
          tstop<-c(tail(ldata_list[[i]][,visitTime],-1),unique_stime[i]+0.01)
        }
      new_event2<-c(rep(0,(nrow(ldata_list[[i]])-1)),unique_event[i])
      ldata_list[[i]]$tstop<-tstop
      ldata_list[[i]]$new_event2<-new_event2
    }
    ldata_list<-Reduce('rbind',ldata_list)
    return(ldata_list)
  }
  new_sdata<-st_function(ldata=ldata,sdata=new_sdata,id=id,visitTime=visitTime,stime=tname,event=ename)
  ldata<-cbind(ldata,tstart=new_sdata[,"tstart"],tstop=new_sdata[,"tstop"],new_event2=new_sdata[,"new_event2"])



  comList<-na.omit(ldata)

  ##############################################################################
  # use lme model on each imputed data
  # then impute time dependent coxph td() values
  # using the model parameters obtained from the last model
  ##############################################################################
  coefs <- data.frame()
  vars <- data.frame()
  fitted_model<-list()

  datcom <-comList
  #datcom <- datcom[datcom$ind == 2, ]
  if (model == "Cox") {
    fit <- suppressWarnings(coxph(as.formula(paste("Surv(time=tstart,time2=tstop,event=new_event2)~",
                                                   paste(c(tdnames, fnames), collapse = "+"))),
                                  data = datcom,x=TRUE))
    fitted_model[[1]]<-fit
    coefs <- rbind(coefs, summary(fit)$coef[, c("coef")])
    vars <- rbind(vars, (summary(fit)$coef[, c("se(coef)")])^2)
    nam <- attr(fit$coef, "names")
  } else {
    fit <- suppressWarnings(aalen(as.formula(paste("Surv(time=tstart,time2=tstop,event=new_event2)~",
                                                   paste(paste("const(", c(tdnames, fnames),
                                                               ")", sep = ""), collapse = "+"))), data = datcom,
                                  robust = F, n.sim = 0, silent = 0))
    fitted_model[[1]]<-fit
    coefs <- rbind(coefs, t(fit$gamma))
    vars <- rbind(vars, diag(fit$var.gamma))
    nam <- attr(fit$gamma, "dimnames")[[1]]
    nam <- str_replace(nam, "const\\(", "")
    nam <- str_replace(nam, "\\)", "")
  }

  res <- summ_cox(fit)
  result<-list()
  result$res<-res
  result$timeDep<-timeDep
  result$coefs<-coefs
  result$vars<-vars
  result$nam<-nam
  result$fitted_model<-fitted_model
  result$imp.data<-comList
  result$miss.data<-ldata[,tdlist]
  result$complete_data<-comList
  result$tname<-tname
  result$ename<-ename
  result$method<-"CC"
  class(result)<-"jmCC"
  return(result)
}

