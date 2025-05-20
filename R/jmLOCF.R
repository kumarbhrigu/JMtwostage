#' Two stage Joint Model with LOCF(Last observation carry forward method)
#'
#' This function fits a joint model using the Last Observation Carried Forward (LOCF) method for time-dependent biomarkers with missing values. It handles longitudinal data and survival data jointly using a Cox proportional hazards model (or an Aalen additive model). The function imputes missing values for time-dependent biomarkers using the LOCF method.
#'
#' @param ldata A data frame containing the longitudinal data. This should include individual-level data with repeated measures over time.
#' @param sdata A data frame containing the survival data, including the time-to-event and event status variables.
#' @param id A string specifying the column name in both `ldata` and `sdata` identifying the individual subject (e.g., patient ID).
#' @param visitTime A string specifying the column name in `ldata` for the time points at which observations were made.
#' @param timeDep A vector of column names in `ldata` specifying the time-dependent biomarkers.
#' @param coxModel A formula representing the Cox proportional hazards model. It is required to specify the structure of the model with the time-dependent biomarker.
#' @param model A string specifying the type of survival model to use. The default is `"Cox"`, but `"Aalen"` can be specified for an Aalen additive model.
#'
#' @return A list containing the following components:
#' \item{res}{A summary of the fitted model.}
#' \item{coefs}{A data frame containing the coefficients from the fitted model.}
#' \item{vars}{A data frame containing the variance-covariance matrix of the coefficients.}
#' \item{nam}{The names of the coefficients.}
#' \item{fitted_model}{The fitted model object (either `coxph` or `aalen`).}
#' \item{imp.data}{A data frame containing the imputed longitudinal data.}
#' \item{complete_data}{A data frame containing the complete data with time-dependent biomarker imputations.}
#' \item{tname}{The name of the time variable used in the model.}
#' \item{ename}{The name of the event variable used in the model.}
#' \item{method}{The method used for imputations ("LOCF").}
#'
#' @details
#' This function implements a joint modeling framework for longitudinal and survival data using the LOCF method for imputing missing time-dependent biomarkers. It supports both the Cox proportional hazards model (`model="Cox"`) and the Aalen additive model (`model="Aalen"`).
#'
#' The longitudinal data is first pre-processed to create time-dependent covariates, and then the missing biomarker values are imputed using LOCF. The survival data is used to fit a joint model, where the survival model accounts for the time-dependent covariates.
#' @export
#' @references
#' Goodrich, B., et al. "rstanarm: Bayesian applied regression modeling via Stan. R package version 2.17. 4." Online< http://mc-stan. org (2018).
#'
#' Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
#' @examples
#'  library(survival)
#'  model_jmLOCF<-jmLOCF(ldata=long_data,sdata=surv_data,
#'        timeDep=c("marker_1","marker_2","marker_3"),
#'        visitTime="Time",
#'        coxModel=Surv(survival_time,survival_status)~Z_1+Z_2+td(marker_1)+
#'                      td(marker_2)+td(marker_3),
#'        model="Cox",
#'        id="ID")
#'  model_jmLOCF
#' @import survival
#' @importFrom utils tail
#' @importFrom dplyr group_by mutate
#' @importFrom stringr str_replace
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
jmLOCF <- function(ldata,
                    sdata,
                    id,
                    visitTime,
                    timeDep,
                    coxModel=NULL,
                    model="Cox") {
  call<-match.call()
  m<-match.call(expand.dots = FALSE)
  #preparing data
  tdlist<-timeDep
  if(sum(tdlist%in%names(ldata))!=length(timeDep)){
  stop("some biomarker is not included in the longitudinal data")
  }
  #imp_length<-length(impModel)
  #if(length(tdlist)!=length(impModel)){
  #stop("you have missed to put some imputation model")
  #}
  sVar<-all.vars(coxModel)
  special<-c("td")
  Terms<-terms(coxModel,special,data=sdata)
  m$formula<-Terms
  m[[1]]<-as.name("model.frame")

  #sdata, m = 5,id,visitTime,timeDep,impModel = NULL, coxModel = NULL
  m$ldata<-m$M<-m$visitTime<-m$timeDep<-m$coxModel<-m$model<-NULL
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
  # impVar<-list()
  # for(i in 1:imp_length){
  #   impVar[[i]]<-longmodeldata(data=ldata,longformula = impModel[[i]])
  # }
  # Recover names of all variables in model
  tname <- as.character(Terms[[2]][2])
  ename <- as.character(Terms[[2]][3])
  fnames <- attr(Terms, "term.labels")[-(attr(Terms, "specials")$td -1)]
  tdnames <- attr(Terms, "term.labels")[(attr(Terms, "specials")$td -1)]
  tdnames <- str_replace(tdnames, "td\\(", "")
  tdnames <- str_replace(tdnames, "\\)", "")
  ##############################################################################
  #prepare the data for Imputation
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



  impList<-list()
  for(i in 1:length(timeDep)){
    impList[[i]]<-LOCF(ldata,marker=tdnames[[i]])
  }


  comList<-ldata
  for(i in 1:length(timeDep)){
    comList[,c(timeDep[i])]<-impList[[i]][,tdnames[[i]]]
  }

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
    }
    else {
      fit<-suppressWarnings(aalen(as.formula(paste("Surv(time=tstart,time2=tstop,event=new_event2)~",
                                                     paste(paste("const(",c(tdnames,fnames),
                                                                 ")",sep=""),collapse="+"))),data=datcom,
                                    robust=F,n.sim =0,silent=0))
      fitted_model[[1]]<-fit
      coefs <- rbind(coefs, t(fit$gamma))
      vars <- rbind(vars, diag(fit$var.gamma))
      nam <- attr(fit$gamma, "dimnames")[[1]]
      nam <- str_replace(nam, "const\\(", "")
      nam <- str_replace(nam, "\\)", "")
    }
  res<-summ_cox(fit)

  result<-list()
  result$res<-res
  result$timeDep<-timeDep
  result$coefs<-coefs
  result$vars<-vars
  result$nam<-nam
  result$fitted_model<-fitted_model
  result$imp.data<-comList
  result$complete_data<-comList
  result$tname<-tname
  result$ename<-ename
  result$method<-"LOCF"
  class(result)<-"jmLOCF"
  return(result)
}
