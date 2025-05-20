#' Two stage Joint Model with 2S Method(single imputation)
#'
#' The `jm2s` function performs single imputation for time-dependent covariates (TDCs) in survival analysis. It imputes missing values for longitudinal biomarkers and fits a survival model (Cox or Aalen) using the imputed dataset. The function returns the results of the survival model fitting, including coefficients, variances, and imputed dataset.
#'
#' @param ldata A data frame containing the longitudinal data (biomarker measurements over time) with columns including the subject ID, visit times, and time-dependent covariates.
#' @param sdata A data frame containing the survival data, including the time-to-event information and event status.
#' @param id The name of the variable representing the unique individual ID in both `ldata` and `sdata`.
#' @param visitTime The name of the variable representing the visit time in the longitudinal data (`ldata`).
#' @param timeDep A character vector specifying the time-dependent variables (biomarkers) in `ldata` that will be imputed.
#' @param impModel A list of imputation models, one for each time-dependent variable in `timeDep`. Each model should be specified using a formula.
#' @param coxModel A survival model formula for Cox proportional hazards regression, such as `Surv(time, event) ~ x1 + x2 + td(x3)`.
#' @param model A string specifying the type of survival model to fit. Options are `'Cox'` (default) or `'Aalen'`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{res}{A summary of the results from the survival model.}
#'     \item{coefs}{A data frame containing the model coefficients.}
#'     \item{vars}{A data frame containing the variances of the model coefficients.}
#'     \item{nam}{A character vector containing the names of the model coefficients.}
#'     \item{fitted_model}{A list containing the fitted models for each imputation.}
#'     \item{imp.data}{A list of imputed data for the time-dependent variables.}
#'     \item{miss.data}{The missing data for the time-dependent variables in the original dataset.}
#'     \item{complete_data}{A list of completed datasets after imputation.}
#'     \item{tname}{The name of the time variable in the survival data.}
#'     \item{ename}{The name of the event variable in the survival data.}
#'     \item{method}{The method used for the analysis. Always "2s".}
#'   }
#' @export
#' @importFrom utils tail
#' @import survival
#' @references
#' Goodrich, B., et al. "rstanarm: Bayesian applied regression modeling via Stan. R package version 2.17. 4." Online< http://mc-stan. org (2018).
#'
#' Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
#' @examples
#' # Example usage:
#' library(survival)
#' model_jm2s<-jm2s(ldata=long_data,sdata=surv_data,
#'             timeDep=c("marker_1","marker_2","marker_3"),
#'             impModel=list(marker_1~Z_1+Z_2+Time+(1|ID),
#'                           marker_2~Z_1+Z_2+Time+(1|ID),
#'                           marker_3~Z_1+Z_2+Time+(1|ID)),
#'             visitTime="Time",
#'             coxModel=Surv(survival_time,survival_status)~Z_1+Z_2+
#'             td(marker_1)+td(marker_2)+td(marker_3),
#'             model="Cox",
#'             id="ID")
#' model_jm2s
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
jm2s <- function(ldata,
                 sdata,
                 id,
                 visitTime,
                 timeDep,
                 impModel=NULL,
                 coxModel=NULL,
                 model="Cox"){
  call<-match.call()
  m<-match.call(expand.dots=FALSE)
  #preparing data
  tdlist<-timeDep

  if(sum(tdlist%in%names(ldata))!=length(timeDep)){
    stop("some biomarker is not included in the longitudinal data")
  }
  imp_length<-length(impModel)
  if(length(tdlist)!=length(impModel)){
    stop("you have missed to put some imputation model")
  }
  sVar<-all.vars(coxModel)
  special<-c("td")
  Terms<-terms(coxModel,special,data=sdata)
  m$formula<-Terms
  m[[1]]<-as.name("model.frame")

  #sdata, m = 5,id,visitTime,timeDep,impModel = NULL, coxModel = NULL
  m$ldata<-m$visitTime<-m$impModel<-m$timeDep<-m$coxModel<-m$model<-NULL
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
  impVar<-list()
  for(i in 1:imp_length){
    impVar[[i]]<-longmodeldata(data=ldata,longformula = impModel[[i]])
  }
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
    impList[[i]]<-impY(data1=ldata,impmodel=impModel[[i]],nimp=1,idvar=id)
  }

  comList<-ldata
  for(j in 1:length(timeDep)){
     comList[,c(timeDep[j])]<-impList[[j]]$imp.y.list
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
  result$imp.data<-impList
  result$miss.data<-ldata[,tdlist]
  result$complete_data<-comList
  result$tname<-tname
  result$ename<-ename
  result$method<-"2s"
  class(result)<-"jm2s"
  return(result)
}

