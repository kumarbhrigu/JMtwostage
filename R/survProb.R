#' Function for caluculating Baseline cummulative hazard
#'
#' This function calculates the baseline hazard, baseline cummulative hazard, at different unique event times
#'
#' @param fitted.object model object fitted through JMtwostage
#' @param nimp imputation number
#' @param conf_level confindence level
#'
#' @return list with baseline hazard, baseline cummulative hazard.
#' @export
#'
#' @examples
#' library(survival)
#' model_wMI<- jmwMI(ldata=long_data,sdata=surv_data,
#'                   timeDep=c("marker_1","marker_2","marker_3"),
#'                   impModel=list(marker_1~Z_1+Z_2+Time+(1|ID),
#'                                 marker_2~Z_1+Z_2+Time+(1|ID),
#'                                 marker_3~Z_1+Z_2+Time+(1|ID)),
#'                   ipwModel=list(~Z_1+Time+(1|ID),
#'                                 ~Z_1+Time+(1|ID),
#'                                 ~Z_1+Time+(1|ID)),
#'                   visitTime="Time",
#'                   coxModel=Surv(survival_time,survival_status)~Z_1+Z_2+
#'                            td(marker_1)+td(marker_2)+td(marker_3),
#'                   model="Cox",id="ID")
#' survProb(model_wMI,nimp=1)
survProb<- function(fitted.object,nimp=1,conf_level=0.95){
  nimp<-nimp
  n.fitted.object<-fitted.object$fitted_model[[nimp]]
  sbeta<-as.vector(n.fitted.object$coefficients)
  sbeta_se<-as.vector(sqrt(diag(vcov(n.fitted.object))))

  if(fitted.object$method=="MI"||fitted.object$method=="wMI"){
    event_time<-fitted.object$complete_data[[nimp]][,fitted.object$tname]
    event_status<-fitted.object$complete_data[[nimp]][,fitted.object$ename]
  }else{
    event_time<-fitted.object$complete_data[,fitted.object$tname]
    event_status<-fitted.object$complete_data[,fitted.object$ename]
  }

  X<-n.fitted.object$x
  # Function to estimate baseline hazard, cumulative hazard, survival, and survival CI
  # Step 1: Compute risk scores
  lin_pred <- X %*% sbeta  # Linear predictors
  risk_scores <- exp(lin_pred)  # Risk scores (exp(X beta))
  # Step 2: Sort data by event time
  sorted_indices <- order(event_time)
  event_time <- event_time[sorted_indices]
  event_status <- event_status[sorted_indices]
  risk_scores <- risk_scores[sorted_indices]
  # Step 3: Calculate baseline hazard and cumulative baseline hazard
  unique_times <- unique(event_time[event_status == 1])
  baseline_hazard <- numeric(length(unique_times))
  cum_baseline_hazard <- numeric(length(unique_times))

  for (i in seq_along(unique_times)) {
    time <- unique_times[i]
    at_risk <- event_time >= time
    events_at_time <- (event_time == time) & (event_status == 1)
    baseline_hazard[i] <- sum(events_at_time) / sum(risk_scores[at_risk])
    cum_baseline_hazard[i] <- sum(baseline_hazard[1:i])
  }

  # Return results as a list
  list(
    baseline_hazard = baseline_hazard,
    cum_baseline_hazard = cum_baseline_hazard,
    unique_times = unique_times
  )
}
