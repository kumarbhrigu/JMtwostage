#' Landmark Time Cut for Longitudinal and Survival Data
#'
#' This function filters and processes survival and longitudinal data based on a specified landmark time.
#' It adjusts the time variables in the data (i.e., survival time and longitudinal visit time) relative to the landmark time.
#' The data are filtered to include only observations with survival times greater than or equal to the landmark time.
#'
#' @param ldata A data frame containing the longitudinal data. It must include the following columns:
#'   - `id`: A variable identifying subjects (e.g., subject ID).
#'   - `visitTime`: The time variable when the longitudinal measurement was recorded.
#'   - `stime`: The survival time variable.
#'   - `status`: The event status variable (e.g., 1 for event, 0 for censoring).
#'   - `timeDep`: A column containing time-dependent covariates.
#'   - `timeFixed`: A column containing time-fixed covariates.
#'
#' @param sdata A data frame containing the survival data. It must include the following columns:
#'   - `id`: A variable identifying subjects (same as in `ldata`).
#'   - `visitTime`: The time variable when the longitudinal measurement was recorded.
#'   - `stime`: The survival time variable.
#'   - `status`: The event status variable (same as in `ldata`).
#'   - `timeDep`: A column containing time-dependent covariates.
#'   - `timeFixed`: A column containing time-fixed covariates.
#'
#' @param id A character string specifying the name of the variable in both `ldata` and `sdata` that identifies subjects.
#'
#' @param stime A character string specifying the name of the survival time variable in both `ldata` and `sdata`.
#'
#' @param status A character string specifying the name of the event status variable in both `ldata` and `sdata`.
#'
#' @param visitTime A character string specifying the name of the visit time variable in both `ldata` and `sdata`.
#'
#' @param lmtime A numeric value specifying the landmark time. This is used to filter the data and adjust the time variables.
#'
#' @param timeDep A character string specifying the name of the time-dependent covariate(s) in both `ldata` and `sdata`.
#'
#' @param timeFixed A character string specifying the name of the time-fixed covariate(s) in both `ldata` and `sdata`.
#'
#' @return A list containing two components:
#'   - `Longdata`: A data frame containing the longitudinal data (`ldata`) filtered and adjusted relative to the landmark time.
#'   - `Survdata`: A data frame containing the survival data (`sdata`) filtered and adjusted relative to the landmark time.
#'
#' @details
#' The function checks that the landmark time is within the range of survival times in `sdata`. If the `lmtime` is out of bounds,
#' an error message is returned. Both the survival and longitudinal data are filtered to include only records where the survival
#' time is greater than or equal to the landmark time. After filtering, the function adjusts the survival and visit times by
#' subtracting the landmark time and ensures that negative longitudinal visit times are set to zero.
#'
#' @examples
#' # Example usage
#' result <- cutLMdata(
#'   ldata = long_data,
#'   sdata = surv_data,
#'   id = "ID",
#'   stime = "survival_time",
#'   status = "survival_status",
#'   visitTime = "Time",
#'   lmtime = 2,
#'   timeDep = c("marker_1","marker_2","marker_3"),
#'   timeFixed = c("Z_1","Z_2")
#' )
#' # Accessing the processed data:
#' head(result$Longdata)
#' head(result$Survdata)
#'
#' @seealso
#' `Surv` function for creating survival objects.
#'
#' @importFrom stats aggregate
#' @export
cutLMdata<-function(ldata,sdata,id,stime,status,visitTime,lmtime,timeDep,timeFixed){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  if(lmtime>max(sdata[,stime]))
    stop("landamrk time is greater the maximum survival time")
  if(lmtime<min(sdata[,stime]))
    stop("landmark time is less than the minimum survival time")
  #preparing data
  ldata<-ldata[,c(id,visitTime,stime,status,timeDep,timeFixed)]
  sdata<-sdata[,c(id,visitTime,stime,status,timeDep,timeFixed)]
  #cut survival data with respect to the landmark time
  sdata<-sdata[sdata[,stime]>=lmtime,]
  ldata<-ldata[ldata[,stime]>=lmtime,]
  sdata[,stime]<-sdata[,stime]-lmtime
  ldata[,stime]<-ldata[,stime]-lmtime
  ldata[,visitTime]<-ldata[,visitTime]-lmtime
  ldata[,visitTime]<-ifelse(ldata[,visitTime]<0,0,ldata[,visitTime])
  return(list(Longdata=ldata,Survdata=sdata))
}
