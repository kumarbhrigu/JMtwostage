#' Plot Longitudinal and Survival Predictions from a Two stage Joint Model
#'
#' This function generates a plot combining longitudinal marker trajectories
#' and survival probability estimates from a fitted joint model.
#'
#' @param fitted.object A fitted object of JMtwosatge package.
#' @param newdata A data frame containing new covariate values for prediction. Defaults to `NULL`.
#' @param nimp An integer specifying which imputed dataset to use (for multiple imputation models). Defaults to `1`.
#' @param id A specific ID for which survival and marker trajectories should be plotted. Defaults to `NULL`.
#' @param marker Logical; if `TRUE`, includes longitudinal markers in the plot. Defaults to `TRUE`.
#'
#' @return A plot displaying the survival function along with longitudinal marker trajectories.
#' @details The function extracts survival estimates from the joint model
#' and overlays them with longitudinal marker values for an individual,
#' if `id` is provided. It supports models fitted using multiple imputation.
#'
#' @examples
#' library(survival)
#' model_wMI<- jmwMI(ldata=long_data,
#'                   sdata=surv_data,
#'                   timeDep=c("marker_1","marker_2","marker_3"),
#'                   impModel=list(marker_1~Z_1+Z_2+Time+(1|ID),
#'                                 marker_2~Z_1+Z_2+Time+(1|ID),
#'                                 marker_3~Z_1+Z_2+Time+(1|ID)),
#'                   ipwModel=list(~Z_1+Time+(1|ID),
#'                                 ~Z_1+Time+(1|ID),
#'                                 ~Z_1+Time+(1|ID)),
#'                   visitTime="Time",
#'                   coxModel=Surv(survival_time,survival_status)~Z_1+Z_2+
#'                   td(marker_1)+td(marker_2)+td(marker_3),
#'                   model="Cox",id="ID")
#' plot_jmtwostage(fitted.object = model_wMI,id=4)
#' @import survival
#' @import grDevices
#' @import graphics
#' @importFrom stats smooth.spline
#' @export
plot_jmtwostage<-function(fitted.object,newdata=NULL,nimp=1,id=NULL,marker=TRUE){
  m<-fitted.object
  if(fitted.object$method=="MI"||fitted.object$method=="wMI"){
    fitted.model<-m$fitted_model[[nimp]]
    com_data<-m$complete_data[[nimp]]
  } else {
    fitted.model<-m$fitted_model[[1]]
    com_data<-m$complete_data
  }

  if(is.null(id)==FALSE&&is.null(newdata)==TRUE){
    newdata<-com_data[com_data$ID==id,]
  }

  if(is.null(id)==TRUE&&is.null(newdata)==FALSE){
    newdata<-newdata
  }


  if(is.null(newdata)==FALSE){
    surv_fit<-survfit(fitted.model,newdata=as.data.frame(newdata)[nrow(newdata),])
  }else{
    surv_fit<-survfit(fitted.model)
  }

  surv_fit_surv<-surv_fit$surv
  surv_fit_up<-surv_fit$upper
  surv_fit_low<-surv_fit$lower
  surv_time<-surv_fit$time
  data_surv<-data.frame(Time=surv_time,Survival=surv_fit_surv,LL=surv_fit_low,
                        UL=surv_fit_up)

  smoothed_surv<-smooth.spline(data_surv$Time, data_surv$Survival, spar = 0.6)
  smoothed_LL<-smooth.spline(data_surv$Time, data_surv$LL, spar = 0.7)
  smoothed_UL<-smooth.spline(data_surv$Time, data_surv$UL, spar = 0.7)

  if(marker==TRUE&is.null(id)==FALSE){
    if(nrow(newdata)>1){
      marker_data<-newdata[,c("Time",fitted.object$timeDep)]
    }else{
      marker_data<-NULL
    }
  }

  if(is.null(marker_data)==TRUE){
    warning("Required number of marker data is not available for the longitudinal plot. Still here is the survival plot.")
  }


  if(is.null(marker_data)==FALSE){
    oldpar1 <- par(no.readonly = TRUE)
    K<-ifelse(length(fitted.object$timeDep)>3,3,length(fitted.object$timeDep))
    m<-cbind(1:K,rep(K+1,K))
    widths<-c(0.3,0.35)
    layout(m,widths=widths)
    par(mar=c(0,5,0,0),oma=c(4,0,3,0),cex.axis=1.1,
        font.axis=2,font.lab=2,font.main=2)
    plot(x=marker_data$Time,y=marker_data$marker_1,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
    points(x=marker_data$Time,y=marker_data$marker_1,pch=8,col='red')
    box(lwd=1.3)
    mtext("marker_1",line=2.5,side=2,font=2)

    par(mar = c(0, 5, 0, 0))
    plot(x=marker_data$Time,y=marker_data$marker_2,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
    points(x=marker_data$Time,y=marker_data$marker_2,pch=8,col='red')
    #axis(1, at = xticks)
    box(lwd=1.3)
    mtext("marker_2",line=2.5,side=2,font=2)
    #xticks<-pretty(c(pred_list[[3]]$Time))

    par(mar = c(0,5, 0, 0))
    plot(x=marker_data$Time,y=marker_data$marker_3,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
    points(x=marker_data$Time,y=marker_data$marker_3,pch=8,col='red')
    #axis(1, at = xticks)
    box(lwd=1.3)
    mtext("marker_3",line=2.5,side=2,font=2)
    axis(1,las=1)

    par(mar=c(0,0,0,5))
    plot(x=smoothed_surv$x,y=smoothed_surv$y,col='black',type='l',las=1,
         yaxt='n',xaxs='i',xaxt='n',ylab='',ylim=c(0,1))

    # Add confidence interval shading
    polygon(c(smoothed_LL$x, rev(smoothed_LL$x)),
            c(smoothed_UL$y, rev(smoothed_LL$y)),
            col = rgb(0, 0, 0.6, 0.4),border = NA)
    lines(x=smoothed_surv$x,y=smoothed_surv$y,col='blue',lwd=2)
    # Add a bold boundary line to the plot
    axis(1,las=1)
    axis(4,las=1)

    # Customize axis text and titles in bold
    mtext(paste("Longitudinal and Survival trajectories for ID=",id), 3,
          line = 1, outer = TRUE, font = 2, cex = 1.3)
    mtext("Time", 1,
          line = 2.5, outer = TRUE,font=2)
    mtext("Survival Prediction", 4,
          line = 2.5,font=2)
    box(lwd=1.3)
    on.exit(par(oldpar1))

  }else{
    oldpar2 <- par(no.readonly = TRUE)
    plot(NA, xlim=range(smoothed_surv$x), ylim=c(0,1), xlab='', ylab='', type='n')

    # Add confidence interval polygon
    polygon(c(smoothed_LL$x, rev(smoothed_LL$x)),
            c(smoothed_UL$y, rev(smoothed_LL$y)),
            col=rgb(0, 0, 0.6, 0.4), border=NA)

    # Add smoothed survival trajectory
    lines(smoothed_surv$x, smoothed_surv$y, col='blue', lwd=2)

    # Add axis labels and title
    mtext("Time", side=1, line=2, font=2)
    mtext("Survival Probability", side=2, line=2, font=2)
    mtext(paste("Survival trajectory for ID:", id), side=3, line=1.5, font=2, cex=1.3)
    box(lwd=1.3)
    on.exit(par(oldpar2))

  }

}
