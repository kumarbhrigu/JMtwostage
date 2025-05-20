#' Perform Multiple Imputation for a Longitudinal Response Variable Using Linear Mixed Effects Model
#'
#' This function performs multiple imputation for a longitudinal response variable using a linear mixed-effects model. It generates
#' imputed values for the missing response variable based on the model parameters estimated from the data.
#'
#' @param data1 A data frame containing the data for the longitudinal model.
#' @param impmodel A formula specifying the imputation model for the response variable, typically a linear mixed-effects model
#' (e.g., `y ~ x1 + x2 + (1|id)`).
#' @param nimp The number of imputations to perform.
#' @param idvar A string specifying the name of the variable identifying the subjects (ID variable).
#'
#' @return A list containing:
#' \item{imp.y.list}{A list of imputed response values for each imputation.}
#' \item{imp.data.list}{A list of data frames containing the imputed data for each imputation.}
#'
#' @details The function fits a linear mixed-effects model to the data, generates the model's parameters (fixed effects, random effects,
#' residual variance), and uses the model to impute missing values for the response variable. It incorporates the variability of the
#' model parameters across multiple imputations.
#' @import lme4
#' @import stats
#' @import insight
#' @import MASS
#' @keywords internal
#' @export
impY<-function(data1,impmodel,nimp,idvar){
  lme.imp.y<-lmer(impmodel,data=data1)
  #beta.mu.imp<-lme.imp.y@beta
  #beta.sigma.imp<-vcov(lme.imp.y)
  nj<-as.data.frame(table(data1[idvar]))$Freq
  m<-length(unique(data1[,idvar]))
  gamma.imp.y<-lme.imp.y@beta
  gamma.cov.y<-vcov(lme.imp.y)
  sigma.imp.y<-get_variance_residual(lme.imp.y)
  sigma.imp.b<-get_variance(lme.imp.y)$var.random

  gammalist<-mvrnorm(n=nimp,mu=gamma.imp.y,Sigma=gamma.cov.y)
  if(nimp==1){gammalist<-as.matrix(gammalist)}
  imp.var<-all.vars(nobars(impmodel))
  imp.y<-data1[,imp.var[1]]
  imp.x<- as.matrix(cbind(1,data1[imp.var[-1]]))
  r.y<-ifelse(is.na(data1[,imp.var[1]])==T,0,1)
  id.y<-rep(1:m,times=nj)
  nj.obs<-aggregate(r.y~id.y,FUN=sum,data=data.frame(r.y,id.y))$r.y
  obs.resid<-ifelse(r.y,(imp.y-imp.x%*%as.matrix(gamma.imp.y,ncol=1)),0)
  int.b<-(nj.obs/sigma.imp.y+1/sigma.imp.b)^(-1)*aggregate(obs.resid~id.y,FUN=sum,data=data.frame(obs.resid,id.y))$obs.resid/sigma.imp.y

  if(nimp!=1){
    imp.y.list<-list()
    imp.data.list<-list()
    for(i in 1:nimp){
      imp.y.list[[i]]<-c(imp.x%*%as.vector(gammalist[i,]))+int.b[id.y]

      imp.data.list[[i]]<-data1
      imp.y.list[[i]]<-ifelse(r.y==0,imp.y.list[[i]],data1[,imp.var[1]])

      imp.data.list[[i]][,imp.var[1]]<-ifelse(r.y==0,imp.y.list[[i]],data1[,imp.var[1]])
    }
  }else{
    imp.y.list<-c(imp.x%*%as.vector(gammalist))+int.b[id.y]

    imp.data.list<-data1
    imp.y.list<-ifelse(r.y==0,imp.y.list,data1[,imp.var[1]])
    imp.data.list[,imp.var[1]]<-ifelse(r.y==0,imp.y.list,data1[,imp.var[1]])
  }
  #int.nu<-Reduce('+',int.nulist)/M
  result<-list()
  result$imp.y.list<-imp.y.list
  result$imp.data.list<-imp.data.list
  return(result)
}
