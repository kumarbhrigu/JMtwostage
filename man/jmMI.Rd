% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jmMI.R
\name{jmMI}
\alias{jmMI}
\title{Two stage Joint Model with MI(Multiple imputation method)}
\usage{
jmMI(
  ldata,
  sdata,
  M = 5,
  id,
  visitTime,
  timeDep,
  impModel = NULL,
  coxModel = NULL,
  model = "Cox"
)
}
\arguments{
\item{ldata}{A data frame containing the longitudinal data. This should include variables for subject ID, time, and biomarkers.}

\item{sdata}{A data frame containing the survival data. This should include variables for the subject ID, survival time, and event indicator.}

\item{M}{The number of imputations to perform. Default is 5.}

\item{id}{A character string specifying the variable in `ldata` and `sdata` representing the subject ID.}

\item{visitTime}{A character string specifying the variable in `ldata` representing the visit times.}

\item{timeDep}{A character vector specifying the time-dependent covariates in the longitudinal data (`ldata`).}

\item{impModel}{A list of models for imputation of the time-dependent covariates. Each model should be a formula.}

\item{coxModel}{A formula specifying the Cox model for the survival data.}

\item{model}{A character string specifying the model to use for the survival analysis. Options are "Cox" (default) or "Aalen".}
}
\value{
A list containing the following components:
\item{res}{A summary of the results from the multiple imputation procedure.}
\item{coefs}{A matrix of the estimated coefficients from the survival model for each imputation.}
\item{vars}{A matrix of the variances of the estimated coefficients for each imputation.}
\item{nam}{A character vector of the names of the variables in the survival model.}
\item{fitted_model}{A list of the fitted survival models for each imputation.}
\item{imp.data}{A list of the imputed data sets.}
\item{miss.data}{The missing values in the time-dependent covariates.}
\item{complete_data}{A list of the complete imputed data sets.}
\item{tname}{The name of the survival time variable.}
\item{ename}{The name of the event indicator variable.}
\item{method}{A character string indicating the method used ("MI").}
}
\description{
This function implements multiple imputation for a joint model that handles time-dependent covariates with missing variables. It applies a Cox or Aalen model to the imputed data, combining results from multiple imputations using Rubin's rule.
}
\details{
This function is designed to handle joint modeling for longitudinal and survival data with time-dependent covariates. It allows for the imputation of missing data in the time-dependent covariates using the provided imputation models. The function can then apply a Cox or Aalen model to the imputed data and return the results.
}
\examples{
library(survival)
model_jmMI <- jmMI(ldata=long_data,sdata=surv_data,
              timeDep=c("marker_1","marker_2","marker_3"),
              impModel=list(marker_1~Z_1+Z_2+Time+(1|ID),
                           marker_2~Z_1+Z_2+Time+(1|ID),
                           marker_3~Z_1+Z_2+Time+(1|ID)),
              visitTime="Time",
              coxModel=Surv(survival_time,survival_status)~Z_1+Z_2+
              td(marker_1)+td(marker_2)+td(marker_3),
              model="Cox",
              id="ID")
model_jmMI
}
\references{
Goodrich, B., et al. "rstanarm: Bayesian applied regression modeling via Stan. R package version 2.17. 4." Online< http://mc-stan. org (2018).

Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
}
\author{
Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
}
