#' Prepare Data for Longitudinal Model with Random Effects
#'
#' This function prepares data for a longitudinal model by extracting the necessary components from the formula and splitting the data
#' according to the random effect terms. The function handles both fixed and random effects, including random slopes and intercepts,
#' and returns the processed data in a list suitable for model fitting.
#'
#' @param data A data frame containing the variables for the model.
#' @param longformula A formula specifying the longitudinal model with both fixed and random effects. The formula should include
#' the fixed effects on the right-hand side of the `~` and random effects enclosed in parentheses, for example:
#' `albumin ~ trt + age + day + (1|id) + (day|id)`.
#' @param idvar (Optional) A string specifying the name of the ID variable in the dataset. If not provided, the function assumes
#' that the ID variable is present in the formula.
#'
#' @return A list containing the following elements:
#' \item{id}{A vector with the ID variable from the dataset.}
#' \item{long.y}{A vector of the dependent variable (Y) for the longitudinal model.}
#' \item{long.x}{A matrix of the fixed effects (X) from the model formula.}
#' \item{long.z}{A list of model matrices for the random effect terms (Z).}
#' \item{varzfirst}{A list of the first variable names in the random effect terms.}
#' \item{varzsecond}{A list of the second variable names in the random effect terms.}
#'
#' @importFrom Matrix bdiag
#' @keywords internal
#' @export
longmodeldata<-function(data,longformula,idvar=NULL){
  #data<-PBC;formula<-albumin~trt+age+day+(1|id)+(day|id)
  formula1<-deparse(longformula)
  rdmatches<-gregexpr("\\([^()]+\\)", formula1)
  extracted_rdterms <- regmatches(formula1, rdmatches)
  fixedmodel<-gsub("\\s*\\+?\\s*\\([^()]+\\)", "", formula1)
  var.y<-strsplit(fixedmodel,"\\~")[[1]][1]
  var.x<-strsplit(fixedmodel,"\\~")[[1]][2]
  rdpattern <- "\\((.*?)\\|(.*?)\\)"
  var.zfirst<-list();var.zsecond<-list();var.z<-list()
  for(i in 1:length(extracted_rdterms[[1]])){
    match <- regexec(rdpattern,extracted_rdterms[[1]][i] )
    var.zfirst[[i]]<-regmatches(extracted_rdterms[[1]][i], match)[[1]][2]
    var.zsecond[[i]]<-regmatches(extracted_rdterms[[1]][i], match)[[1]][3]
    look.dim<-data.frame(table(data[,gsub("\\s+","",var.zsecond[[i]])]))$Freq
    zdatasplit<-split(data,data[,gsub("\\s+","",var.zsecond[[i]])])
    zdata.modelMat<-lapply(zdatasplit,function(x){model.matrix(as.formula(paste0("~",var.zfirst[[i]])),x)})
    var.z[[i]]<-bdiag(zdata.modelMat)
  }
  framelist<-list()
  framelist$id<-data[,idvar]
  framelist$long.y<-data[,gsub("\\s+","",var.y)]
  framelist$long.x<-model.matrix(as.formula(paste0("~",var.x)),data)
  framelist$long.z<-var.z
  framelist$varzfirst<-var.zfirst
  framelist$varzsecond<-var.zsecond
  return(framelist)
}
