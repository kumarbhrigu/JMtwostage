#' Last Observation Carried Forward (LOCF) Imputation
#'
#' This function performs Last Observation Carried Forward (LOCF) imputation for missing values in a specified column (marker) of a dataset.
#' The missing values are replaced by the last observed value before them in the dataset.
#'
#' @param dataset A data frame or tibble that contains the data. The dataset should have at least one column with missing values to apply LOCF imputation.
#' @param marker The name (or index) of the column in the dataset where the missing values should be replaced with the last observed value.
#'
#' @return A data frame with the missing values in the `marker` column replaced by the last observed value before them.
#'
#' @examples
#' data1 <- data.frame(id = 1:5, value = c(1, NA, 3, NA, 5))
#' LOCF(data1, "value")
#' @export
LOCF<-function(dataset,marker)
{ if(nrow(dataset)>1)
{for(k in 2:(nrow(dataset)))
{ if(is.na(dataset[k,marker]))  dataset[k,marker]<-dataset[k-1,marker]
}}
  return(dataset)
}
