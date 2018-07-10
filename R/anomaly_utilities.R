#' @section Anomaly Utilities:
#'
#' This section includes functions used in the identification of anomalies in groundwater quality data.
#'
#' This section includes
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' anomaly_statistical_methods
#' @description
#' Function to calculate the anomaly threshold using two statistical methods:
#' 1. mean + 2*std
#' 2. Quantile 90
#' @param gdata A geochemical dataset object
#' @param method A character string with the method to define the threshold
#' @return
#' This function returns a list with the following entries
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family anomaly functions
#' @importFrom stats sd
#' @importFrom stats quantile
#' @export
anomaly_statistical_methods <- function(gdata, method  = c("mean", "quantile")){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  #
  threshold <- NULL
  if(method == "mean"){
    mean.df <- apply(gdata$dataset, 2, mean)
    sd.df <- apply(gdata$dataset, 2, sd)
    threshold <- mean.df + 2*sd.df
  }
  else if(method == "quantile"){
    q90 <- apply(gdata$dataset, 2, quantile, probs = 0.9)
    threshold <- q90
  }
  else{
    stop("ERROR: unknown method")
  }
  res <- list(threshold = threshold)
}
#' @title
#' anomaly_graphical_method
#' @description
#' Function to create the QQplots required in the identification of homogeneous statistical
#' groups.
#' @param gdata A geochemical dataset object
#' @return
#' This function returns a ggplot2 object with the QQplots of each major ion
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family anomaly functions
#' @importFrom ggplot2 geom_qq
#' @importFrom ggplot2 ggplot
#' @export
anomaly_graphical_method <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  #
  dataset <- gdata$dataset
  Elements <- NULL
  p1 <- ggplot() + geom_qq(aes(sample = Elements), data = dataset)
  return(p1)
}
#' @title
#' anomaly_concentration_area_method
#' @description
#' Identification of anomalous areas using the concentration-area approach
#' @param gdata A geochemical dataset object
#' @param int.method A character string with the interpolation method. Currently three methods
#' are supported:
#' \itemize{
#' \item akima
#' \item idw
#' \item kriging
#' }
#' @return
#' This function returns a list with the following entries
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family anomaly functions
#' @export
anomaly_concentration_area_method <- function(gdata, int.method = c("akima",
                                                                    "idw",
                                                                    "kriging") ){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  #

}
#
anomaly_robust_kriging <- function(gdata){

}
#
anomaly_factorial_kriging <- function(gdata){

}
#
anomaly_polish_median_kriging <- function(gdata){

}
