#' @section Mixing Utilities:
#'
#' This section is devoted to functions used in the study of mixing processes from hydrogeochenmical information.
#'
#' This section includes the definition of three functions used to create the Durov diagram:
#' m3_model
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' m3_mixing_model
#' @description
#' Function that implements the M# mixing model.
#' @param gdata A geochemical_dataset object
#' @param end.members Vector with the end members of the solution
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item res1
#' \item res2
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family mixing functions
#' @importFrom stats princomp lm
#' @importFrom graphics plot
#' @references
#' Laaksoharju, M., Skarrman, C., \& Skarrman, E. (1999). Multivariate mixing and mass
#' balance (M3) calculations, a new tool for decoding hydrogeochemical information.
#' Applied Geochemistry, 14(7), 861–871. http://doi.org/10.1016/S0883-2927(99)00024-4
#' @export
m3_mixing_model <- function(gdata, end.members){
  p <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  # Step 1. Multivariate Analysis
  # Scaling dataset
  dataset <- scale(gdata$dataset[meql_ions], center = TRUE, scale = TRUE)
  nsamples <- nrow(dataset)
  n.endmembers <- length(end.members)
  res.pca <- princomp(dataset, cor = TRUE, scores = TRUE)
  res.pca.def <- res.pca$scores[,1:2]
  # Step 2. Calculate mixing ratios using end members
  end.members.df <- gdata$dataset[end.members,]
  end.members.mat <- t(as.matrix(end.members.df))
  FFt <- t(end.members.mat)%*%end.members.mat
  Ie <- matrix(1.0, nrow = n.endmembers, ncol = 1)
  tmp1 <- cbind(FFt,Ie)
  A <- rbind(tmp1, c(t(Ie), 0))
  for(i in 1:nsamples){
    b <- rbind(t(end.members)%*%res.pca.def[i,], 1)
    mixing.ratios <- solve(A+1e-6*diag(n.endmembers),b)
    res.mixing.ratios[i,] <- mixing.ratios[1:5]
  }
  # Step 3. Mass Balance Calculations
  ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl","SO4")
  X <- matrix(0.0, nrow = nsamples, ncol = 2)
  for(i in 1:length(ions)){
    current.ion <- ions[i]
    y <- dataset[current.ion]
    for(isample in 1:nsamples){
      X[isample,] <- res.mixing.ratios[isample,]%*%end.members.mat
    }
  }
  current.df <- data.frame(y = y, PC1 = X[,1], PC2 = X[,2])
  res.lm <- lm(y ~ PC1+PC2, data = current.df)
  res <- list(mixing.ratios = res.mixing.ratios, res.pca = res.pca.def,
              res.lm = res.lm)
  return(res)
}
#' @title
#' select_end_members
#' @description
#' Function to estimate the end members of a solution
#' @param gdata A geochemical_dataset object
#' @return
#' This function returns
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family mixing functions
select_end_members <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  dataset <- scale(gdata$dataset, center = TRUE, scale = TRUE)
  res.pca <- princomp(dataset, cor = TRUE, scores = TRUE)
  plot(res.pca)
}
