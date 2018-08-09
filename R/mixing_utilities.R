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
#' Function that implements the M3 mixing model. This model is based on the projection of
#' the chemical information in a low-dimensional space defined using Principal Component
#' Analysis (PCA). The compositions of pre-defined end-members are also projected and from
#' these coordinates, the mixing ratios of each end-member are calculated using a constrained
#' least-squares procedure (mixing ratios add up to 1).
#' @param gdata A geochemical_dataset object with the major ions only.
#' @param end.members If provided, this is a numeric vector with the indices of the end members
#' to be used in the mixing calculations. If it is not specified then the end members are
#' defined using the convex hull of the first two principal components.
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item mixing.ratios: Matrix with the mixing ratios estimated using constrained
#' least-squares.
#' \item res.pca: Object returned by the princomp function. This object contains the
#' eigenvectors, eigenvalues, scores and other results of the PCA.
#' \item res.lm: List with the results of the multiple linear regression applied on each
#' ion included in the gdata geochemical_dataset object. Each entry of this list is a lm
#' object for each ion.
#' \item names.lm: Character vector with the names of the ions used in the mixing analysis.
#' \item X:
#' \item dataset:
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family mixing functions
#' @importFrom stats princomp lm residuals predict
#' @importFrom graphics plot
#' @importFrom MASS stepAIC
#' @importFrom grDevices chull
#' @references
#' Laaksoharju, M., Skarrman, C., \& Skarrman, E. (1999). Multivariate mixing and mass
#' balance (M3) calculations, a new tool for decoding hydrogeochemical information.
#' Applied Geochemistry, 14(7), 861–871. http://doi.org/10.1016/S0883-2927(99)00024-4
#' @export
m3_mixing_model <- function(gdata, end.members = NULL){
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
  res.pca <- princomp(dataset, cor = TRUE, scores = TRUE)
  res.pca.def <- res.pca$scores[,1:2]
  n.endmembers <- 0
  if(is.null(end.members)){
    pos.end.members <- chull(res.pca.def[,1], res.pca.def[,2])
    n.endmembers <- length(pos.end.members)
  }
  else{
    pos.end.members <- end.members
    n.endmembers <- length(end.members)
  }
  # Step 2. Calculate mixing ratios using end members
  end.members.mat <- t(res.pca.def[pos.end.members,])
  FFt <- t(end.members.mat)%*%end.members.mat
  Ie <- matrix(1.0, nrow = n.endmembers, ncol = 1)
  tmp1 <- cbind(FFt,Ie)
  A <- rbind(tmp1, c(t(Ie), 0))
  res.mixing.ratios <- matrix(0.0, nrow = nsamples, ncol = n.endmembers)
  for(i in 1:nsamples){
    b <- rbind(t(end.members.mat)%*%res.pca.def[i,], 1)
    mixing.ratios <- solve(A+1e-6*diag(n.endmembers+1),b)
    res.mixing.ratios[i,] <- mixing.ratios[1:n.endmembers]
  }
  # Correct mixing ratios for end members
  for(i in 1:n.endmembers){
    current.ratio <- rep(0.0, n.endmembers)
    current.ratio[i] <- 1
    res.mixing.ratios[pos.end.members[i],] <- current.ratio
  }
  #
  # Step 3. Mass Balance Calculations
  ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl","SO4")
  X <- matrix(0.0, nrow = nsamples, ncol = 2)
  residuals <- matrix(0.0, nrow = nsamples, ncol = length(ions))
  predicted <- matrix(0.0, nrow = nsamples, ncol = length(ions))
  lm.models <- list()
  for(i in 1:length(ions)){
    current.ion <- ions[i]
    y <- dataset[,i]
    for(isample in 1:nsamples){
      X[isample,] <- res.mixing.ratios[isample,]%*%t(end.members.mat)
    }
    current.df <- data.frame(y = y, PC1 = X[,1], PC2 = X[,2])
    res.lm <- lm(y ~ PC1 + PC2, data = current.df)
    lm.models[[i]] <- res.lm
    predicted[,i] <- predict(res.lm)
    residuals[,i] <- residuals(res.lm)
  }
  #
  res <- list(std.dataset = dataset,
              mixing.ratios = res.mixing.ratios,
              res.pca = res.pca.def,
              end.members = t(end.members.mat),
              predicted = predicted,
              residuals = residuals,
              lm.models = lm.models)
  return(res)
}
#
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
#' @title
#' uncertainty_m3_mixing_model
#' @description
#' Function to quantify the uncertainty in the mixing ratios due to the variability in
#' the concentrations of the end.members using Monte Carlo simulation.
#' @param gdata A geochemical_dataset object
#' @param end.members A geochemical_dataset object
#' @param lower A matrix with the lower limits of the concentrations of the end.members.
#' @param upper A matrix with the upper limits of the concentrations of the end.members.
#' @param nreal The number of realizations.
#' @param seed Random seed.
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item a
#' \item b
#' }
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family mixing functions
#' @importFrom stats rlnorm
#' @references
#' Laaksoharju, M., Skarrman, C., \& Skarrman, E. (1999). Multivariate mixing and mass
#' balance (M3) calculations, a new tool for decoding hydrogeochemical information.
#' Applied Geochemistry, 14(7), 861–871. http://doi.org/10.1016/S0883-2927(99)00024-4
uncertainty_m3_mixing_model <- function(gdata, end.members, lower, upper,
                                        nreal = 100, seed = 12345){
  set.seed(seed)
  if(class(gdata) != "geochemical_dataset" |
     class(end.members) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  #
  if(class(lower) != "matrix" | class(upper) != "matrix"){
    stop("ERROR: the lower and upper input variables must be numeric vectors")
  }
  #
  n.endmember <- nrow(lower)
  nvar <- ncol(lower)
  # Define alpha and beta parameters
  alpha <- (log(upper)+log(lower))/2
  beta <- (log(upper)-log(lower))/(2*2.576)
  # Define mu and sigma
  mu <- exp(alpha + (beta^2)/2)
  sigma <- mu*sqrt(exp(beta^2)-1)
  # Generate random end.members
  current.end.members <- array(0.0, dim = c(n.endmember, nvar, nreal))
  for(ireal in 1:nreal){
    for(iend in 1:n.endmember){
      for(ivar in 1:nvar){
        current.end.members[iend,ivar,] <- rlnorm(nreal, meanlog = mu[iend, ivar],
                                                  sdlog = sigma[iend, ivar])
      }
    }
  }
  # Estimate the mixing ratios
  return(current.end.members)
}
#' @title
#' plot_m3_mixing_results
#' @description
#' Function to plot different results obtained from the application of the M3 mixing model
#' @param gdata A geochemical_dataset object
#' @param mixing.res A list with ther results of the m3_mixing_model function
#' @param type A character string specifying the plot type to be created. Currently the
#' supported options are:
#' \itemize{
#' \item concentration
#' \item mixing.ratio
#' }
#' @param element A character string with the name of the ion
#' @return
#' This function returns a ggplot2 object with the requested plot.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family mixing functions
#' @importFrom ggplot2 ggplot geom_point geom_polygon ggtitle
plot_m3_mixing_results <- function(gdata, mixing.res,
                                   type = c("concentration",
                                            "mixing.ratio",
                                            "residual"),
                                   element){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  x <- NULL
  y <- NULL
  p <- NULL
  res.pca.df <- NULL
  current.title <- NULL
  ions <- list("Ca" = 1, "Mg" = 2, "Na" = 3, "K" = 4, "HCO3" = 5,
               "CO3" = 6, "Cl" = 7, "SO4" = 8)
  end.members.df <- data.frame(x = mixing.res$end.members[,1],
                               y = mixing.res$end.members[,2])
  #
  if(type == "concentration"){
    res.pca.df <- data.frame(PC1 = mixing.res$res.pca[,1],
                             PC2 = mixing.res$res.pca[,2],
                             concentration = mixing.res$residuals)
    #
    current.title <- paste0("Concentration: ", element)
    p <- ggplot() + geom_point(aes_string(x = "PC1", y = "PC2",
                                          color = paste0("residuals.",
                                                         ions[[element]])),
                               data = res.pca.df, size = 3) +
      scale_color_gradientn(colors=rainbow(10)) +
      geom_polygon(aes(x = x, y = y), data = end.members.df,
                   fill = NA, colour = "black") +
      theme_bw() +
      ggtitle(current.title)
  }
  else if(type == "mixing.ratio"){
    if(class(element) != "numeric"){
      stop("ERROR: the requested end member must be a number")
    }
    res.pca.df <- data.frame(PC1 = mixing.res$res.pca[,1],
                             PC2 = mixing.res$res.pca[,2],
                             mixing.ratio = mixing.res$residuals)
    #
    current.title <- paste0("Mixing Ratio: End Member", as.character(element))
    p <- ggplot() + geom_point(aes_string(x = "PC1", y = "PC2",
                                          color = paste0("mixing.ratio.",
                                                         as.character(element))),
                               data = res.pca.df, size = 3) +
      scale_color_gradientn(colors=rainbow(10)) +
      geom_polygon(aes(x = x, y = y), data = end.members.df,
                   fill = NA, colour = "black") +
      theme_bw() +
      ggtitle(current.title)
  }
  else if(type == "residual"){
    res.pca.df <- data.frame(PC1 = mixing.res$res.pca[,1],
                             PC2 = mixing.res$res.pca[,2],
                             residuals = mixing.res$residuals)
    #
    current.title <-
    p <- ggplot() + geom_point(aes_string(x = "PC1", y = "PC2",
                                   color = paste0("residuals.",
                                                  ions[[element]])),
                               data = res.pca.df, size = 3) +
      scale_color_gradientn(colors=rainbow(10)) +
      geom_polygon(aes(x = x, y = y), data = end.members.df,
                   fill = NA, colour = "black") +
      theme_bw() +
      ggtitle(element)
  }
  return(p)
}
#' @title
#' constrained_lm
#' @description
#' Constrained least-squares function. The constraint included ensures that the sum of the
#' estimated coefficients is equal to 1 and the coefficients are positive.
#' @param y A column matrix with the y values
#' @param X The design matrix of the problem
#' @param tol A numeric value with the tolerance specified to avoid numerical instability
#' during the solution of the least-squares problem.
#' @return
#' This function returns a numeric vector with the values of the estimated coeffcicients.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family mixing functions
constrained_lm <- function(y, X, tol = 1e-12){
  nx <- nrow(X)
  FFt <- t(X)%*%X
  Ie <- matrix(1.0, nrow = nx, ncol = 1)
  FFtconstrained <- cbind(FFt, Ie)
  Xmod <- rbind(FFtconstrained, c(t(Ie), 0))
  b <- rbind(t(X)%*%y, 1)
  res <- solve(Xmod+tol*diag(nx+1),b)
  return(res)
}
#' @title
#' mix_model
#' @description
#' Function to estimate mixing ratios with uncertain end.members
#' @param gdata.gd A geochemical_dataset
#' @param end.members.gd A geochemical_dataset
#' @return
#' A numeric vector
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family mixing functions
#' @importFrom pracma inv
mix_model <- function(gdata.gd, end.members.gd){
  if(class(gdata.gd) != "geochemical_dataset" |
     class(end.members.gd) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  dataset.mat <- as.matrix(gdata.gd$dataset)
  ndat <- nrow(dataset.mat)
  nvar <- ncol(dataset.mat)
  end.members.mat <- as.matrix(end.members.gd$dataset)
  nend <- nrow(end.members.mat)
  #1. Estimate initial mixing ratios using constrained lm
  mixing.ratios.mat <- matrix(0.0, nrow = ndat, ncol = nend)
  for(idat in 1:ndat){
    y <- matrix(dataset.mat[idat,], nrow = nvar, ncol = 1)
    mixing.ratios.mat[idat,] <- constrained_lm(y, end.members.mat, tol = 1e-12)
  }
  #2. Estimation of mu_s assuming mixing ratios known
  Gamma <- cbind(mixing.ratios.mat, -1*diag(ndat))
  Cs <- inv(Gamma%*%t(Gamma) + 1e-12*diag(ndat))
  zs <- dataset.mat
  lambdas <- -Cs%*%Gamma%*%zs
  sample.conc <-  zs - Gamma%*%lambdas
  #3. Estimate mixing ratios (Delta) assuming mu_s known

}
