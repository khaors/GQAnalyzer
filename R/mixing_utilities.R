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
#' @param gdata A geochemical_dataset object
#' @param end.members A geochemical_dataset with the chemical compositions of the
#' end-members of the solution.
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
#' @importFrom stats princomp lm residuals
#' @importFrom graphics plot
#' @importFrom MASS stepAIC
#' @importFrom grDevices chull
#' @references
#' Laaksoharju, M., Skarrman, C., \& Skarrman, E. (1999). Multivariate mixing and mass
#' balance (M3) calculations, a new tool for decoding hydrogeochemical information.
#' Applied Geochemistry, 14(7), 861–871. http://doi.org/10.1016/S0883-2927(99)00024-4
#' @export
m3_mixing_model <- function(gdata, end.members){
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset" |
     class(end.members) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  # Step 1. Multivariate Analysis
  # Scaling dataset = samples + end.members
  ions <- names(gdata$dataset)
  nvar <- ncol(gdata$dataset)
  samples.conc <- gdata$dataset[meql_ions]
  end.members.conc <- end.members$dataset[meql_ions]
  dataset <- rbind(samples.conc, end.members.conc)
  dataset <- scale(dataset, center = TRUE, scale = TRUE)
  nsamples <- nrow(dataset)
  n.endmembers <- nrow(end.members.conc)
  res.pca <- princomp(dataset, cor = TRUE, scores = TRUE)
  res.pca.def <- res.pca$scores
  # Step 2. Calculate mixing ratios using end members
  beginp <- nrow(samples.conc)+1
  endp <- nsamples
  end.members.pca <- res.pca.def[beginp:endp,]
  end.members.mat <- t(end.members.pca)
  end.members.pchull <- chull(end.members.pca[,1], end.members.pca[,2])
  end.members.chull <- end.members.pca[end.members.pchull,1:2]
  end.members.chull <- rbind(end.members.chull, end.members.chull[1,])
  end.members.chull.df <- as.data.frame(end.members.chull)
  end.members.chull.df$Type <- c(seq(1,n.endmembers, 1), 1)
  names(end.members.chull.df) <- c("PC1", "PC2", "Type")
  FFt <- t(end.members.mat)%*%end.members.mat
  Ie <- matrix(1.0, nrow = n.endmembers, ncol = 1)
  tmp1 <- cbind(FFt,Ie)
  A <- rbind(tmp1, c(t(Ie), 0))
  res.mixing.ratios <- matrix(0.0, nrow = nsamples, ncol = n.endmembers)
  for(i in 1:nsamples){
    b <- rbind(t(end.members.mat)%*%res.pca.def[i,], 1)
    mixing.ratios <- solve(A+1e-12*diag(n.endmembers+1),b)
    res.mixing.ratios[i,] <- mixing.ratios[1:n.endmembers]
  }
  # Step 3. Mass Balance Calculations
  X <- matrix(0.0, nrow = nsamples, ncol = nvar)
  for(isample in 1:nsamples){
    X[isample,] <- res.mixing.ratios[isample,]%*%t(end.members.mat)
  }
  #
  res.lm <- list()
  names.lm <- ions
  pc.names <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8')
  ion.residuals <- matrix(0.0, nrow = nsamples, ncol = nvar)
  rel.ion.residuals <- matrix(0.0, nrow = nsamples, ncol = nvar)
  for(i in 1:length(ions)){
    y <- dataset[1:nsamples,i]
    current.df <- data.frame(X[1:nsamples,])
    names(current.df) <- pc.names
    current.df$y <- y
    current.lm <- lm(y ~ ., data = current.df)
    current.slm <- stepAIC(current.lm, direction = "both", trace = 0)
    res.lm[[i]] <- current.slm
    ion.residuals[,i] <- residuals(current.slm)
    rel.ion.residuals[,i] <- 100*ion.residuals[,i]/y
  }
  res <- list(mixing.ratios = res.mixing.ratios, res.pca = res.pca.def,
              res.lm = res.lm, names.lm = names.lm, X = X, dataset = dataset,
              mass.balance = ion.residuals,
              relative.mass.balance = rel.ion.residuals,
              end.members = end.members.conc,
              original.samples = samples.conc,
              end.members.convexhull.df = end.members.chull.df)
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
#' carrera_mixing_model
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
carrera_mixing_model <- function(gdata.gd, end.members.gd){
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
