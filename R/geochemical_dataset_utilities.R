#' GQAnalyzer: A package for Hydrogeochemical plotting
#'
#' This package provides functions to create the following plots:
#' \itemize{
#' \item Piper plot
#' \item Stiff diagrams
#' \item Durov plot
#' \item Ternary
#' \item Radial
#' \item Multirectangular
#' \item Maucha
#' }
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' geochemical_dataset
#' @description
#' Function to create a geochemical_dataset.
#' @param name Name of the geochemical dataset
#' @param data Data.frame with the concentration of anions and cations (see details)
#' @return
#' An object of the geochemical_dataset class with the following member variables:
#' \itemize{
#' \item name: Name of the geochemical dataset
#' \item dataset: data.frame with the corresponding measured concentrations
#' \item meqL: data.frame with the meqL of the measured concentrations
#' \item anions: data.frame with the meqL of the major anions
#' \item cations: data.frame with the meqL of the major cations
#' \item meqL.min: Numeric vector with the minimum values of the meqL
#' \item meqL.max: Numeric vector with the maximum values of the meqL
#' \item meqL.schoeller: data.frame with the concentrations of the major ions organized to
#' create tehe Schoeller plot.
#' }
#' @details
#' A geochemical dataset is a data.frame with the values of the concentrations of the anions
#' and cations obtained during hydrogeochemical sampling. Internally this object calculates the
#' meqL of the measured concentrations of the major ions. The concentrations of the ions must be
#' specified in mg/L.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
geochemical_dataset <- function(name = "GeochemicalDataset", data){
  if(class(data) != "data.frame"){
    stop('ERROR: A data.frame is required as input')
  }
  ndat <- nrow(data)
  nvar <- ncol(data)
  major_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  durov_var <- c("pH", "TDS")
  dnames <- colnames(data)
  print(dnames)
  for(i in 1:8){
    if(!major_ions[i] %in% dnames){
      msg <- paste0('ERROR: a major ion ', major_ions[i], ' is not present in the input data.frame')
      stop(msg)
    }
  }
  # Check for durov variables
  for(i in 1:2){
    if(!durov_var[i] %in% dnames){
      warning("Warning: Durov plot can not be created with current data.frame")
    }
  }
  # Convert concentrations into meqL and obtain min and max of results
  tmp <- convert_meql(gdata = data)
  colnames(tmp$meqL) <- major_ions
  meqL.mn <- apply(tmp$meqL, 2, min)
  meqL.mx <- apply(tmp$meqL, 2, max)
  meqL.df <- as.data.frame(tmp$meqL+1e-3)
  #names(meqL.df[,9]) <- "ID"
  #print(meqL.df)
  #
  meqL.df1 <- meqL.df %>% gather(key = "major_ions", value = "concentration")
  # Calculate the fractions of conservative mixing between fresh and sea water.
  tmp1 <- calculate_convervative_mixing(gdata = data)
  mix.df <- as.data.frame(tmp1$mix)
  names(mix.df) <- major_ions
  #
  result <- list(name = name,
                 dataset = data,
                 meqL = tmp$meqL,
                 anions = tmp$anions,
                 cations = tmp$cations,
                 meqL.min = meqL.mn,
                 meqL.max = meqL.mx,
                 meqL.schoeller = meqL.df1,
                 sea.water.fraction = tmp1$f.sea,
                 fresh.water.fraction = 1.0-tmp1$f.sea,
                 mix.mmol = tmp1$mix.mmol,
                 mix = mix.df)
  class(result) <- "geochemical_dataset"
  invisible(result)
  return(result)
}
#' @title
#' [.geochemical_dataset
#' @description
#' Function to access an element of a geochemical_dataset.
#' @param x A geochemical_dataset
#' @param i An index
#' @param ... Additional parameters (for consistency purposes).
#' @return
#' This function returns a specific entry of a geochemical_dataset as a geochemical_dataset
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom dplyr full_join
#' @export
`[.geochemical_dataset` <- function(x, i, ...){
  if(class(x) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  #
  ndat <- nrow(x$dataset)
  nvar <- ncol(x$dataset)
  #
  if(i < 0 || i> ndat){
    stop('ERROR: the requested sample is out of limits')
  }
  #
  IDsamples <- NULL
  #IDsamples <- x$ID[i,]
  nsamples <- length(i)
  tmp <- data.frame(ions = character(), concentration = double())
  #
  res <- list(dataset = x$dataset[i,],
              meqL = matrix(x$meqL[i,], nrow = nsamples),
              anions = matrix(x$anions[i,], nrow = nsamples),
              cations = matrix(x$cations[i,], nrow = nsamples),
              meqL.min = x$meqL.min,
              meqL.max = x$meqL.max) #,
              #meqL.schoeller = tmp)
  class(res) <- "geochemical_dataset"
  invisible(res)
  return(res)
}
#' @title
#' convert_meql
#' @description
#' Function to calculate the miliequivalent of a given composition.
#' @param gdata A data.frame with the concentrations of the major ions
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item anions: A matrix with the values of the normalized meqL of the anions.
#' \item cations: A matrix with the values of the normalized meqL of the cations.
#' \item meqL: A matrix with the values of the normalized meqL of the meqL.
#' \item mmol: A matrix with the concentration values in mmol/L.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export convert_meql
convert_meql <- function(gdata){
  #Ca, Mg, Na, K, HCO3, CO3, Cl, SO4
  gmol <- c(40.078, 24.305, 22.989768, 39.0983, 61.01714, 60.0092, 35.4527, 96.0636)
  eqmol <- c(2., 2., 1., 1., 1., 2., 1., 2.)
  ndat <- length(gdata$Ca)
  meqL <- matrix(0.0, nrow = ndat, ncol = 8)#HCO3	CO3	Cl	SO4	F
  meqL[,1] <- gdata$Ca*eqmol[1]/gmol[1]
  meqL[,2] <- gdata$Mg*eqmol[2]/gmol[2]
  meqL[,3] <- gdata$Na*eqmol[3]/gmol[3]
  meqL[,4] <- gdata$K*eqmol[4]/gmol[4]
  meqL[,5] <- gdata$HCO3*eqmol[5]/gmol[5]
  meqL[,6] <- gdata$CO3*eqmol[6]/gmol[6]
  meqL[,7] <- gdata$Cl*eqmol[7]/gmol[7]
  meqL[,8] <- gdata$SO4*eqmol[8]/gmol[8]
  sumcations <- rowSums(meqL[,1:4])
  sumanions <- rowSums(meqL[,5:8])
  # Coordinates of cations
  cations <- matrix(0.0, nrow = ndat, ncol = 3)
  cations[,1] <- meqL[,1]/sumcations
  cations[,2] <- meqL[,2]/sumcations
  cations[,3] <- (meqL[,3]+meqL[,4])/sumcations
  # Coordinates of anions
  anions <- matrix(0.0, nrow = ndat, ncol = 3)
  anions[,1] <- (meqL[,5]+meqL[,6])/sumanions
  anions[,3] <- meqL[,7]/sumanions
  anions[,2] <- meqL[,8]/sumanions
  #
  mmol <- matrix(0.0, nrow = ndat, ncol = 8)
  mmol[,1] <- gdata$Ca/gmol[1]
  mmol[,2] <- gdata$Mg/gmol[2]
  mmol[,3] <- gdata$Na/gmol[3]
  mmol[,4] <- gdata$K/gmol[4]
  mmol[,5] <- gdata$HCO3/gmol[5]
  mmol[,6] <- gdata$CO3/gmol[6]
  mmol[,7] <- gdata$Cl/gmol[7]
  mmol[,8] <- gdata$SO4/gmol[8]
  #
  results <- list(anions = anions, cations = cations, meqL = meqL, mmol = mmol)
  return(results)
}
#' @title
#' calculate_conservative_mixing
#' @description
#' Function to calculate the mixing ratio of a perfect-conservative mixture between fresh and
#' sea water from a given geochemical dataset.
#' @param gdata A data.frame with the concentrations of the major ions
#' @param fresh.water The chemical composition of fresh water
#' @param sea.water The chemical composition of sea water
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item f.sea: A vector with the mixing ratio
#' \item mix: Matrix with the results of conservative mixing between fresh and sea water (in mmol/L)
#' \item react: Matrix with the concentrations due to reactions (in mmol/L)
#' }
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family base functions
#' @export calculate_convervative_mixing
calculate_convervative_mixing <- function(gdata, fresh.water = NULL, sea.water = NULL){
  #
  gmol <- c(40.078, 24.305, 22.989768, 39.0983, 61.01714, 60.0092, 35.4527, 96.0636)
  # mmol/L
  sea.water.comp <- c(10.7, 55.1, 485.0, 10.6, 2.4, 0, 566, 29.4)
  # mmo/L
  fresh.water.comp <- c(3.0, 0, 0, 0, 6, 0, 0, 0)
  #
  res.meqL <- convert_meql(gdata)
  mmol <- res.meqL$mmol
  ndat <- nrow(mmol)
  f.sea <- mmol[,7] / 566.0
  gdata.mix <- matrix(0.0, nrow = ndat, ncol = ncol(mmol))
  for(i in 1:nrow(mmol)){
    gdata.mix[i,] <- f.sea[i]*sea.water.comp + (1-f.sea[i])*fresh.water.comp
  }
  #
  gdata.react <- mmol - gdata.mix
  #
  gdata.mix.conc <- matrix(0.0, nrow = ndat, ncol = 8)
  gdata.react.conc <- matrix(0.0, nrow = ndat, ncol = 8)
  for(i in 1:8){
    gdata.mix.conc[,i] <- gdata.mix[,i]*gmol[i]
    gdata.react.conc[,i] <- gdata.react[,i]*gmol[i]
  }
  #
  res <- list(f.sea = f.sea, mix.mmol = gdata.mix, react.mmol = gdata.react,
              mix = gdata.mix.conc, react = gdata.react.conc)
  return(res)
}
#' @title
#' plot.geochemical_dataset
#' @description
#' Generic function to create different hydrogeochemical plots
#' @param x,y A geochemical_dataset object
#' @param ... Additional parameters to the plot function
#' @param type A character string with the type of hydrogeochemical plot to be created.
#' Currently the following plots are supported:
#' \itemize{
#' \item piper
#' \item ternary
#' \item durov
#' \item schoeller
#' \item stiff
#' \item multirectangular
#' \item radial
#' \item maucha
#' }
#' @param measure A character variable specifying the type of measure to be used in the plot.
#' Currently the types supported include:
#' \itemize{
#' \item 'conc': concentrations
#' \item 'meql': miliequivalent
#' }
#' @param vars Character vector. Variables to use in the plot. In the case of the Durov plot, these variables
#' are already defined. Used only for compatibility with the plot function.
#' @param color Character variable that specifies the variable to color the data inside the plot.
#' @param Size Character variable that specifies the variable to define the size of the data inside the plot.
#' @return
#' This function returns a ggplot2 object with the corresponding hydrogeochemical plot.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family base functions
#' @export
plot.geochemical_dataset <- function(x, y = NULL, ..., type = c('piper',
                                                                'ternary',
                                                                'durov',
                                                                'schoeller',
                                                                'stiff',
                                                                'multirectangular',
                                                                'radial'),
                                     measure = c('conc', 'meql'),
                                     vars = NULL,
                                     color = NULL,
                                     Size = NULL){
  base_cmd <- "plot_"
  current.cmd <- paste0(base_cmd, type)
  additional.args <- list(...)
  current.args <- list(x = x, measure = measure, vars = vars, color = color,
                       Size = Size, additional.args = additional.args)
  p <- do.call(current.cmd, args = current.args)
  return(p)
}
#' @title
#' plot_schoeller
#' @description
#' Function to create a Schoeller plot
#' @param x A geochemical dataset
#' @param measure The type of measure
#' @param vars variables
#' @param color color
#' @param Size  size
#' @return
#' A ggplot2 object with the corresponding Schoeller plot
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @export
#' @importFrom ggplot2 ggplot scale_x_discrete scale_y_log10
plot_schoeller <- function(x, measure = c('conc', 'meql'),
                           vars = NULL, color = NULL,
                           Size = NULL){
  gdata <- x
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  schoeller.df <- gdata$meqL.schoeller
  nsamples <- nrow(schoeller.df)/8
  if(is.null(color)){
    stop('ERROR: This function requires a color variable for grouping and also for color')
  }
  #print(schoeller.df)
  #
  p <- NULL
  tmp <- gdata$dataset[color]
  schoeller.df[color] <- tmp[,1]
  if(is.null(Size)){
    p <- ggplot(data =schoeller.df) + geom_point(aes_string(x = "major_ions",
                                                            y = "concentration",
                                                            color = color),
                                                 size = 3) +
      scale_x_discrete() +
      geom_line(aes_string(x = "major_ions", y = "concentration",
                           color = color, group = color),
                data = schoeller.df) +
      scale_y_log10()
  }
  else{
    tmp <- gdata$dataset[Size]
    schoeller.df[Size] <- tmp[,1]
    if(class(Size) == "numeric"){
      p <- ggplot(data =schoeller.df) + geom_point(aes_string(x = "major_ions",
                                                              y = "concentration",
                                                              color = color,
                                                              size = Size)) +
        scale_x_discrete() +
        geom_line(aes_string(x = "major_ions", y = "concentration",
                             color = color, group = color),
                  data = schoeller.df) +
        scale_y_log10()
    }
  }
  p <- p + theme_bw()
  return(p)
}
