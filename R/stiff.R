#' @section Stiff Diagram:
#'
#' The functions defined to create a Stiff diagram are:
#' stiff_tranform, ggplot_stiff
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' stiff_transform
#' @description
#' Function to calculate the coordinates of a geochemical dataset in a Stiff diagram
#' @param gdata A geochemical_dataset object
#' @return
#' A data.frame with the corresponding coordinates
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family stiff functions
#' @export stiff_transform
stiff_transform <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset object is required as input')
  }
  var <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  meqL <- gdata$meqL
  #print(meqL)
  meqL.max <- gdata$meqL.max
  meqL.min <- gdata$meqL.min
  meqL1 <- matrix(0.0, nrow = nrow(meqL), ncol = ncol(meqL))
  for(i in 1:8){
    if(abs(meqL.max[i]) > 1.0e-10){
      meqL1[,i] <- meqL[,i]/meqL.max[i]
    }
    else{
      meqL1[,i] <- meqL[,i]
    }
  }
  meqL1[,1:4] <- -1*meqL1[,1:4]
  results <- as.data.frame(meqL1)
  names(results) <- var
  nsamples <- nrow(meqL)
  xc <- matrix(0.0, nrow = 1, ncol = (7*nsamples))
  xc[1,1:7] <- c(results$Mg[1], results$Ca[1], results$Na[1]+results$K[1], results$SO4[1],
                 results$HCO3[1] + results$CO3[1], results$Cl[1], results$Mg[1])
  if(nsamples>1){
    for(isample in 2:nsamples){
      beginp <- 1+7*(isample-1)
      endp <- beginp+7-1
      xc[1,beginp:endp] <- c(results$Mg[isample], results$Ca[isample],
                             results$Na[isample]+results$K[isample],
                             results$SO4[isample], results$HCO3[isample] + results$CO3[isample],
                             results$Cl[isample], results$Mg[isample])
    }
  }
  #yc <- matrix(0.0, nrow=nrow(meqL), ncol = 7)
  ybase <- c(0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.1)
  #yc[1,] <- ybase
  yc <- rep(ybase, nrow(meqL))
  ions_base <- c('Mg', 'Ca', 'Na+K', 'SO4', 'HCO3+CO3', 'Cl', 'Mg')
  ions <-rep(ions_base, nsamples)
  coords <- data.frame(xc = t(xc), yc = yc, ions = ions)
  return(coords)
}
#' @title
#' ggplot_stiff
#' @description
#' Function to create a Stiff diagram
#' @return
#' A ggplot2 plot with the Stiff diagram template
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family stiff functions
#' @importFrom ggplot2 ggplot geom_segment geom_text theme_bw theme
#' @export ggplot_stiff
ggplot_stiff <- function(){
  p <- ggplot() +
    ## Main axis
    geom_segment(aes(x = 0.0,y = 0.075, xend = 0.0, yend = 0.225)) +
    geom_segment(aes(x = -1,y = 0.225, xend = 1, yend = 0.225)) +
    # Cations
    geom_text(aes(x = -1.1,y = 0.1, label = "Mg"), size = 4) +
    geom_text(aes(x = -1.1,y = 0.15, label = "Ca"), size = 4) +
    geom_text(aes(x = -1.1,y = 0.2, label = "Na+K"), size = 4) +
    # Anions
    geom_text(aes(x = 1.1,y = 0.1, label = "Cl"), size = 4) +
    geom_text(aes(x = 1.1,y = 0.15, label = "HCO3+CO3"), size = 4) +
    geom_text(aes(x = 1.1,y = 0.2, label = "SO4"), size = 4) +
    #
    geom_segment(aes(x = -1, y = 0.1, yend = 0.1, xend = 1), linetype = "dashed",
                 size = 0.25, colour = "grey50") +
    geom_segment(aes(x = -1, y = 0.15, yend = 0.15, xend = 1), linetype = "dashed",
                 size = 0.25, colour = "grey50") +
    geom_segment(aes(x = -1, y = 0.2, yend = 0.2, xend = 1), linetype = "dashed",
                 size = 0.25, colour = "grey50") +
    #
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.ticks = element_blank())

  return(p)
}
#' @title
#' plot_stiff
#' @description
#' Function to create a Stiff plot for a given sample
#' @param x A geochemical_dataset object
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
#' @param labels Character variable that specifies the labels to be used in the current plot
#' @param additional.args A list with additional arguments
#' @return
#' This function returns a ggplot2 object with the Durov plot.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family stiff functions
#' @importFrom ggplot2 geom_polygon
#' @export
plot_stiff <- function(x, measure = c('conc', 'meql'),
                       vars = NULL, color = NULL,
                       Size = NULL, labels = NULL, additional.args = NULL){
  gdata <- x
  xc <- NULL
  yc <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  stiff.df <- stiff_transform(gdata)
  p <- ggplot_stiff()
  p <- p + geom_polygon(aes(x = xc, y = yc), data=stiff.df,
                        color = "red", fill = "blue", alpha = 0.5)
  return(p)
}
