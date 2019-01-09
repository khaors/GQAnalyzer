#' @section Multirectangular Diagram:
#'
#' The functions defined to create a Multirectangular diagram are:
#' multirectangular_transform, ggplot_multirectangular
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' multirectangular_transform
#' @description
#' Function to calculate the coordinates of a geochemical_dataset in a multirectangular diagram.
#' @param gdata A geochemical_dataset object
#' @return
#' A data.frame with the coordinates
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family multirectangular functions
#' @export multirectangular_transform
multirectangular_transform <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset object is required as input')
  }
  #
  cations_var <- c("Ca", "Mg", "Na", "K")
  anions_var <- c("HCO3", "CO3", "Cl", "SO4")
  anions <- gdata$anions
  cations <- gdata$cations
  anions.pmx <- apply(anions, 1, which.max)
  cations.pmx <- apply(cations, 1, which.max)
  anions.mx <- apply(anions, 1, max)
  cations.mx <- apply(cations, 1, max)
  #
  ndat <- nrow(gdata$dataset)
  class.anions <- c("Carbonatada/Bicarbonatada", "Clorurada", "Sulfatada")
  class.cations <- c("Calcica", "Magnesica", "Sodica/Potasica")
  classification <- vector("character", length = ndat)
  for(i in 1:ndat){
    classification[i] <- paste0(class.anions[anions.pmx[i]], "-", class.cations[cations.pmx[i]])
  }
  #
  danions <- anions.pmx - 1
  dcations <- cations.pmx - 1
  cy <- anions.mx + danions
  cx <- cations.mx + dcations

  tmp <- NULL
  results <- data.frame(cx = cx, cy = cy, da = anions.pmx, dc = cations.pmx,
                        class = classification)
  return(results)
}
#' @title
#' ggplot_multirectangular
#' @description
#' Function to create the template of a multirectangular diagram.
#' @return
#' A ggplot2 plot with the template of a multirectangular diagram.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family multirectangular functions
#' @importFrom ggplot2 ggplot geom_segment theme_gray theme coord_fixed
#' @export
ggplot_multirectangular <- function(){
  p <- ggplot() +
    # Large square
    geom_segment(aes(x = 0,y = 0, xend = 3, yend = 0)) +
    geom_segment(aes(x = 3,y = 0, xend = 3, yend = 3)) +
    geom_segment(aes(x = 3,y = 3, xend = 0, yend = 3)) +
    geom_segment(aes(x = 0,y = 3, xend = 0, yend = 0)) +
    # lines to subdivide large square
    geom_segment(aes(x = 1,y = 0, xend = 1, yend = 3)) +
    geom_segment(aes(x = 2,y = 0, xend = 2, yend = 3)) +
    # lines to subdivide large square
    geom_segment(aes(x = 0,y = 1, xend = 3, yend = 1)) +
    geom_segment(aes(x = 0,y = 2, xend = 3, yend = 2)) +
    # Add tick marks to each square in x direction
    geom_text(aes(c(0.025, 0.5, 0.975), rep(-.1,3), label = c(0,50,100)),size=3) +
    geom_text(aes(c(1.025, 1.5, 1.975), rep(3.1, 3), label = c(0,50,100)),size=3) +
    geom_text(aes(c(2.025, 2.5, 2.975), rep(-.1,3), label = c(0,50,100)),size=3) +
    #
    geom_text(aes(rep(-.1,3), c(0.025, 0.5, 0.975), label = c(0,50,100)),size=3) +
    geom_text(aes(rep(3.1, 3), c(1.025, 1.5, 1.975), label = c(0,50,100)),size=3) +
    geom_text(aes(rep(-.1,3),c(2.025, 2.5, 2.975), label = c(0,50,100)),size=3) +
    #
    theme_bw() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_blank(), axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
  coord_equal()
  #
  return(p)
}
#' @title 
#' add_labels_multirectangular
#' @description 
#' Function to add the labels to a multirectangular diagram
#' @param color A character string specifying the color of the labels
#' @param Size A numeric value specifying the size of the text label
#' @return 
#' This function returns a geom_text to be added to the multirectangular
#'  diagram
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family multirectanguar functions
#' @importFrom ggplot2 geom_text
#' @export
add_labels_multirectangular <- function(color = NULL, 
                                        Size = NULL){
  color1 <- "black"
  Size1 <- 4
  if(!missing(color)){
    color1 <- color
  }
  #
  if(!missing(Size)){
    Size1 <- Size
  }
  #
  res <- list() 
  # Add labels to each square in x direction
  res[[1]] <- geom_text(aes(x = 0.5, y = -0.2, label = "%Ca(meq/L)"), 
                        size = Size1, colour = color1)
  res[[2]] <- geom_text(aes(x = 1.5, y = -0.2, label = "%Mg(meq/L)"), 
                        size = Size1, colour = color1)
  res[[3]] <- geom_text(aes(x = 2.5, y = -0.2, label = "%Na+K(meq/L)"), 
                        size = Size1, colour = color1)
    # Add labels to each square in y direction
  res[[4]] <- geom_text(aes(x = -0.2, y = 0.5, label = "%HCO3+CO3(meq/L)",
                            angle = 90), size = Size1, colour = color1)
  res[[5]] <- geom_text(aes(x = -0.2, y = 1.5, label = "%Cl(meq/L)", 
                            angle = 90), size = Size1, colour = color1)
  res[[6]] <- geom_text(aes(x = -0.2, y = 2.5, label = "%SO4(meq/L)", 
                            angle = 90), size = Size1, colour = color1)
  return(res)  
}
#

#' @title
#' plot_multirectangular
#' @description
#' Function to create the Multirectangular plot
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
#' @param additional.args A list with additional arguments 
#' @return
#' This function returns a ggplot2 object with the Multirectangular plot.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family multirectangular functions
#' @importFrom ggplot2 ggplot geom_point coord_fixed scale_color_gradientn scale_colour_gradientn scale_color_discrete aes_string
#' @importFrom grDevices rainbow
#' @export
#' @references
#' Ahmad, N., Sen, Z., & Ahmad, M. (2003). Ground Water Quality Assessment Using Multi-Rectangular Diagrams.
#' Ground Water, 41(6), 828–832. http://doi.org/10.1111/j.1745-6584.2003.tb02423.x
plot_multirectangular <- function(x, measure = c('conc', 'meql'),
                                  vars = NULL, color = NULL,
                                  Size = NULL, 
                                  additional.args = NULL){
  gdata <- x
  cx <- NULL
  cy <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  multirec.df <- multirectangular_transform(gdata)
  p <- ggplot_multirectangular() + add_labels_multirectangular()
  if(is.null(color)){
    if(is.null(Size)){
      p <- p + geom_point(aes(x = cx, y = cy), data = multirec.df, size = 3)
    }
    else{
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes(x = cx, y = cy), data = multirec.df, 
                            size = Size)
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        multirec.df[Size] <- tmp[,1]
        p <- p + geom_point(aes_string(x = "cx", y = "cy", size = Size),
                            data = multirec.df)
      }
    }
  }
  else {
    tmp <- gdata$dataset[color]
    multirec.df[color] <- tmp[,1]
    if(is.null(Size)){
      p <- p + geom_point(aes_string(x = "cx", y = "cy", color = color),
                          data = multirec.df, size = 3)
    }
    else{
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes_string(x = "cx", y = "cy", color = color),
                            data = multirec.df, size = Size)
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        multirec.df[Size] <- tmp[,1]
        p <- p + geom_point(aes_string(x = "cx", y = "cy", color = color,
                            size = Size), data = multirec.df)
      }
    }
  }
  #
  if(!is.null(color)){
    tmp <- gdata$dataset[color]
    if(class(tmp[,1]) == "numeric"){
      p <- p + scale_colour_gradientn(colours = rainbow(10))
    }
    else if(class(tmp[,1]) == "factor"){
      p <- p + scale_color_discrete()
    }
  }
  #
  return(p)
}

