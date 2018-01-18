#' @section Radial Diagram:
#'
#' The functions defined to create a Radial diagram are:
#' radial_tranform, ggplot_radial
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' radial_transform
#' @description
#' Function to calculate the coordinates of a geochemical_dataset in a Radial Diagram
#' @param gdata A geochemical_dataset object
#' @return
#' A data.frame with the corresponding coordinates
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family radial functions
#' @export radial_transform
radial_transform <- function(gdata){
  theta <- (2*pi/7)
  theta1 <- (pi/2)+c(seq(0,6),0)*theta
  yy <- sin(theta1)
  xx <- cos(theta1)
  r <- max(gdata$meqL.max)
  #
  current_sample <- gdata
  current_meqL <- current_sample$meqL/r
  current_coords <- c(current_meqL[1,1], current_meqL[1,2], current_meqL[1,3],
                      current_meqL[1,4],  current_meqL[1,5],
                      current_meqL[1,6] + current_meqL[1,7],  current_meqL[1,8],
                      current_meqL[1,1])
  current_coords <- rev(current_coords)
  ccoords.x <- current_coords*xx
  ccoords.y <- current_coords*yy
  coords.df <- data.frame(xc = ccoords.x, yc = ccoords.y)
  return(coords.df)
}
#' @title
#' ggplot_radial
#' @description
#' Function to create the template of a radial diagram
#' @return
#' A ggplot2 with the template of a Radial diagram
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family radial functions
#' @importFrom ggplot2 ggplot geom_segment theme_gray theme coord_fixed
#' @export ggplot_radial
ggplot_radial <- function(){
  xc <- NULL
  yc <- NULL
  ions <- NULL
  theta <- (2*pi/7)
  theta1 <- (pi/2)+c(seq(0,6),0)*theta
  yy <- sin(theta1)
  xx <- cos(theta1)
  #
  yyl <- 1.05*sin(theta1)
  xxl <- 1.05*cos(theta1)
  #
  vars <- c("Ca", "SO4", "Cl", "HCO3+CO3", "K", "Na", "Mg", "")
  coords.labels.df <- data.frame(xc = xxl, yc = yyl, ions = vars)
  #
  p <- ggplot() +
    geom_segment(aes(x=xx[1],y=yy[1]),xend=0,yend=0) +
    geom_segment(aes(x=xx[2],y=yy[2]),xend=0,yend=0) +
    geom_segment(aes(x=xx[3],y=yy[3]),xend=0,yend=0) +
    geom_segment(aes(x=xx[4],y=yy[4]),xend=0,yend=0) +
    geom_segment(aes(x=xx[5],y=yy[5]),xend=0,yend=0) +
    geom_segment(aes(x=xx[6],y=yy[6]),xend=0,yend=0) +
    geom_segment(aes(x=xx[7],y=yy[7]),xend=0,yend=0) +
    #
    geom_segment(aes(x=xx[1],y=yy[1]),xend=xx[2],yend=yy[2]) +
    geom_segment(aes(x=xx[2],y=yy[2]),xend=xx[3],yend=yy[3]) +
    geom_segment(aes(x=xx[3],y=yy[3]),xend=xx[4],yend=yy[4]) +
    geom_segment(aes(x=xx[4],y=yy[4]),xend=xx[5],yend=yy[5]) +
    geom_segment(aes(x=xx[5],y=yy[5]),xend=xx[6],yend=yy[6]) +
    geom_segment(aes(x=xx[6],y=yy[6]),xend=xx[7],yend=yy[7]) +
    geom_segment(aes(x=xx[7],y=yy[7]),xend=xx[8],yend=yy[8]) +
    # Add labels
    geom_text(aes(x = xc, y = yc, label = ions), data = coords.labels.df) +
    #
    coord_fixed() +
    #
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

  return(p)
}
#' @title
#' plot_radial
#' @description
#' Function to create a radial plot from a single sample
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
#' @return
#' This function returns a ggplot2 object with the Durov plot.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family radial functions
#' @export
plot_radial <- function(x, measure = c('conc', 'meql'),
                        vars = NULL, color = NULL,
                        Size = NULL){
  gdata <- x
  xc <- NULL
  yc <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  radial.df <- radial_transform(gdata)
  p <- ggplot_radial() + geom_polygon(aes(x = xc, y = yc), data = radial.df,
                                      fill = "blue", alpha = 0.4)
  return(p)
}

