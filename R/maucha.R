#' @section Maucha Diagram:
#'
#' The functions defined to create a Multirectangular diagram are:
#' maucha_transform, ggplot_maucha
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' maucha_transform
#' @description
#' Function to calculate the coordinates of a given sample in a Maucha plot
#' @param gdata A geochemical_dataset object
#' @return
#' This function returns a data.frame with the following variables:
#' \itemize{
#' \item xcirc
#' \item ycirc
#' }
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family maucha functions
#' @export
maucha_transform <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
       stop('ERROR: a geochemical_dataset object is required as input')
  }
  colours <-  c("yellow2", "purple", "black", "cyan", "green", "blue", "grey", "red")
  mlabels <- c("Na", "K", "x", "TAL", "Cl", "SO4", "Mg", "Ca")
  ion <- gdata$meqL
  total_anions <- apply(gdata$anions, 2, sum)
  total_cations <- apply(gdata$cations, 2, sum)
  total_ion <- sum(total_anions) + sum(total_cations)
  # Create the circle
  sector_angle <- 360/16
  A <- total_ion
  sector <- sector_angle * ( pi / 180 )
  R <- ( (A / 16) *2 / sin(sector) ) ^ 0.5
  circleSector <- sector / 4
  x.circ <- array(1:64)
  y.circ <- array(1:64)
  for (i in 1:64 ) {
    x.circ[i] <- R * cos(circleSector*i)
    y.circ[i] <- R * sin(circleSector*i)
  }
  x.circ <- c(x.circ, x.circ[1])
  y.circ <- c(y.circ, y.circ[1])
  # Plot Maucha Symbol
  count <- 1
  polygons <- vector("list", length = 8)
    for (i in 1:8 ) { # construct the kites, one for each variable
      ion.conc <- ion[i]
      #For each ion, construct a
      Ai <- ion.conc / total_ion * A
      #radial shape made up of kite-
      h <- count - 1
      #shaped quadrilaterals, then
      j <- count + 1
      #offset it by xm, ym.
      x.c <- R * cos(sector*h)
      y.c <- R * sin(sector*h)
      #x.d,y.d
      #x.b,y.b
      x.d <- R * cos(sector*j)
      y.d <- R * sin(sector*j)
      a <- Ai / ( R * sin(sector) )
      x.b <- a * cos(sector*count)
      y.b <- a * sin(sector*count)
      x.t <- 1.5 * max(a,R) * cos(sector*count)
      y.t <- 1.5 * max(a,R) * sin(sector*count)
      x <- c(0, x.d, x.b, x.c, 0)
      y <- c(0, y.d, y.b, y.c, 0)
      polygons[[i]] <- list(x = x, y = y, color = colours[[i]], label = mlabels[i])
      #polygon( x, y, col=colours[[i]], border = NA )
      #x.t <- x.t
      #y.t <- y.t
      #if(label) print(paste(x.t,y.t,label,mlabels[i]))
      #if(label) text(x.t,y.t,mlabels[i],cex=0.3)
      count <- count + 2
    }
  #
  res <- list(x.circ = x.circ, y.circ = y.circ, polygons = polygons)
  return(res)
}
#' @title
#' ggplot_maucha
#' @description
#' Function to create the template of the Maucha plot.
#' @return
#' This function returns a gggplot2 object with the template of the Maucha plot.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family maucha functions
#' @importFrom ggplot2 ggplot geom_line coord_equal theme_bw theme
#' @export
ggplot_maucha <- function(){
  R <- 0.8082583
  dx <- R/100
  x <- seq(-R, R, by = dx)
  y <- NULL
  y1 <- sqrt(R^2-x^2)
  y2 <- -1*y1
  circ.df <- data.frame(x = x, y = y1)
  circ.df1 <- data.frame(x = x, y = -1*y1)
  p1 <- ggplot() + geom_line(aes(x = x, y = y), data = circ.df, col = "black") +
    geom_line(aes(x = x, y = y), data = circ.df1, col = "black") +
    coord_equal() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  return(p1)
}
#' @title
#' plot_maucha
#' @description
#' Function to create a Maucha plot.
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
#' @family maucha functions
#' @importFrom ggplot2 ggplot geom_polygon
#' @export
plot_maucha <- function(x, measure = c('conc', 'meql'),
                        vars = NULL, color = NULL,
                        Size = NULL){
  gdata <- x
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  p <- ggplot_maucha()
  y <- NULL
  maucha.df <- maucha_transform(gdata)
  for(i in 1:8){
    pol.df <- as.data.frame(maucha.df$polygons[[i]])
    p <- p + geom_polygon(aes(x = x, y = y), data = pol.df, fill = pol.df$color[1],
                          color = "black")
  }
  return(p)
}
