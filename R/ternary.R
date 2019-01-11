#' @section Ternary Plot:
#'
#' The functions defined to create a Ternary plot are:
#' ternary_transform, ggplot_ternary
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' ternary_transform
#' @description
#' Function to calculate the coordinates of a point in a ternary diagram.
#' @param a,b,c Concentration values
#' @return
#' This function returns a data.frame with the coordinates inside a ternary diagram.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family ternary functions
#' @export ternary_transform
ternary_transform <- function(a, b, c){
  sumall <- a + b + c
  a_norm <- a/sumall
  b_norm <- b/sumall
  c_norm <- c/sumall
  xc <- 0.5*(2*b+c)/sumall
  yc <- 0.5*sqrt(3)*c/sumall
  results <- as.data.frame(list(xc = xc, yc = yc))
  names(results) <- c("xc", "yc")
  return(results)
}
#' @title
#' ggplot_ternary
#' @description
#' Function to create the ternary diagram
#' @return
#' This function returns the base of a ternary diagram
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family ternary functions
#' @importFrom ggplot2 ggplot geom_segment theme_bw theme coord_equal
#' @export ggplot_ternary
ggplot_ternary <- function(){
  #
  x1 <- NULL
  x2 <- NULL
  y1 <- NULL
  y2 <- NULL
  #
  p <- ggplot() +
    ## left hand ternary plot
    geom_segment(aes(x=0,y=0, xend=1, yend=0)) +
    geom_segment(aes(x=0,y=0, xend=.50, yend=.86603)) +
    geom_segment(aes(x=.50,y=.86603, xend=1, yend=0)) +
    #
    #
    coord_equal(ratio = 0.9) +
    #
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  #
  return(p)
}
#' @title 
#' add_grid_lines_ternary
#' @description 
#' Function to add the grid lines to a ternary diagram
#' @param color A character string specifying the color of the grid lines in a 
#' ternary plot
#' @param Size A numeric value specifying the size of the grid labels in a 
#' ternary plot
#' @return
#' This function returns a list with the geom components used to create a grid 
#' in a ternary diagram
#' @author 
#' Oscar Garcia-Cabrejo, \email{khaors@gmail,com}
#' @family ternary functions 
#' @importFrom ggplot2 geom_segment geom_text
#' @export
add_grid_lines_ternary <- function(color = NULL, Size = NULL){
  #
  x1 <- NULL
  x2 <- NULL
  y1 <- NULL
  y2 <- NULL
  grid1p1 <- NULL
  grid1p2 <- NULL
  grid1p3 <- NULL
  #
  grid1p1 <- data.frame(x1 = c(.20,.40,.60,.80), x2= c(.10,.20,.30,.40),
                        y1 = c(0,0,0,0), y2 = c(.173206,.346412,.519618, .692824))
  grid1p2 <- data.frame(x1 = c(.20,.40,.60,.80), x2= c(.60,.70,.80,.90),
                        y1 = c(0,0,0,0), y2 = c(.692824, .519618,.346412,.173206))
  grid1p3 <- data.frame(x1 = c(.10,.20,.30,.40), x2= c(.90,.80,.70,.60),
                        y1 = c(.173206, .346412, .519618, .692824),
                        y2 = c(.173206, .346412, .519618, .692824))
  #
  color1 <- "grey50"
  Size1 <- 3
  if(!missing(color)){
    color1 <- color
  }
  if(!missing(Size)){
    Size1 <- Size
  }
  #
  res <- list()
  #
  res[[1]] <- geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), 
                           data=grid1p1, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[2]] <- geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), 
                           data=grid1p2, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[3]] <- geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), 
                           data=grid1p3, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[4]] <- geom_text(aes(c(.20, .40, .60, .80), 
                            c(-.05, -.05, -.05, -.05), 
                            label=c(20, 40, 60, 80)), size=Size1)
  res[[5]] <- geom_text(aes(c(.35, .25, .15, .05), 
                            grid1p2$y2, label=c(20, 40, 60, 80)), size=Size1)
  res[[6]] <- geom_text(aes(c(.95, .85, .75, .65),
                            grid1p3$y2, label=c(20, 40, 60, 80)), size=Size1)
  #
  return(res)
}
#' @title 
#' add_labels_ternary
#' @description 
#' Function to add the labels of a ternary diagram
#' @param vars A character vector with the names of the variables
#' @param middle A logical variable specifying if the labels are located in the 
#' central part of the ternary axis (default=TRUE)
#' @param color A character string specifying the color of the grid lines in a 
#' ternary plot
#' @param Size A numeric value specifying the size of the grid labels in a 
#' ternary plot
#' @return 
#' This function returns a list with the geom_text component that specifies the 
#' labels of a ternary plot
#' @author 
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family ternary functions
#' @importFrom ggplot2 geom_text
#' @export 
add_labels_ternary <- function(vars, middle = TRUE, color = NULL, Size = NULL){
  poslabels <- NULL
  x <- NULL
  y <- NULL
  angle <- NULL
  if(middle){
    poslabels <- list(x = c(.13,.83,0.5), 
                      y = c(.45, .45, -0.1),
                      vars = c(vars[3],vars[2],vars[1]),
                      angle = c(60, -60, 0.0))
  }
  else{
    poslabels <- list(x = c(-0.1,.5,1.1), 
                      y = c(0, .87, 0.0),
                      vars = c(vars[3],vars[2],vars[1]),
                      angle = c(0, 0, 0.0))
  }
  poslabels <- as.data.frame(poslabels)
  #
  res <- list()
  res[[1]] <- geom_text(aes(x = x, y = y, label = vars, angle = angle), 
                        data = poslabels)
  return(res)
}
#' @title
#' plot_ternary
#' @description
#' Function to create a ternary plot for a geochemical_dataset
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
#' This function returns a ggplot2 object with the Durov plot.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family ternary functions
#' @importFrom ggplot2 geom_point scale_color_gradientn scale_color_discrete
#' @export
plot_ternary <- function(x, measure = c('conc', 'meql'),
                         vars = NULL, color = NULL,
                         Size = NULL, 
                         additional.args = NULL){
  gdata <- x
  p <- NULL
  angle <- NULL
  xc <- NULL
  yc <- NULL
  y <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  if(missing(vars) | is.null(vars)){
    stop('ERROR: Variable names must be specified to create a ternary diagram')
  }
  pos <- vector('numeric', 3)
  ternary.df <- NULL
  if(measure == "meql"){
    for(ivar in 1:3){
      test <- vars[ivar]==meql_ions
      pos[ivar] <- which.max(test)
    }
    ternary.df <- ternary_transform(gdata$meqL[,pos[1]],
                                    gdata$meqL[,pos[2]],
                                    gdata$meqL[,pos[3]])
  }
  else if(measure == "conc"){
    for(ivar in 1:3){
      test <- vars[ivar]==conc_ions
      pos[ivar] <- which.max(test)
    }
    ternary.df <- ternary_transform(gdata$dataset[,pos[1]],
                                    gdata$dataset[,pos[2]],
                                    gdata$dataset[,pos[3]])
  }
  #
  p <- ggplot_ternary() + add_grid_lines_ternary() + add_labels_piper(vars)
  if(is.null(color)){
    if(is.null(Size)){
      p <- p + geom_point(aes(x = xc, y = yc), data = ternary.df, size = 3)
    }
    else{
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes(x = xc, y = yc), data = ternary.df, size = Size)
      }
      else if(class(Size) == "character"){
        ternary.df[Size] <- gdata$dataset[Size]
        p <- p + geom_point(aes(x = xc, y = yc, size = Size), data = ternary.df)
      }
    }    }
  else{
    tmp <- gdata$dataset[color]
    ternary.df[color] <- tmp[,1]
    if(is.null(Size)){
      p <- p + geom_point(aes_string(x = "xc", y = "yc", color = color),
                          data = ternary.df, size = 3)
    }
    else {
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes_string(x = "xc", y = "yc", color = color),
                            data = ternary.df, size = Size)
      }
      else if(class(Size) == "character"){
        ternary.df[Size] <- gdata$dataset[Size]
        p <-  p + geom_point(aes_string(x = "xc", y = "yc", color = color,
                                        size = Size), data = ternary.df)
      }
    }
  }
  p <- p + coord_fixed()

  if(!is.null(color)){
    tmp <- gdata$dataset[color]
    if(class(tmp[,1]) == "numeric"){
      p <- p + scale_colour_gradientn(colours = rainbow(20))
    }
    else if(class(tmp[,1]) == "factor"){
      p <- p + scale_color_discrete()
    }
  }
  return(p)
}
