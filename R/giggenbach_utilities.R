#' @section Giggenbach Diagrams:
#'
#' The functions defined to create the Giggenbach diagrams are:
#'
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' giggenbach_ternary_transform
#' @description
#' This function calculates the coordinates of the ternary plot proposed by Giggenbach (1988)
#' to identify the samples that can be used in the geothermometer calculations
#' @param gdata A geochemical_dataset object
#' @return
#' This function returns a data.frame with the following entries
#' \itemize{
#' \item cNa: concentrations of Na divided by 1000
#' \item cK: concentration of K divided by 100
#' \item cMg: square root of concentration of Mg
#' }
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family giggenbach_functions
#' @export
giggenbach_ternary_transform <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset object is required as input')
  }
  cations_var <- c("Ca", "Mg", "Na", "K")
  anions_var <- c("HCO3", "CO3", "Cl", "SO4")
  #
  cNa <- gdata$dataset$Na/1000
  cK <- gdata$dataset$K/100
  cMg <- sqrt(gdata$dataset$Mg)
  #
  ternary <- ternary_transform(cK, cMg, cNa)
  results <- ternary
}
#' @title
#' ggplot_giggenbach
#' @description
#' Function to create the base of the Giggenbach's geothermometer plot
#' @return
#' This function returns the base of the Giggenbach's geothermometer plot
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family giggenbach functions
#' @importFrom ggplot2 ggplot geom_segment
#' @export ggplot_giggenbach
ggplot_giggenbach <- function(){
  xc <- NULL
  yc <- NULL
  xi <- NULL
  yi <- NULL
  p <- ggplot_ternary() +   
    geom_segment(aes(x = xi, y = yi, xend = xc, yend=yc), 
                 data = giggenbach_lines1, 
                 color = "gray") +
    geom_segment(aes(x = xi, y = yi, xend = xc, yend=yc), 
                 data = giggenbach_lines2, color = "gray")
  return(p)
}
#' @title 
#' add_lines_giggenbach
#' @description 
#' Auxiliary function to add the lines of the inmature and mature waters 
#' in a Giggenbach's geothermometer plot.
#' @param color A character string specifying the color of the lines separating 
#' mature and inmature waters.
#' @param Size A numeric value specifying the size of the temperature labels 
#' in the Giggenbach's geothermometer plot 
#' @param dy A numeric value specifying the vertical displacement of the 
#' temperature labels in the Giggenbach's geothermometer plot 
#' @return 
#' This function adds lines to define separate fields in the Giggenbach's 
#' geothermometer plot that has been previously created with the function 
#' ggplot_giggenbach
#' @author Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family giggenbach functions
#' @importFrom ggplot2 geom_path geom_point geom_text
#' @export 
add_lines_giggenbach <- function(color = NULL, Size = NULL, dy = NULL){
  xc <- NULL
  yc <- NULL
  if(is.null(color)){
    color <- "red"
  }
  if(is.null(Size)){
    Size <- 3
  }
  if(is.null(dy)){
    dy <- 0.025
  }
  res <- list()
  res[[1]] <- geom_path(aes(x = xc, y = yc), data = giggenbach_limits, 
                         color = color)
  res[[2]] <- geom_point(aes(x = xc, y = yc), data = ternary.df2a, 
              color = color)
  res[[3]] <- geom_text(aes(x = xc, y = yc + dy, label = as.character(tdat$Temp)), 
              data = ternary.df2a, size = Size, col = color)
  return(res)
}
#' @title
#' plot_giggenbach
#' @description
#' Function to create the ternary plot proposed by Giggenbach(1988) to identify the samples
#' to be used in geothermometer calculations
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
#' This function returns a ggplot2 object with the Giggenbach ternary plot.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family giggenbach functions
#' @export
plot_giggenbach <- function(x, measure = c('conc', 'meql'),
                            vars = NULL, color = NULL,
                            Size = NULL, labels = NULL,
                            additional.args = NULL){
  gdata <- x
  y <- NULL
  cx <- NULL
  cy <- NULL
  yc <- NULL
  xc <- NULL
  angle <- NULL
  poslabels <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  giggenbach.df <- giggenbach_ternary_transform(gdata)
  ternary.df <- ternary_transform(giggenbach.df$cNa,
                                  giggenbach.df$cK,
                                  giggenbach.df$cMg)
  p <- ggplot_ternary()
  if(is.null(color)){
    if(is.null(Size)){
      p <- p + geom_point(aes(x = xc, y = yc), data = ternary.df, size = 3)
    }
    else{
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes(x = xc, y = yc), data = ternary.df, 
                            size = Size)
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
  p <- p + geom_text(aes(x = x, y = y, label = vars, angle = angle),
                     data = poslabels) +
    coord_fixed()

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
