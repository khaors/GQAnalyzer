#' @section Compositional-ILR Diagram:
#'
#' The functions defined to create a diagram with an ILR transform are:
#' ilr_transform, plot_ilr
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' ilr_transform
#' @description
#' Function to calculate the coordinates of a geochemical_dataset in a multirectangular diagram.
#' @param gdata A geochemical_dataset object
#' @return
#' A data.frame with the coordinates
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @source
#' Shelton, Jenna L., Engle, Mark A., Buccianti, Antonella, Blondes, Madalyn S. (2018).
#' The isometric log-ratio (ilr)-ion plot: A proposed alternative to the Piper diagram.
#' ournal of Geochemical Exploration, Vol 190, Issue March.
#' @family ilr functions
#' @export ilr_transform
ilr_transform <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset object is required as input')
  }
  #
  cations_var <- c("Ca", "Mg", "Na", "K")
  anions_var <- c("HCO3", "CO3", "Cl", "SO4")
  anions <- gdata$anions
  cations <- gdata$cations
  #
  tol <- 1e-3
  term11 <- (gdata$dataset$Ca+tol)*(gdata$dataset$Mg+tol)
  term12 <- gdata$dataset$Na+gdata$dataset$K+2*tol
  z1 <- sqrt(2/3)*log(sqrt(term11)/(term12))
  term21 <- (gdata$dataset$Ca+tol)/(gdata$dataset$Mg+tol)
  z2 <- (1/sqrt(2))*log(term21)
  term31 <- sqrt((gdata$dataset$Cl+tol)*(gdata$dataset$SO4+tol))
  term32 <- (gdata$dataset$HCO3+gdata$dataset$CO3+2*tol)
  z3 <- sqrt(2/3)*log(term31/term32)
  z4 <- (1/sqrt(2))*log((gdata$dataset$Cl+tol)/(gdata$dataset$SO4+tol))
  #
  res <- data.frame(z1 = z1, z2 = z2, z3 = z3, z4 = z4)
  return(res)
}
#' @title
#' plot_ilr
#' @description
#' This function creates the compositional plot
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
#' @family ilr functions
#' @importFrom ggplot2 ggplot geom_point geom_smooth coord_equal geom_qq
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @export plot_ilr
plot_ilr <- function(x, measure = c('conc', 'meql'),
                     vars = NULL, color = NULL,
                     Size = NULL){
  gdata <- x
  z1 <- NULL
  z2 <- NULL
  z3 <- NULL
  z4 <- NULL
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  ilr.df <- ilr_transform(x)
  p <- NULL
  p1 <- NULL
  p2 <- NULL
  p3 <- NULL
  p4 <- NULL
  if(is.null(color)){
    if(is.null(Size)){
      print(c("nsize","ncolor"))
      p1 <- ggplot() + geom_point(aes(x = z2, y = z4),
                               data = ilr.df, size = 3) +
       theme_bw()
      p2 <- ggplot() + geom_point(aes(x = z1, y = z4),
                                data = ilr.df, size = 3) +
        theme_bw()
      p3 <- ggplot() + geom_point(aes(x = z2, y = z3),
                                  data = ilr.df, size = 3) +
        theme_bw()
      p4 <- ggplot() + geom_point(aes(x = z1, y = z3),
                                  data = ilr.df, size = 3) +
        theme_bw()
      p <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
    }
    else{
      if(class(Size) == "numeric"){
        p1 <- ggplot() + geom_point(aes(x = z2, y = z4),
                                  data = ilr.df, size = Size) +
          theme_bw()
        p2 <- ggplot() + geom_point(aes(x = z1, y = z4),
                                    data = ilr.df, size = Size) +
          theme_bw()
        p3 <- ggplot() + geom_point(aes(x = z2, y = z3),
                                    data = ilr.df, size = Size) +
          theme_bw()
        p4 <- ggplot() + geom_point(aes(x = z1, y = z3),
                                    data = ilr.df, size = Size) +
          theme_bw()
        p <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        ilr.df[Size] <- tmp[,1]
        p1 <- ggplot() + geom_point(aes_string(x = "z2", y = "z4",
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p2 <- ggplot() + geom_point(aes_string(x = "z1", y = "z4",
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p3 <- ggplot() + geom_point(aes_string(x = "z2", y = "z3",
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p4 <- ggplot() + geom_point(aes_string(x = "z1", y = "z3",
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
      }
    }
  }
  else {
    tmp <- gdata$dataset[color]
    ilr.df[color] <- tmp[,1]
    if(is.null(Size)){
      p1 <- ggplot() + geom_point(aes_string(x = "z2", y = "z4",
                                           color = color),
                                data = ilr.df, size = 3) +
        theme_bw()
      p2 <- ggplot() + geom_point(aes_string(x = "z1", y = "z4",
                                           color = color),
                                data = ilr.df, size = 3) +
        theme_bw()
      p3 <- ggplot() + geom_point(aes_string(x = "z2", y = "z3",
                                           color = color),
                                data = ilr.df, size = 3)+
        theme_bw()
      p4 <- ggplot() + geom_point(aes_string(x = "z1", y = "z3",
                                           color = color),
                                data = ilr.df, size = 3)+
        theme_bw()
      if(class(tmp[,1]) == "numeric"){
        p1 <- p1 + scale_color_gradientn(colours = rainbow(10))
        p2 <- p2 + scale_color_gradientn(colours = rainbow(10))
        p3 <- p3 + scale_color_gradientn(colours = rainbow(10))
        p4 <- p4 + scale_color_gradientn(colours = rainbow(10))
      }
      else if(class(tmp[,1] == "factor")){
        p1 <- p1 + scale_color_discrete()
        p2 <- p2 + scale_color_discrete()
        p3 <- p3 + scale_color_discrete()
        p4 <- p4 + scale_color_discrete()
      }
      p <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
    }
    else{
      if(class(Size) == "numeric"){
        p1 <- ggplot() + geom_point(aes_string(x = "z2", y = "z4",
                                             color = color),
                                  data = ilr.df, size = Size) +
          theme_bw()
        p2 <- ggplot() + geom_point(aes_string(x = "z1", y = "z4",
                                             color = color),
                                  data = ilr.df, size = Size) +
          theme_bw()
        p3 <- ggplot() + geom_point(aes_string(x = "z2", y = "z4",
                                             color = color),
                                  data = ilr.df, size = Size) +
          theme_bw()
        p4 <- ggplot() + geom_point(aes_string(x = "z1", y = "z4",
                                             color = color),
                                  data = ilr.df, size = Size) +
          theme_bw()
        p <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        ilr.df[Size] <- tmp[,1]
        p1 <- ggplot() + geom_point(aes_string(x = "z2", y = "z4",
                                             color = color,
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p2 <- ggplot() + geom_point(aes_string(x = "z1", y = "z4",
                                             color = color,
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p3 <- ggplot() + geom_point(aes_string(x = "z2", y = "z3",
                                             color = color,
                                             size = Size),
                                  data = ilr.df) +
          theme_bw()
        p4 <- ggplot() + geom_point(aes_string(x = "z1", y = "z3",
                                             color = color,
                                             size = Size),
                                  data = ilr.df)+
          theme_bw()
        if(class(tmp[,1]) == "numeric"){
          p1 <- p1 + scale_color_gradientn(colours = rainbow(10))
          p2 <- p2 + scale_color_gradientn(colours = rainbow(10))
          p3 <- p3 + scale_color_gradientn(colours = rainbow(10))
          p4 <- p4 + scale_color_gradientn(colours = rainbow(10))
        }
        else if(class(tmp[,1] == "factor")){
          p1 <- p1 + scale_color_discrete()
          p2 <- p2 + scale_color_discrete()
          p3 <- p3 + scale_color_discrete()
          p4 <- p4 + scale_color_discrete()
        }
        p <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
       }
    }
  }
  # #
  return(p)
}
