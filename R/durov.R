#' @section Durov Diagram:
#'
#' The Durov diagram is a composite diagram that includes a central rectangular diagram of
#'
#' This section includes the definition of three functions used to create the Durov diagram:
#' durov_transform, durov_transform1, ggplot_durov
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' durov_transform1
#' @description
#' Function to calculate the coordinates of a given composition in a Durov Diagram.
#' @param gdata A geochemical_dataset object
#' @return
#' A data.frame with the coordinates
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family durov functions
#' @export durov_transform1
durov_transform1 <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical dataset object is required as input')
  }
  # Main square: Ca vs HCO3+CO3
  cx_mnsquare <- gdata$cations[,1]
  cy_mnsquare <- gdata$anions[,1]
  # pH Coordinates: y in (-1,0)
  pH.mn <- min(gdata$dataset$pH)
  pH.mx <- max(gdata$dataset$pH)
  #print(c(pH.mn,pH.mx))
  pH_scaled <- -1+(gdata$dataset$pH-pH.mn)/((pH.mx-pH.mn))
  cx_lwsquare <- gdata$cations[,1]
  cy_lwsquare <- pH_scaled
  # TDS coordinates: x in (1,2) y in (0,1)
  TDS.mn <- min(gdata$dataset$TDS)
  TDS.x <- max(gdata$dataset$TDS)
  TDS_scaled <- (gdata$dataset$TDS-TDS.mn)/(TDS.x-TDS.mn)
  cx_rsquare <- TDS_scaled + 1
  cy_rsquare <- gdata$anions[,1]
  # Upper ternary plot
  piper.df <- piper_transform(gdata) #transform_piper_data(gdata)
  cx_usquare <- piper.df$an_x-1.2
  cy_usquare <- 1 + piper.df$an_y
  # Left ternary plot
  cx_lsquare <- -piper.df$cat_y
  cy_lsquare <- piper.df$cat_x
  #
  results <- list(cx_mnsquare = cx_mnsquare, cy_mnsquare = cy_mnsquare,
                  cx_lsquare = cx_lsquare, cy_lsquare = cy_lsquare,
                  cx_lwsquare = cx_lwsquare, cy_lwsquare = cy_lwsquare,
                  cx_usquare = cx_usquare, cy_usquare = cy_usquare,
                  cx_rsquare = cx_rsquare, cy_rsquare = cy_rsquare)
  results <- as.data.frame(results)
  return(results)
}
#' @title
#' durov_transform
#' @description
#' Calculate the coordinates in the Durov plot.
#' @param Ca A numeric vector with the Ca concentration
#' @param Mg A numeric vector with the Mg concentration
#' @param Na A numeric vector with the Na concentration
#' @param K A numeric vector with the K concentration
#' @param HCO3 A numeric vector with the HCO3 concentration
#' @param CO3 A numeric vector with the CO3 concentration
#' @param Cl A numeric vector with the Cl concentration
#' @param SO4 A numeric vector with the SO4 concentration
#' @param pH A numeric vector with the pH concentration
#' @param TDS A numeric vector with the TDS concentration
#' @return
#' A list with the following entries
#' \itemize{
#' \item cx_mnsquare, cy_mnsquare
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family Durov functions
#' @export durov_transform
durov_transform <- function(Ca, Mg, Na, K, HCO3, CO3, Cl, SO4, pH, TDS){
  #
  tri_data <- anion_cation_transform(Ca, Mg, Na, K, HCO3, CO3, Cl, SO4)
  # Main square: Ca vs HCO3+CO3
  cx_mnsquare <- tri_data$cations[,1]
  cy_mnsquare <- tri_data$anions[,1]
  # pH Coordinates: y in (-1,0)
  pH.mn <- min(pH)
  pH.mx <- max(pH)
  pH_scaled <- -1+(pH-pH.mn)/((pH.mx-pH.mn))
  cx_lwsquare <- tri_data$cations[,1]
  cy_lwsquare <- pH_scaled
  # TDS coordinates: x in (1,2) y in (0,1)
  TDS.mn <- min(TDS)
  TDS.x <- max(TDS)
  TDS_scaled <- (TDS-TDS.mn)/(TDS.x-TDS.mn)
  cx_rsquare <- TDS_scaled + 1
  cy_rsquare <- tri_data$anions[,1]
  # Upper ternary plot
  cx_usquare <- tri_data$an_x-1.2
  cy_usquare <- 1 + tri_data$an_y
  # Left ternary plot
  cx_lsquare <- -tri_data$cat_y
  cy_lsquare <- tri_data$cat_x
  #
  results <- list(cx_mnsquare = cx_mnsquare, cy_mnsquare = cy_mnsquare,
                  cx_lsquare = cx_lsquare, cy_lsquare = cy_lsquare,
                  cx_lwsquare = cx_lwsquare, cy_lwsquare = cy_lwsquare,
                  cx_usquare = cx_usquare, cy_usquare = cy_usquare,
                  cx_rsquare = cx_rsquare, cy_rsquare = cy_rsquare)
}
#' @title
#' ggplot_durov
#' @description
#' Function
#' @return
#' A ggplot with the Durov diagram
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family durov functions
#' @importFrom ggplot2 ggplot geom_segment geom_text theme_bw theme aes
#' @export ggplot_durov
ggplot_durov <- function(){
  #
  x1 <- NULL
  x2 <- NULL
  y1 <- NULL
  y2 <- NULL
  grid1p1 <- data.frame(x1 = c(.20,.40,.60,.80), x2= c(.10,.20,.30,.40),
                        y1 = c(1,1,1,1), y2 = c(1.173206,1.346412,1.519618,1.692824))
  grid1p2 <- data.frame(x1 = c(.20,.40,.60,.80), x2= c(.60,.70,.80,.90),
                        y1 = c(1,1,1,1), y2 = c(1.692824,1.519618,1.346412,1.173206))
  grid1p3 <- data.frame(x1 = c(.10,.20,.30,.40), x2= c(.90,.80,.70,.60),
                        y1 = c(1.173206, 1.346412, 1.519618, 1.692824),
                        y2 = c(1.173206, 1.346412, 1.519618, 1.692824))
  #
  grid2p1 <- data.frame(y1 = c(.20,.40,.60,.80), y2= c(.10,.20,.30,.40),
                        x1 = c(0,0,0,0), x2 = c(-.173206,-.346412,-.519618,-.692824))
  grid2p2 <- data.frame(y1 = c(.20,.40,.60,.80), y2= c(.60,.70,.80,.90),
                        x1 = c(0,0,0,0), x2 = c(-.692824, -.519618,-.346412,-.173206))
  grid2p3 <- data.frame(y1 = c(.10,.20,.30,.40), y2= c(.90,.80,.70,.60),
                        x1 = c(-.173206, -.346412, -.519618, -.692824),
                        x2 = c(-.173206, -.346412, -.519618, -.692824))

  #
  gridmp1 <- data.frame(x1 = c(0, 0, 0, 0), x2= c(1,1,1,1),
                        y1 = c(.2,.4,.6,.8), y2 = c(.2,.4,.6,.8))
  gridmp2 <- data.frame(x1 = c(.2,.4,.6,.8), x2= c(.2,.4,.6,.8),
                        y1 = c(0, 0, 0, 0), y2 = c(1,1,1,1))
  #
  gridlp1 <- data.frame(x1 = c(0, 0, 0, 0), x2= c(1,1,1,1),
                        y1 = c(-.2,-.4,-.6,-.8), y2 = c(-.2,-.4,-.6,-.8))
  gridlp2 <- data.frame(x1 = c(.2,.4,.6,.8), x2= c(.2,.4,.6,.8),
                        y1 = c(0, 0, 0, 0), y2 = c(-1,-1,-1,-1))
  #
  gridrp1 <- data.frame(x1 = c(1, 1, 1, 1), x2= c(2,2,2,2),
                        y1 = c(.2,.4,.6,.8), y2 = c(.2,.4,.6,.8))
  gridrp2 <- data.frame(x1 = c(1.2,1.4,1.6,1.8), x2= c(1.2,1.4,1.6,1.8),
                        y1 = c(0, 0, 0, 0), y2 = c(1,1,1,1))
  #
  p <- ggplot() + geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0)) +
    geom_segment(aes(x = 1, y = 0, xend = 1, yend = 1)) +
    geom_segment(aes(x = 1, y = 1, xend = 0, yend = 1)) +
    geom_segment(aes(x = 0, y = 1, xend = 0, yend = 0)) +
    # TDS box
    geom_segment(aes(x = 1, y = 0, xend = 2, yend = 0)) +
    geom_segment(aes(x = 2, y = 0, xend = 2, yend = 1)) +
    geom_segment(aes(x = 2, y = 1, xend = 1, yend = 1)) +
    # pH box
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = -1)) +
    geom_segment(aes(x = 0, y = -1, xend = 1, yend = -1)) +
    geom_segment(aes(x = 1, y = -1, xend = 1, yend = 0)) +
    # Upper ternary plot
    geom_segment(aes(x = 0, y = 1, xend = 0.5, yend = (1.+ sqrt(3)/2) )) +
    geom_segment(aes(x = 0.5, y = (1.0 + sqrt(3)/2), xend = 1, yend = 1)) +
    # Left ternary plot
    geom_segment(aes(x = 0, y = 1, xend = -sqrt(3)/2, yend = 0.5)) +
    geom_segment(aes(x = -sqrt(3)/2, y = 0.5, xend = 0, yend = 0)) +
    # Grid lines of the upper ternary plot
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=grid1p1,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=grid1p2,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=grid1p3,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    # Grid lines of the left ternary plot
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=grid2p1,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=grid2p2,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=grid2p3,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    # Grid lines of the central square
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=gridmp1,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=gridmp2,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    # Grid lines of the lower square
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=gridlp1,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=gridlp2,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    # Grid lines of the right square
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=gridrp1,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data=gridrp2,
                 linetype = "dashed", size = 0.25, colour = "grey50") +
    # labels for upper ternary plot
    geom_text(aes(c(.20, .40, .60, .80), c(1-.05, 1-.05, 1-.05, 1-.05),
                  label=c(80, 60, 40, 20)), size=3) +
    geom_text(aes(c(.35, .25, .15, .05), grid1p2$y2, label=c(80, 60, 40, 20)), size=3) +
    geom_text(aes(c(.95, .85, .75, .65), grid1p3$y2, label=c(80, 60, 40, 20)), size=3) +
    # labels for left ternary plot
    geom_text(aes(c(-.173206,-.346412,-.519618,-.692824), c(.10,.20,.30,.40)-0.05,
                  label=c(20, 40, 60, 80)), size = 3) +
    geom_text(aes(rep(0.05, 4), c(.20,.40,.60,.80), label = c(80, 60, 40, 20)),
              size = 3) +
    geom_text(aes(x = c(-.692824, -.519618,-.346412,-.173206),
                  y = c(.60,.70,.80,.90) + 0.05,
                  label = c(20, 40, 60, 80)), size = 3) +
    #
    geom_text(aes(x =-0.075, y = 1, label = "Ca"), size = 3) +
    geom_text(aes(x =-0.05, y = 1.05, label = "SO4"), size = 3) +
    geom_text(aes(x = 1.15, y = 1.05, label = "HCO3 + CO3"), size = 3) +
    geom_text(aes(x = 0.5, y = 1.9, label = "Cl"), size = 3) +
    geom_text(aes(x = -.1, y = -0.05, label = "Na + K"), size = 3) +
    geom_text(aes(x = -0.92, y = 0.5, label = "Mg"), size = 3) +
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
#' plot_durov
#' @description
#' Function to create the Durov plot
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
#' @family durov functions
#' @importFrom ggplot2 ggplot geom_point coord_fixed scale_color_gradientn scale_colour_gradientn scale_color_discrete aes_string
#' @importFrom grDevices rainbow
#' @export
plot_durov <- function(x, measure = c('conc', 'meql'),
                       vars = NULL, color = NULL,
                       Size = NULL){
  gdata <- x
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  durov.df <- durov_transform1(gdata)
  p <- ggplot_durov()
  cx_mnsquare <- NULL
  cy_mnsquare <- NULL
  cx_lwsquare <- NULL
  cy_lwsquare <- NULL
  cx_rsquare <- NULL
  cy_rsquare <- NULL
  cx_usquare <- NULL
  cy_usquare <- NULL
  cx_lsquare <- NULL
  cy_lsquare <- NULL
  if(is.null(color)){
    if(is.null(Size)){
      p <- p +  geom_point(aes(x=cx_mnsquare, y = cy_mnsquare),
                           data = durov.df) +
        geom_point(aes(x = cx_lwsquare, y = cy_lwsquare), data = durov.df) +
        geom_point(aes(x = cx_rsquare, y = cy_rsquare), data = durov.df) +
        geom_point(aes(x = cx_usquare, y = cy_usquare), data = durov.df) +
        geom_point(aes(x = cx_lsquare, y = cy_lsquare), data = durov.df) +
        coord_fixed()
    }
    else{
      if(class(Size) == "numeric"){
        p <- p +  geom_point(aes(x=cx_mnsquare, y = cy_mnsquare),
                             data = durov.df, size = Size) +
          geom_point(aes(x = cx_lwsquare, y = cy_lwsquare), data = durov.df,
                     size = Size) +
          geom_point(aes(x = cx_rsquare, y = cy_rsquare), data = durov.df,
                     size = Size) +
          geom_point(aes(x = cx_usquare, y = cy_usquare), data = durov.df,
                     size = Size) +
          geom_point(aes(x = cx_lsquare, y = cy_lsquare), data = durov.df,
                     size = Size) +
          coord_fixed()
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        durov.df[Size] <- tmp[,1]
        p <- p +  geom_point(aes_string(x = "cx_mnsquare",
                                        y = "cy_mnsquare",
                                        size = Size),
                             data = durov.df) +
          geom_point(aes_string(x = "cx_lwsquare", y = "cy_lwsquare",
                                size = Size),
                     data = durov.df) +
          geom_point(aes_string(x = "cx_rsquare", y = "cy_rsquare",
                                size = Size),
                     data = durov.df) +
          geom_point(aes_string(x = "cx_usquare", y = "cy_usquare",
                                size = Size),
                     data = durov.df) +
          geom_point(aes_string(x = "cx_lsquare", y = "cy_lsquare",
                                size = Size),
                     data = durov.df) +
          coord_fixed()
      }
    }
  }
  else{
    tmp <- gdata$dataset[color]
    durov.df[color] <- tmp[,1]
    if(is.null(Size)){
      p <- p +  geom_point(aes_string(x = "cx_mnsquare",
                                      y = "cy_mnsquare",
                                      color = color),
                           data = durov.df,
                           size = 3) +
        geom_point(aes_string(x = "cx_lwsquare", y = "cy_lwsquare", color = color),
                   data = durov.df, size =3) +
        geom_point(aes_string(x = "cx_rsquare", y = "cy_rsquare", color = color),
                   data = durov.df, size = 3) +
        geom_point(aes_string(x = "cx_usquare", y = "cy_usquare", color = color),
                   data = durov.df, size = 3) +
        geom_point(aes_string(x = "cx_lsquare", y = "cy_lsquare", color = color),
                   data =durov.df, size = 3) +
        coord_fixed()
    }
    else{
      if(class(Size) == "numeric"){
        p <- p +  geom_point(aes_string(x = "cx_mnsquare",
                                        y = "cy_mnsquare",
                                        color = color),
                             data = durov.df,
                             size = Size) +
          geom_point(aes(x = cx_lwsquare, y = cy_lwsquare), data = durov.df,
                     size = Size) +
          geom_point(aes(x = cx_rsquare, y = cy_rsquare), data = durov.df,
                     size = Size) +
          geom_point(aes(x = cx_usquare, y = cy_usquare), data = durov.df,
                     size = Size) +
          geom_point(aes(x = cx_lsquare, y = cy_lsquare), data = durov.df,
                     size = Size) +
          coord_fixed()
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        durov.df[Size] <- tmp[,1]
        p <- p +  geom_point(aes_string(x = "cx_mnsquare",
                                        y = "cy_mnsquare",
                                        color = color,
                                        size = Size),
                             data = durov.df) +
          geom_point(aes_string(x = "cx_lwsquare", y = "cy_lwsquare",
                                color = color, size = Size),
                     data = durov.df) +
          geom_point(aes_string(x = "cx_rsquare", y = "cy_rsquare",
                                color = color, size = Size),
                     data = durov.df) +
          geom_point(aes_string(x = "cx_usquare", y = "cy_usquare",
                                color = color, size = Size),
                     data = durov.df) +
          geom_point(aes_string(x = "cx_lsquare", y = "cy_lsquare",
                                color = color, size = Size),
                     data =durov.df) +
          coord_fixed()
      }
    }
  }
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

