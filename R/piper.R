#' @section Piper Diagram:
#'
#' The functions defined to create a Piper diagram are:
#' anion_cation_tranform, transform_piper_data, ggplot_piper
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' anion_cation_transform
#' @description
#' Function to calculate the coordinates of anions and cations used in the Piper diagram.
#' @param Ca A numeric vector with the Ca concentration
#' @param Mg A numeric vector with the Mg concentration
#' @param Na A numeric vector with the Na concentration
#' @param K A numeric vector with the K concentration
#' @param HCO3 A numeric vector with the HCO3 concentration
#' @param CO3 A numeric vector with the CO3 concentration
#' @param Cl A numeric vector with the Cl concentration
#' @param SO4 A numeric vector with the SO4 concentration
#' @param name A numeric vector with the Ca concentration
#' @return
#' A list with the following entries:
#' \itemize{
#' \item cat_x,cat_y: numeric vectors with the x and y coordinates used to create the
#' ternary diagram of the cations.
#' \item an_x, an_ynumeric vectors with the x and y coordinates used to create the ternary
#'  diagram of the anions.
#' \item d_x, d_y
#' \item  meqL: A numeric vector with the meqL of the input concentrations.
#' \item sumanions: Sum of the meqL of the anions
#' \item anions: A matrix with the values of the normalized meqL of the anions.
#' \item cations: A matrix with the values of the normalized meqL of the anions.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family piper functions
#' @export anion_cation_transform
anion_cation_transform <- function(Ca, Mg, Na, K, HCO3, CO3, Cl, SO4, name = NULL){
  h <- 0.5*tan(pi/3)
  offset <- 0.1
  offsety <- offset*tan(pi/3)
  gmol <- c(40.078, 24.305, 22.989768, 39.0983, 61.01714, 60.0092, 35.4527, 96.0636)
  eqmol <- c(2., 2., 1., 1., 1., 2., 1., 2.)
  ndat <- length(Ca)
  meqL <- matrix(0.0, nrow = ndat, ncol = 8)
  meqL[,1] <- Ca*eqmol[1]/gmol[1]
  meqL[,2] <- Mg*eqmol[2]/gmol[2]
  meqL[,3] <- Na*eqmol[3]/gmol[3]
  meqL[,4] <- K*eqmol[4]/gmol[4]
  meqL[,5] <- HCO3*eqmol[5]/gmol[5]
  meqL[,6] <- CO3*eqmol[6]/gmol[6]
  meqL[,7] <- Cl*eqmol[7]/gmol[7]
  meqL[,8] <- SO4*eqmol[8]/gmol[8]
  sumcations <- rowSums(meqL[,1:4])
  sumanions <- rowSums(meqL[,5:8])
  # Coordinates of cations
  cations <- matrix(0.0, nrow = ndat, ncol = 3)
  cations[,1] <- meqL[,1]/sumcations
  cations[,2] <- meqL[,2]/sumcations
  cations[,3] <- (meqL[,3]+meqL[,4])/sumcations
  cat_x <- 0.5*( 2*cations[,3] + cations[,2] )
  cat_y <- h*cations[,2]
  # Coordinates of anions
  anions <- matrix(0.0, nrow = ndat, ncol = 3)
  anions[,1] <- (meqL[,5]+meqL[,6])/sumanions
  anions[,3] <- meqL[,7]/sumanions
  anions[,2] <- meqL[,8]/sumanions
  an_x <- 1 + 2*offset + 0.5*( 2*anions[,3] + anions[,2] )
  an_y  <- h*anions[,2]
  # Coordinates inside diamons
  d_x <- an_y/(4*h) + 0.5*an_x - cat_y/(4*h) + 0.5*cat_x
  d_y <- 0.5*an_y + h*an_x + 0.5*cat_y - h*cat_x
  #
  results <- list(cat_x = cat_x, cat_y = cat_y, an_x = an_x, an_y = an_y,
                  d_x = d_x, d_y = d_y, meqL = meqL, sumanions = sumanions,
                  anions = anions, cations = cations)
  return(results)
}
#' @title
#' transform_piper_data
#' @description
#' Function to transform the concentrations of Mg, Ca, Cl and SO4 to local coordinates of the Piper diagram.
#' @param Mg A numeric vector with the concentrations of Mg
#' @param Ca A numeric vector with the concentrations of Ca
#' @param Cl A numeric vector with the concentrations of Cl
#' @param SO4 A numeric vector with the concentrations of SO4
#' @param name name?
#' @return
#' A data.frame with the following variables
#' \itemize{
#' \item observation:
#' \item x:
#' \item y:
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family piper functions
#' @export transform_piper_data
transform_piper_data <- function(Mg, Ca, Cl, SO4, name = NULL){
  if(is.null(name)){
    name = rep(1:length(Mg),3)
  } else {
    name = rep(name,3)
  }
  y1 <- Mg * 0.86603
  x1 <- 100*(1-(Ca/100) - (Mg/200))
  y2 <- SO4 * 0.86603
  x2 <-120+(100*Cl/100 + 0.5 * 100*SO4/100)
  # Function to define the coordinates of the point
  new_point <- function(x1, x2, y1, y2, grad=1.73206){
    b1 <- y1-(grad*x1)
    b2 <- y2-(-grad*x2)
    M <- matrix(c(grad, -grad, -1,-1), ncol=2)
    intercepts <- as.matrix(c(b1,b2))
    t_mat <- -solve(M) %*% intercepts
    data.frame(x=t_mat[1,1], y=t_mat[2,1])
  }
  np_list <- lapply(1:length(x1), function(i) new_point(x1[i], x2[i], y1[i], y2[i]))
  npoints <- do.call("rbind", np_list)
  data.frame(observation = name, x = c(x1, x2, npoints$x), y=c(y = y1, y2, npoints$y))
}
#' @title
#' piper_transform
#' @description
#' Function to calculate the coordinates of set of samples in a piper diagram.
#' @param gdata A geochemical_dataset
#' @return
#' This function returns a data.frame with the following variables:
#' \itemize{
#' \item cat_x,cat_y: coordinates of the cations in the ternary diagram (scaled between 0 and 1).
#' \item an_x,an_y: coordiantes of the anions in the ternary diagram (scaled between 0 and 1).
#' \item d_x,d_y: coordinates of the samples in the central diamond (scaled between 0 and 1).
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family piper functions
#' @export
piper_transform <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset object is required as input')
  }
  h <- 0.5*tan(pi/3)
  offset <- 0.1
  #
  cat_x <- 0.5*( 2*gdata$cations[,3] + gdata$cations[,2] )
  cat_y <- h*gdata$cations[,2]
  #
  an_x <- 1 + 2*offset + 0.5*( 2*gdata$anions[,3] + gdata$anions[,2] )
  an_y  <- h*gdata$anions[,2]
  # Coordinates inside diamons
  d_x <- an_y/(4*h) + 0.5*an_x - cat_y/(4*h) + 0.5*cat_x
  d_y <- 0.5*an_y + h*an_x + 0.5*cat_y - h*cat_x
  #
  results <- list(cat_x = cat_x, cat_y = cat_y, an_x = an_x, an_y = an_y,
                  d_x = d_x, d_y = d_y)
  return(as.data.frame(results))
}

#' @title
#' ggplot_piper
#' @description
#' Function to create the ternary diagrams that compose the Piper diagram
#' @return
#' A ggplot2 object
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family piper functions
#' @importFrom ggplot2 ggplot geom_segment geom_text theme_bw theme element_blank coord_equal
#' @export ggplot_piper
ggplot_piper <- function() {
  x1 <- NULL
  x2 <- NULL
  y1 <- NULL
  y2 <- NULL
  #
  p <- ggplot() +
    ## left hand ternary plot
    geom_segment(aes(x=0,y=0, xend=100, yend=0)) +
    geom_segment(aes(x=0,y=0, xend=50, yend=86.603)) +
    geom_segment(aes(x=50,y=86.603, xend=100, yend=0)) +
    ## right hand ternary plot
    geom_segment(aes(x=120,y=0, xend=220, yend=0)) +
    geom_segment(aes(x=120,y=0, xend=170, yend=86.603)) +
    geom_segment(aes(x=170,y=86.603, xend=220, yend=0)) +
    ## Upper diamond
    geom_segment(aes(x=110,y=190.5266, xend=60, yend=103.9236)) +
    geom_segment(aes(x=110,y=190.5266, xend=160, yend=103.9236)) +
    geom_segment(aes(x=110,y=17.3206, xend=160, yend=103.9236)) +
    geom_segment(aes(x=110,y=17.3206, xend=60, yend=103.9236)) +
    #
    coord_equal(ratio = 1) +
    #
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  return(p)
}
#' @title
#' add_grid_lines_piper
#' @description
#' Function to add the grid lines to a Piper plot
#' @param color A character string specifying the color of the grid lines in 
#' Piper plot
#' @param Size A numeric value specifying the size of the grid labels in a 
#' Piper plot
#' @return
#' This function returns a list with the geom components used to create 
#' a grid in a Piper diagram
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family piper functions
#' @importFrom ggplot2 geom_segment geom_text
#' @export 
add_grid_lines_piper <- function(color = NULL, Size = NULL){
  color1 <- "grey50"
  Size1 <- 3
  if(!missing(color)){
    color1 <- color
  }
  if(!missing(Size)){
    Size1 <- Size
  }
  #
  x1 <- NULL
  x2 <- NULL
  y1 <- NULL
  y2 <- NULL
  #
  grid1p1 <- data.frame(x1 = c(20,40,60,80), x2= c(10,20,30,40), 
                        y1 = c(0,0,0,0), 
                        y2 = c(17.3206,34.6412,51.9618, 69.2824))
  grid1p2 <- data.frame(x1 = c(20,40,60,80), x2= c(60,70,80,90),
                        y1 = c(0,0,0,0), 
                        y2 = c(69.2824, 51.9618,34.6412,17.3206))
  grid1p3 <- data.frame(x1 = c(10,20,30,40), x2= c(90,80,70,60), 
                        y1 = c(17.3206,34.6412,51.9618, 69.2824), 
                        y2 = c(17.3206,34.6412,51.9618, 69.2824))
  #
  grid2p1 <- grid1p1
  grid2p1$x1 <- grid2p1$x1+120
  grid2p1$x2 <- grid2p1$x2+120
  grid2p2 <- grid1p2
  grid2p2$x1 <- grid2p2$x1+120
  grid2p2$x2 <- grid2p2$x2+120
  grid2p3 <- grid1p3
  grid2p3$x1 <- grid2p3$x1+120
  grid2p3$x2 <- grid2p3$x2+120
  # Gridlines diamons
  grid3p1 <- data.frame(x1 = c(100, 90, 80, 70),
                        y1 = c(34.6412, 51.9618, 69.2824, 86.603),
                        x2 = c(150, 140, 130, 120),
                        y2 = c(121.2442,138.5648,155.8854,173.2060))
  grid3p2 <- data.frame(x1 = c(70, 80, 90, 100), 
                        y1 = c(121.2442,138.5648,155.8854,173.2060),
                        x2 = c(120, 130, 140, 150),
                        y2 = c(34.6412, 51.9618, 69.2824, 86.603))
  
  ## Add grid lines to the plots
  res <- list()
  res[[1]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid1p1, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[2]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid1p2, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[3]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid1p3, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[4]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid2p1, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[5]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid2p2, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[6]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid2p3, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[7]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid3p1, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[8]] <- geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), 
                           data=grid3p2, linetype = "dashed", 
                           size = 0.25, colour = color1)
  res[[9]] <- geom_text(aes(c(155,145,135,125), grid2p2$y2, 
                            label=c(20, 40, 60, 80)), size = Size1, 
                        colour = color1)
  res[[10]] <- geom_text(aes(c(215,205,195,185),grid2p3$y2, 
                  label=c(20, 40, 60, 80)), size = Size1, 
                  colour = color1)
  res[[11]] <- geom_text(aes(c(140,160,180,200), c(-5,-5,-5,-5), 
                  label=c(20, 40, 60, 80)), size = Size1, 
                  colour = color1)
  res[[12]] <- geom_text(aes(grid3p1$x1-5,grid3p1$y1, 
                  label=c(80, 60, 40, 20)), size = Size1, 
                  colour = color1)
  res[[13]] <- geom_text(aes(grid3p1$x2+5,grid3p1$y2, 
                  label=c(20, 40, 60, 80)), size = Size1, 
                  colour = color1)
  res[[14]] <- geom_text(aes(grid3p2$x1-5,grid3p2$y1, 
                  label=c(20, 40, 60, 80)), size = Size1, 
                  colour = color1)
  res[[15]] <- geom_text(aes(grid3p2$x2+5,grid3p2$y2, 
                  label=c(80, 60, 40, 20)), size=Size1, 
                  colour = color1)
  res[[16]] <- geom_text(aes(c(20,40,60,80),c(-5,-5,-5,-5),
                              label=c(80, 60, 40, 20)), size = Size1,
                          colour = color1)
  res[[17]] <- geom_text(aes(c(35,25,15,5), grid1p2$y2,
                              label=c(80, 60, 40, 20)), size = Size1,
                          colour = color1)
  res[[18]] <- geom_text(aes(c(95,85,75,65), grid1p3$y2,
                              label=c(80, 60, 40, 20)), size = Size1,
                          colour = color1)
  
  #
  return(res)
}
#' @title 
#' add_labels_piper
#' @description 
#' Function to add the labels to a ternary diagram
#' @param labels A character vector with the labels to be added to the 
#' Piper diagram 
#' @param middle A logical flag indicating if the labels must be located on the 
#' middle part of the axis. 
#' @param color A character string specifying the color of the labels
#' @param Size A numeric value specifying the size of the text label
#' @return 
#' This function returns a geom_text to be added to the Piper diagram
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family piper functions
#' @importFrom ggplot2 geom_text
#' @export
add_labels_piper <- function(labels, middle = TRUE,  color = NULL, 
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
  res[[1]] <- geom_text(aes(17, 50, label="Mg2"), angle = 60, size = Size1, 
                        colour = color1, parse = TRUE)
  res[[2]] <- geom_text(aes(82.5, 50, label="Na + K"), angle = -60, 
                        size=Size1, colour = color1, parse = TRUE)
  res[[3]] <- geom_text(aes(50, -10, label="Ca2"), size = Size1, 
                        colour = color1, parse = TRUE)
  #
  res[[4]] <- geom_text(aes(170, -10, label="Cl-phantom()"), size = Size1, 
                        colour = color1, parse = TRUE)
  res[[5]] <- geom_text(aes(205,50, label="SO4"), angle=-60, size = Size1, 
                        colour = color1, parse = TRUE)
  res[[6]] <- geom_text(aes(137.5, 50, label = "Alkalinity~as~HCO3"), 
                        angle=60, size = Size1, colour = color1, 
                        parse = TRUE)
  res[[7]] <- geom_text(aes(72.5, 150, label = "SO4~+~Cl-phantom()"), 
                        angle=60, size = Size1, colour = color1, 
                        parse = TRUE)
  res[[8]] <- geom_text(aes(147.5,150, label="Ca2~+~Mg2"), angle=-60, 
              size = Size1, colour = color1, parse = TRUE)
  #  
  return(res)  
}
#' @title
#' plot_piper
#' @description
#' Function to create the Piper plot
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
#' @family piper functions
#' @importFrom ggplot2 ggplot geom_point coord_fixed scale_color_gradientn scale_colour_gradientn scale_color_discrete aes_string
#' @importFrom grDevices rainbow
#' @export
plot_piper <- function(x, measure = c('conc', 'meql'),
                       vars = NULL, color = NULL,
                       Size = NULL, 
                       additional.args = NULL){
  gdata <- x
  conc_ions <- colnames(gdata$dataset)
  meql_ions <- c("Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4")
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset is required as input')
  }
  #  if(missing(vars) | is.null(vars)){
  #    stop('ERROR: Variable names must be specified to create a ternary diagram')
  #  }
  p <- ggplot_piper() + add_grid_lines_piper() + add_labels_piper()
  piper.df <- piper_transform(gdata)
  cat_x <- NULL
  cat_y <- NULL
  an_x <- NULL
  an_y <- NULL
  d_x <- NULL
  d_y <- NULL
  if(is.null(color)){
    if(is.null(Size)){
      p <- p + geom_point(aes(x = 100*cat_x, y = 100*cat_y), data = piper.df) +
        geom_point(aes(x = 100*an_x, y = 100*an_y), data= piper.df) +
        geom_point(aes(x = 100*d_x, y = 100*d_y), data = piper.df)
    }
    else{
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes(x = 100*cat_x, y = 100*cat_y), data = piper.df,
                            size = Size) +
          geom_point(aes(x = 100*an_x, y = 100*an_y), data= piper.df,
                     size = Size) +
          geom_point(aes(x = 100*d_x, y = 100*d_y), data = piper.df,
                     size = Size)
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        piper.df[Size] <- tmp[,1]
        p <- p + geom_point(aes(x = 100*cat_x, y = 100*cat_y, size = Size), data = piper.df) +
          geom_point(aes(x = 100*an_x, y = 100*an_y, size = Size), data = piper.df) +
          geom_point(aes(x = 100*d_x, y = 100*d_y, size = Size), data = piper.df)
      }
    }
  }
  else{
    tmp <- gdata$dataset[color]
    piper.df[color] <- tmp[,1]
    if(is.null(Size)){
      p <- p + geom_point(aes_string(x = "100*cat_x", y = "100*cat_y", color = color),
                          data = piper.df, size = 3) +
        geom_point(aes_string(x = "100*an_x", y = "100*an_y", color = color), data= piper.df, size = 3) +
        geom_point(aes_string(x = "100*d_x", y = "100*d_y", color = color), data = piper.df, size = 3)
    }
    else{
      if(class(Size) == "numeric"){
        p <- p + geom_point(aes_string(x = "100*cat_x", y = "100*cat_y",
                                       color = color),
                            data = piper.df, size = Size) +
          geom_point(aes_string(x = "100*an_x", y = "100*an_y", color = color),
                     data= piper.df, size = Size) +
          geom_point(aes_string(x = "100*d_x", y = "100*d_y", color = color),
                     data= piper.df, size = Size)
      }
      else if(class(Size) == "character"){
        tmp <- gdata$dataset[Size]
        piper.df[Size] <- tmp[,1]
        p <- p + geom_point(aes_string(x = "100*cat_x", y = "100*cat_y",
                                       color = color, size = Size),
                            data = piper.df) +
          geom_point(aes_string(x = "100*an_x", y = "100*an_y",
                                color = color, size = Size),
                     data = piper.df) +
          geom_point(aes_string(x = "100*d_x", y = "100*d_y",
                                color = color, size = Size),
                     data = piper.df)
      }
    }
  }
  #
  p <- p + coord_fixed()
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
#' @title
#' piper_classification
#' @description
#' Function to classify the groundwater samples according to the Piper diagram
#' @param d_x A numeric vector with the x coordinates of the samples inside the diamond of the Piper diagram.
#' @param d_y A numeric vector with the x coordinates of the samples inside the diamond of the Piper diagram.
#' @return
#' This function returns a character vector with the corresponding classification
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family piper functions
#' @importFrom sp point.in.polygon
#' @export
piper_classification <- function(d_x, d_y){
  # Upper triangle: calcium chloride
  pol.CaCl.x <- c(110,85,135,110)
  pol.CaCl.y <- c(190.5266,147.2251,147.2251,190.5266)
  # Mixed type 1
  pol.mixed1.x <- c(85,135,60,85)
  pol.mixed1.y <- c(147.2251,147.2251,103.9236,147.2251)
  # Mixed type 2
  pol.mixed2.x <- c(110,85,135,110)
  pol.mixed2.y <- c(103.9236,60.62,60.62,103.9236)
  # Magnesium-Bicarbonate
  pol.MgHCO3.x <- c(60,85,110,85,60)
  pol.MgHCO3.y <- c(103.9236,147.2251,103.9236,60.62,103.9236)
  # NaCl
  pol.NaCl.x <- c(110,135,160,135,110)
  pol.NaCl.y <- c(103.9236,147.2251,103.9236,60.62,103.9236)
  # NaHCO3
  pol.NaHCO3.x <- c(85,135,110,85)
  pol.NaHCO3.y <- c(60.62,60.62,17.3206,60.62)
  #
  classification <- c("Calcium.Chloride", "Mixed1", "Mixed2", "Magnesium.Bicarbonate",
                      "Sodium.Chloride", "Sodium.Bicarbonate")
  #
  ndat <- length(d_x)
  results <- matrix(0, nrow = ndat, ncol = 6)
  for(i in 1:ndat){
    results[i,1] <- point.in.polygon(d_x[i], d_y[i], pol.CaCl.x, pol.CaCl.y)
    results[i,2] <- point.in.polygon(d_x[i], d_y[i], pol.mixed1.x, pol.mixed1.y)
    results[i,3] <- point.in.polygon(d_x[i], d_y[i], pol.mixed2.x, pol.mixed2.y)
    results[i,4] <- point.in.polygon(d_x[i], d_y[i], pol.MgHCO3.x, pol.MgHCO3.x)
    results[i,5] <- point.in.polygon(d_x[i], d_y[i], pol.NaCl.x, pol.NaCl.y)
    results[i,6] <- point.in.polygon(d_x[i], d_y[i], pol.NaHCO3.x, pol.NaHCO3.y)
  }
  #
  pos.classification <- apply(results,1,which.max)
  #
  results.class <- vector("character", length = ndat)
  for(i in 1:ndat){
    results.class[i] <- classification[pos.classification[i]]
  }
  return(results.class)
}
