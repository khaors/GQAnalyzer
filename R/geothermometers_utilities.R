#' @section Geothermometers:
#'
#' The functions defined to calculate the temperatures
#'
#'
#' @docType package
#' @name GQAnalyzer
NULL
#' @title
#' calculate_geothermometers
#' @description
#' Function to calculate several geothermometers
#' @param gdata A geochemical_dataset object
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
calculate_geothermometers <- function(gdata){
  if(class(gdata) != "geochemical_dataset"){
    stop('ERROR: A geochemical_dataset object is required as input')
  }
  cations_var <- c("Ca", "Mg", "Na", "K")
  anions_var <- c("HCO3", "CO3", "Cl", "SO4")
  #
  datanames <- names(gdata$dataset)

  Temp <- gdata$dataset$Temp
  Na <- gdata$dataset$Na
  Mg <- gdata$dataset$Mg
  K <- gdata$dataset$K
  SiO2 <- gdata$dataset$SiO2
  #
  # Silica geothermometers
  # t < 250
  # SiO2 = concentration of Silica in mg/kg
  #
  # For t > 330 Fournier and Potter (1982)
  #

  #
  # Na/K Geothermometers
  #


}
#' @title
#' silica.geothermometers
#' @description
#' This function calculates several SiO2 geothermometers
#' @param SiO2 A numeric vector with the SiO2 concentration
#' @param Temp A numeric vector with the temperature of the samples
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
silica.geothermometers <- function(SiO2, Temp){
  pos.low <- Temp <= 250
  pos.high <- Temp > 250
  ndat <- length(SiO2)
  results <- matrix(0.0, nrow = ndat, ncol = 6)
  Qz.no.steam.loss <- (1309/(5.15-log(SiO2[pos.low])))-273
  Qz.max.steam <- (1522/(5.75-log(SiO2[pos.low])))-273
  Chalcedony <- (1032/(4.69-log(SiO2[pos.low])))-273
  alpha.Cristobalite <- (1000/(4.88-log(SiO2[pos.low])))-273
  beta.Cristobalite <- (781/(4.51-log(SiO2[pos.low])))-273
  amorphous.silica <- (731/(4.52-log(SiO2[pos.low])))-273
  #
  results[pos.low,1] <- Qz.no.steam.loss
  results[pos.low,2] <- Qz.max.steam
  results[pos.low,3] <- Chalcedony
  results[pos.low,4] <- alpha.Cristobalite
  results[pos.low,5] <- beta.Cristobalite
  results[pos.low,6] <- amorphous.silica
  for(i in 1:6){
    results[pos.high,i] <- NA
  }
  results.df <- as.data.frame(results)
  names(results.df) <- c("Qz.no.steam.loss", "Qz.max.steam", "Chalcedony",
                         "alpha.Cristobalite", "beta.Cristobalite",
                         "amorphous.silica")
  return(results.df)
}
#' @title
#' Fournier.Potter.geothermometer
#' @description
#' This function calculates the Fournier-Potter SiO2 geothermometer
#' @param SiO2 A numeric vector with the SiO2 concentration
#' @param Temp A numeric vector with the temperature of the samples
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
Fournier.Potter.geothermometer <- function(SiO2, Temp){
  pos.high <- Temp > 330
  pos.low <- Temp <= 330
  ndat <- length(SiO2)
  results <- vector("numeric", length = ndat)
  k1 <- 4.2198e1
  k2 <- 2.8831e-1
  k3 <- 3.6686e-4
  k4 <- 3.1665e-7
  k5 <- 7.7034e1
  #
  term1 <- -k1+k2*SiO2[pos.high]
  term2 <- -k3*(SiO2[pos.high])**2
  term3 <- k4*(SiO2[pos.high])**3
  term4 <- k5*log((SiO2[pos.high]))
  results[pos.high] <-term1 + term2 + term3 + term4
  results[pos.low] <- NA
  return(results)
}
#' @title
#' Na.K.geothermometers
#' @description
#' This function calculates several Na-K geothermometers
#' @param Na A numerical vector with Na concentrations
#' @param K A numerical vector with K concentrations
#' @param Temp A numerical vector with temperature concentrations
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
Na.K.geothermometers <- function(Na, K, Temp){
  pars <- matrix(0.0, nrow = 7, ncol = 2)
  pars[1,] <- c(855.6, 0.8573) #Truesdell
  pars[2,] <- c(883, 0.780)#Tonani
  pars[3,] <- c(933, 0.993)#Arnosson1
  pars[4,] <- c(1319, 1.699)#Arnosson2
  pars[5,] <- c(1217, 1.483)#Fournier
  pars[6,] <- c(1178, 1.470)#Nieva-Nieva, check second value
  pars[7,] <- c(1390, 1.750)#Giggenbach
  pars[8,] <- c(777, 0.7) #Fournier-Truesdell
  pars[9,] <- c(1289, 0.615)# Verma-Santoyo
  #
  ndat <- length(Na)
  results <- matrix(0.0, nrow = ndat, ncol = 9)
  for(i in 1:9){
    results[,i] <- vapply(Na/K,
                          FUN = function(x,par){res <- (par[1]/(log(x)+par[2]))-273},
                          par = pars[i,])
  }
  #Arnosson-3
  term1 <- log(Na/K)
  term2 <- term1**2
  term3 <- term1**3
  term4 <- term1**4
  results[10,] <- 733.6-770.551*term1+378.189*term2-95.753*term3+9.544*term4
  # Can
  results[11,] <- 1052/(1-exp(1.714*log(Na/K)+0.252))+76
  #DiazGonzales-Santoyo-Reyes1
  results[12,] <- (833/(log(Na/K)+0.894))-273.15
  #DiazGonzales-Santoyo-Reyes2
  results[13,] <- (833/(log(Na/K)+0.908))-273.15
  results.df <- as.data.frame(results)
  names(results.df) <- c("Truesdell", "Tonani",
                         "Arnosson1", "Arnosson2",
                         "Fournier",
                         "Nieva", "Giggengach", "Fournier.Truesdell",
                         "Verma.Santoyo", "Arnosson3", "Can",
                         "DiazGonzales-Santoyo-Reyes1",
                         "DiazGonzales-Santoyo-Reyes2")
  return(results.df)
}
#' @title
#' Na.K.Ca.geothermometer
#' @description
#' This function calculates tshe Na.K.Ca geothermometer
#' @param Na A numerical vector with Na concentrations
#' @param K A numerical vector with K concentrations
#' @param Ca A numerical vector with Ca concentrations
#' @param Temp A numerical vector with temperature concentrations
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
Na.K.Ca.geothermometer <- function(Na, K, Ca, Temp){
  check <- log((sqrt(Ca)/Na)+2.06)
  beta <- -100
  if(check > 0){
    beta <- 4/3
  }
  Temp.formation <- 1647/(log(Na/K)+beta*check+2.47)-273
  if(Temp.formation > 100){
    beta <- 1/3
    Temp.formation <- 1647/(log(Na/K)+beta*check+2.47)-273
  }
  return(Temp.formation)
}
#' @title
#' K.Mg.geothermometer
#' @description
#' This function calculates tshe Na.K.Ca geothermometer
#' @param K A numerical vector with K concentrations
#' @param Mg A numerical vector with Mg concentrations
#' @param Temp A numerical vector with temperature concentrations
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
K.Mg.geothermometer <- function(K, Mg, Temp){
  pars <- matrix(0.0, nrow = 3, ncol = 2)
  pars[1,] <- c(4410, 14)
  pars[2,] <- c(2330, 7.35)
  pars[3,] <- c(1077, 4.033)
  #
  ndat <- length(K)
  results <- matrix(0.0, nrow = ndat, ncol = 3)
  term1 <- log((K**2)/Mg)
  results[,1] <- (4410/(14-term1))-273.15
  results[,2] <- (2230/(7.35-term1))-273.15
  results[,3] <- (1077/(4.033+term1))-273.15
  #
  results.df <- as.data.frame(results)
  names(results.df) <- c("Giggenbach", "Fournier1", "Fournier2")
  return(results.df)
}
#' @title
#' Mg.Li.geothermometer
#' @description
#' This function calculates tshe Na.K.Ca geothermometer
#' @param Li A numerical vector with K concentrations
#' @param Mg A numerical vector with Mg concentrations
#' @param Temp A numerical vector with temperature concentrations
#' @return
#' This function returns a data.frame with the temperatures
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family geothermometer functions
#' @export
Li.Mg.geothermometer <- function(Li, Mg, Temp){
  ndat <- length(Li)
  results <- matrix(0.0, nrow = ndat, ncol = 2)
  #
  term <- log(Li/sqrt(Mg))
  results[,1] <- (2200/(5.47-term))-273.15
  results[,2] <- (1910/(4.63-term))-273.15
  #
  results.df <- as.data.frame(results)
  names(results.df) <- c("Kharaka.Mariner1", "Kharaka.Mariner2")
  return(results.df)
}
#
#Na.Li.geothermometer <- function(Na, Li, Temp){

#}
