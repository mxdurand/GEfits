### Find PPFD at which X percent of Asat is reached
percentLRC <- function(method = "NRH", percent = 0, Pgmax, Rd, phiI0, theta)
{
  Apercent <- Pgmax * percent
  
  
  ## Define Light response functions
  # Non Rectangular hyperbola (most common)
  lcNRH <- function(PARi, Pgmax, Rd, phiI0, theta)
  {
    A <- (((phiI0 * PARi + Pgmax) - sqrt(((phiI0 * PARi + Pgmax)^2) - 4 * theta * phiI0 * PARi * Pgmax)) / (2 * theta)) - Rd
    A <- (A - Apercent)^2
    return(A)
  }
  
  # Rectangular Hyperbola Michaelis-Menten based model
  lcRHMM <- function(PARi, Pgmax, Rd, phiI0)
  {
    A <- (((phiI0 * PARi * Pgmax) / (phiI0 * PARi + Pgmax)) - Rd)
    A <- (A - Apercent)^2
    return(A)
  }
  
  # Hyperbolic tangent based model
  lcHT <- function(PARi, Pgmax, Rd, phiI0)
  {
    A <- Pgmax * tanh((phiI0 * PARi) / Pgmax) - Rd
    A <- (A - Apercent)^2
    return(A)
  }
  
  # Exponential based model
  lcEXP <- function(PARi, Pgmax, Rd, phiI0)
  {
    A <- (Pgmax * (1 - exp((-phiI0 * PARi) / Pgmax))) - Rd
    A <- (A - Apercent)^2
    return(A)
  }
  
  
  if(method == "NRH"){
    lcFUNC <- lcNRH
  } else if(method == "RHMM"){
    lcFUNC <- lcRHMM
  } else if(method == "HT"){
    lcFUNC <- lcHT
  } else if(method == "EXP"){
    lcFUNC <- lcEXP
  } else {
    stop("No other function for now")
  }
  
  
  # fill parameters
  if(!is.null(Pgmax)){formals(lcFUNC)$Pgmax <- Pgmax}
  if(!is.null(Rd)){formals(lcFUNC)$Rd <- Rd}
  if(!is.null(phiI0)){formals(lcFUNC)$phiI0 <- phiI0}
  if(!is.null(theta)){formals(lcFUNC)$theta <- theta}
  
  O <- optimize(f = lcFUNC, interval = c(0,2000))
  Omin <- O$minimum
  return(Omin)
}
