fitLRC <- function(data, varnames = list(Photo = "Photo", PARi = "PARi", Rd = "Rd"), method = "NRH", useRd = F, plot = T, add = F, verbose = T, ...)
{
  ## Define Light response functions
  # Non Rectangular hyperbola (most common)
  lcNRH <- function(formula = T, PARi, Pgmax, Rd, phiI0, theta, ...)
  {
    if(formula){
      return(Photo ~ (((phiI0 * PARi + Pgmax) - sqrt(((phiI0 * PARi + Pgmax)^2) - 4 * theta * phiI0 * PARi * Pgmax)) / (2 * theta)) - Rd)
    } else{
      return((((phiI0 * PARi + Pgmax) - sqrt(((phiI0 * PARi + Pgmax)^2) - 4 * theta * phiI0 * PARi * Pgmax)) / (2 * theta)) - Rd)
    }
  }
  
  # Rectangular Hyperbola Michaelis-Menten based model
  lcRHMM <- function(formula = T, PARi, Pgmax, Rd, phiI0, ...)
  {
    if(formula){
      return(Photo ~ ((phiI0 * PARi * Pgmax) / (phiI0 * PARi + Pgmax)) - Rd)
    } else{
      return((((phiI0 * PARi * Pgmax) / (phiI0 * PARi + Pgmax)) - Rd))
    }
  }
  
  # Hyperbolic tangent based model
  lcHT <- function(formula = T, PARi, Pgmax, Rd, phiI0, ...)
  {
    if(formula){
      return(Photo ~ Pgmax * tanh((phiI0 * PARi) / Pgmax) - Rd)
    } else{
      return(Pgmax * tanh((phiI0 * PARi) / Pgmax) - Rd)
    }
  }
  
  # Exponential based model
  lcEXP <- function(formula = T, PARi, Pgmax, Rd, phiI0, ...)
  {
    if(formula){
      return(Photo ~ (Pgmax * (1 - exp((-phiI0 * PARi) / Pgmax))) - Rd)
    } else{
      return((Pgmax * (1 - exp((-phiI0 * PARi) / Pgmax))) - Rd)
    }
  }
  
  # Get data
  df <- data.frame("Photo" = data[,varnames$Photo], "PARi" = round(data[,varnames$PARi]))
  
  # Select function
  if(method == "NRH"){
    f <- lcNRH()
    lcFUNC <- lcNRH
  } else if(method == "RHMM"){
    f <- lcRHMM()
    lcFUNC <- lcRHMM
  } else if(method == "HT"){
    f <- lcHT()
    lcFUNC <- lcHT
  } else if(method == "EXP"){
    f <- lcEXP()
    lcFUNC <- lcEXP
  } else {
    stop("No other function for now")
  }
  
  # Set Rd value 
  if(useRd){
    Rd <- data[,varnames$Rd]
    if(abs(max(Rd) - min(Rd)) == 0){Rd <- unique(Rd)} else {Rd <- mean(Rd, na.rm = T) ; warning("Rd provided is not a single value, using mean")} # Check if all Rd values are the same, if not use the mean
    Rd <- abs(Rd) # always use positive Rd values
  } else if(nrow(df[df$PARi <= 3,]) > 0){
    Rd <- abs(mean(df[df$PARi <= 3,"Photo"]))
  } else {
    Rd <- 1
  }
  
  # Set upper and lower bounds
  if(method == "NRH"){
    lower_bounds = c(1e-6,1e-6,1e-6)
    upper_bounds = c(20,50,1500)
  } else if(method == "RHMM" | method == "HT" | method == "EXP"){
    lower_bounds = c(1e-6,1e-6)
    upper_bounds = c(20,600)
  } else {
    stop("No other function for now")
  }
  if(useRd){
    lower_bounds <- append(lower_bounds, Rd)
    upper_bounds <- append(upper_bounds, Rd)
    if(verbose == T){
      cat("Using measured Rd = ", Rd, "\n")
    }
  } else {
    lower_bounds <- append(lower_bounds, 1e-6)
    upper_bounds <- append(upper_bounds, 100)
  }
  
  # Find starting value estimation
  if(method == "NRH" | method == "RHMM" | method == "HT" | method == "EXP"){
    
    # First guesses
    Pgmax <- max(df$Photo, na.rm = T) - abs(Rd)
    phiI0 <- unname(coef(lm(Photo~PARi, data = df[df$PARi < 200,]))[2])
    theta <- 1
    
    # function to return sum squares to minimize
    SSfunc <- function(Pgmax, phiI0, Rd, ...)
    {
      vals <- suppressWarnings(lcFUNC(formula = F, PARi = df$PARi, Pgmax = Pgmax, phiI0 = phiI0, Rd = Rd, ...))
      SS <- sum((vals - df$Photo)^2)
      SS[is.nan(SS)] <- 1000
      return(SS)
    }
    
    d = 0.4 ; n = 10
    if(!useRd){
      gg <- expand.grid(Pgmax=seq(Pgmax*(1-d),Pgmax*(1+d),length=n),
                        phiI0=seq(phiI0*(1-d),phiI0*(1+d),length=n),
                        theta=seq(theta*(1-d),theta*(1+d),length=n),
                        Rd=seq(Rd*(1-d),Rd*(1+d),length=n))
      
      m <- with(gg, mapply(SSfunc, Pgmax = Pgmax, phiI0 = phiI0, Rd = Rd, theta = theta))
      Pgmax <- gg$Pgmax[which.min(m)]
      phiI0 <- gg$phiI0[which.min(m)]
      theta <- gg$theta[which.min(m)]
      Rd <- gg$Rd[which.min(m)]
      
      if(verbose == T){
        cat("Starting values are Pgmax = ", round(Pgmax,2), " | phiI0 = ", round(phiI0,2), " | Rd = ", round(Rd,2), " | theta = ", round(theta,2), "(only used in NRH method)", "\n")
      }
    } else {
      gg <- expand.grid(Pgmax=seq(Pgmax*(1-d),Pgmax*(1+d),length=n),
                        phiI0=seq(phiI0*(1-d),phiI0*(1+d),length=n),
                        theta=seq(theta*(1-d),theta*(1+d),length=n))
      
      Rd_real <- Rd
      m <- with(gg, mapply(SSfunc, Pgmax = Pgmax, phiI0 = phiI0, Rd = Rd_real, theta = theta))
      Pgmax <- gg$Pgmax[which.min(m)]
      phiI0 <- gg$phiI0[which.min(m)]
      theta <- gg$theta[which.min(m)]
      
      if(verbose == T){
        cat("Starting values are Pgmax = ", round(Pgmax,2), " | phiI0 = ", round(phiI0,2), " | theta = ", round(theta,2), "(only used in NRH method)", "\n")
      }
    }
  }
  
  if(method == "NRH"){
    fit <- minpack.lm::nlsLM(formula = f, data = df, start = c(phiI0 = phiI0, theta = 1, Pgmax = Pgmax, Rd = Rd), trace = F, lower = lower_bounds, upper = upper_bounds)
    #fit <- onls(formula = f, data = df, start = c(phiI0 = phiI0, theta = 1, Pgmax = Pgmax, Rd = Rd), trace = F, lower = lower_bounds, upper = upper_bounds)
  } else if(method == "RHMM" | method == "HT" | method == "EXP"){
    fit <- minpack.lm::nlsLM(formula = f, data = df, start = c(phiI0 = phiI0, Pgmax = Pgmax, Rd = Rd), trace = F, lower = lower_bounds, upper = upper_bounds)
    #fit <- onls(formula = f, data = df, start = c(phiI0 = phiI0, Pgmax = Pgmax, Rd = Rd), trace = F, lower = lower_bounds, upper = upper_bounds)
  } else {
    stop("No other function for now")
  }
  
  if(plot){
    color = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999") # Color palette
    
    if(add == FALSE){
      par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(5,5,1,1))
      plot(Photo ~ PARi, data = df, pch = 16, cex = 1.4, bty = "L", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
      axis(side = 1, cex.axis = 1, font = 2)
      axis(side = 2, cex.axis = 1, font = 2, las = 2)
      mtext(side = 1, line = 3, cex = 1.2, text = expression(bold(paste("Photosynthetic photon flux density (µmol m"^-2, " s"^-1, ")", sep = ""))))
      mtext(side = 2, line = 2.5, cex = 1.2, text = expression(bold(paste("CO"["2"], " Assimilation (µmol m"^-2, " s"^-1, ")", sep = ""))))
    }
    if(method == "NRH"){
      lines(y = lcFUNC(formula = F, PARi = 0:2000, Pgmax = coef(fit)["Pgmax"], Rd = coef(fit)["Rd"], phiI0 = coef(fit)["phiI0"], theta = coef(fit)["theta"]), x = 0:2000, lwd = 2, col = color[sample(1:length(color),1)])
    } else if(method == "RHMM" | method == "HT" | method == "EXP"){
      lines(y = lcFUNC(formula = F, PARi = 0:2000, Pgmax = coef(fit)["Pgmax"], Rd = coef(fit)["Rd"], phiI0 = coef(fit)["phiI0"]), x = 0:2000, lwd = 2, col = color[sample(1:length(color),1)])
    } 
  }
  
  if(verbose == T){
    cat("RMSE = ", paste(round(sqrt(sum(resid(fit)^2)/length(resid(fit))),3)), "µmol m-2 s-1", "\n")
  }
  return(fit)
}
