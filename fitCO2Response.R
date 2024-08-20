library("plantecophys") # Necessary for Photosyn and findCiTransition functions embedded in the code
library("minpack.lm")


fitCO2gm <- function(data,
                     varnames=list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                     plotit = TRUE,
                     Tcorrect = TRUE,
                     quiet = FALSE,
                     add = FALSE,
                     icol = "red",
                     
                     useRd = FALSE,
                     PPFD = NULL,
                     Tleaf = NULL,
                     gm_guess = NULL,
                     
                     Patm = 100,
                     alpha = 0.24,
                     theta = 0.85,
                     EaV = 82620.87,
                     EdVC = 0,
                     delsC = 645.1013,
                     EaJ = 39676.89,
                     EdVJ = 200000,
                     delsJ = 641.3615,
                     
                     GammaStar = NULL,
                     gstar25 = NULL,
                     Km = NULL)
{
  
  # Get an estimation of Ci transistion and starting values for Vcmax/Jmax/Rd
  fit <- fitCO2(data, varnames, Tcorrect = FALSE, plotit = FALSE, add = FALSE, icol = "red", Patm, quiet, useRd, PPFD, Tleaf, alpha, theta, gmeso = NULL, EaV, EdVC, delsC, EaJ, EdVJ, delsJ, GammaStar, Km)
  transCi <- fit$Ci_transition
  data$Ci <- data[,varnames$Ci] * Patm/100
  
  # Set measured Rd if provided (or warn when provided but not used)
  Rd_meas <- set_Rdmeas(varnames, data, useRd, quiet)
  haveRd <- !is.na(Rd_meas)
  
  if(is.null(gm_guess)){gm_guess = 0.2}
  if(haveRd){
    lower_bounds <- c(1e-15, 1e-15, 1e-15)
    upper_bounds <- c(400, 800, 1)
  } else {
    lower_bounds <- c(1e-15, 1e-15, 1e-15, 1e-15)
    upper_bounds <- c(400, 800, 5, 1)
  }
  
  Tleaf <- if(is.null(Tleaf)){mean(data[,varnames$Tleaf])} else {Tleaf}
  PPFD <- if(is.null(PPFD)){mean(data[,varnames$PPFD])} else {PPFD}
  # NEGATIVE VALUE ARE CHLOROPLASTIC RATES
  gstar25 <- if(is.null(gstar25)){-37.43} else {gstar25} # From Bernacchi 2002
  
  formals(Afunc)$Ci_trans <- transCi
  formals(Afunc)$Tleaf <- Tleaf
  formals(Afunc)$Patm <- Patm
  formals(Afunc)$PPFD <- PPFD
  formals(Afunc)$Alpha <- alpha
  formals(Afunc)$gstar25 <- gstar25
  if(haveRd){formals(Afunc)$Rd <- Rd_meas}
  
  if(!haveRd){
    fitGM <- nlsLM(formula = A ~ Afunc(Ci, Vcmax, Jmax, Rd, gi),
                   data = data,
                   start = list(Vcmax = coef(fit)[1], Jmax = coef(fit)[2], Rd = abs(coef(fit)[3]), gi = gm_guess),
                   lower = lower_bounds, upper = upper_bounds,
                   control=list(minFactor=1/2048,maxiter=1000,tol=0.1), trace = FALSE)
    
    
    Cis <- seq(0,1800,5)
    if(gstar25 < 0){gstar25 <- -gstar25} # Need it positive now
    An_mod <- Afunc(data$Ci, Vcmax = coef(fitGM)[1], Jmax = coef(fitGM)[2], Rd = coef(fitGM)[3], gi = coef(fitGM)[4])
    Ac_mod <- Ac(Ci = Cis, Vcmax = coef(fitGM)[1], Rd = coef(fitGM)[3], gi = coef(fitGM)[4], Tleaf, Patm, gstar25)
    Aj_mod <- Aj(Ci = Cis, Jmax = coef(fitGM)[2], Rd = coef(fitGM)[3], gi = coef(fitGM)[4], Tleaf, Patm, PPFD, Alpha = alpha, gstar25)
  } else {
    fitGM <- nlsLM(formula = A ~ Afunc(Ci, Vcmax, Jmax, Rd_meas, gi),
                   data = data,
                   start = list(Vcmax = coef(fit)[1], Jmax = coef(fit)[2], gi = gm_guess),
                   lower = lower_bounds, upper = upper_bounds,
                   control=list(minFactor=1/2048,maxiter=1000,tol=0.1), trace = FALSE)
    
    Cis <- seq(0,1800,5)
    An_mod <- Afunc(data$Ci, Vcmax = coef(fitGM)[1], Jmax = coef(fitGM)[2], Rd = Rd_meas, gi = coef(fitGM)[3])
    Ac_mod <- Ac(Ci = Cis, Vcmax = coef(fitGM)[1], Rd = Rd_meas, gi = coef(fitGM)[3], Tleaf, Patm, gstar25)
    Aj_mod <- Aj(Ci = Cis, Jmax = coef(fitGM)[2], Rd = Rd_meas, gi = coef(fitGM)[3], Tleaf, Patm, PPFD, Alpha = alpha, gstar25)
  }
  rmse_fit <- sqrt(sum(fitGM$m$resid()^2)/length(fitGM$m$resid()))
  
  if(Tcorrect == TRUE){
    Vcmax_i <- coef(fitGM)[1] / TVcmax(Tleaf,EaV, delsC, EdVC)
    Jmax_i <- coef(fitGM)[2] / TJmax(Tleaf,EaJ, delsJ, EdVJ)
  } else {
    Vcmax_i <- coef(fitGM)[1]
    Jmax_i <- coef(fitGM)[2]
  }
  
  if(plotit == TRUE){
    if(add == FALSE){
      plot(-500, xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(0,1800), ylim = c(-5,max(data$A)*1.1), bty = "L")
      points(A ~ Ci, data, pch = 21, col = "black", bg = "gray85", cex = 1.4)
      axis(side = 1, font = 2)
      axis(side = 2, font = 2, las = 2)
      mtext(side = 1, line = 2.5, text = "Ci (ppm)", font = 2, cex = 1.2)
      mtext(side = 2, line = 2.5, text = expression(bold(paste("A"["net"], " (µmol m"^"-2", " s"^"-1", ")", sep = ""))), cex = 1.2)
    }
    points(Ac_mod ~ Cis, type = "l", lty = 1, lwd = 2, col = icol)
    points(Aj_mod ~ Cis, type = "l", lty = 2, lwd = 2, col = icol)
    if(!haveRd){
      if(Tcorrect == TRUE){
        print(paste0("At 25 °C: ", "Vcmax = ", round(Vcmax_i,2), " | Jmax = ", round(Jmax_i,2), " | Rd = ", round(coef(fitGM)[3],2), " | gm = ", round(coef(fitGM)[4],3), " | RMSE = ", round(rmse_fit,3)))
      } else {
        print(paste0("At ", round(Tleaf,1)," °C: ", "Vcmax = ", round(coef(fitGM)[1],2), " | Jmax = ", round(coef(fitGM)[2],2), " | Rd = ", round(coef(fitGM)[3],2), " | gm = ", round(coef(fitGM)[4],3), " | RMSE = ", round(rmse_fit,3)))
      }
    } else {
      if(Tcorrect == TRUE){
        print(paste0("At 25 °C: ", "Vcmax = ", round(Vcmax_i,2), " | Jmax = ", round(Jmax_i,2), " | Rd = ", round(Rd_meas,2), " | gm = ", round(coef(fitGM)[4],3), " | RMSE = ", round(rmse_fit,3)))
      } else {
        print(paste0("At ", round(Tleaf,1)," °C: ", "Vcmax = ", round(coef(fitGM)[1],2), " | Jmax = ", round(coef(fitGM)[2],2), " | Rd = ", round(Rd_meas,2), " | gm = ", round(coef(fitGM)[4],3), " | RMSE = ", round(rmse_fit,3)))
      }
    }
  }
  l <- list()
  if(!haveRd){
    l$coef <- c(Vcmax_i, Jmax_i, coef(fitGM)[3], coef(fitGM)[4])
  } else {
    l$coef <- c(Vcmax_i, Jmax_i, Rd_meas, coef(fitGM)[3])
  }
  
  l$rmse <- rmse_fit
  l$fit <- fitGM
  return(l)
}


### A-Ci fitting with infinite gm -----------------------------------------------------
fitCO2 <- function(data,
                   varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                   Tcorrect = TRUE,
                   plotit = TRUE,
                   add = FALSE,
                   icol = "red",
                   
                   Patm = 100,
                   quiet = FALSE,
                   
                   useRd = FALSE,
                   PPFD = NULL,
                   Tleaf = NULL,
                   
                   alpha = 0.24,
                   theta = 0.85,
                   gmeso = NULL,
                   EaV = 82620.87,
                   EdVC = 0,
                   delsC = 645.1013,
                   EaJ = 39676.89,
                   EdVJ = 200000,
                   delsJ = 641.3615,
                   
                   GammaStar = NULL,
                   Km = NULL)
{
  ### first are checks and making sure data is given properly
  
  # Defaults
  citransition=NULL
  startValgrid=TRUE
  fitmethod="bilinear"
  algorithm="default"
  fitTPU=FALSE
  alphag=0
  
  
  # Check if Km and gammastar are given
  gstarinput <- !is.null(GammaStar)
  kminput <- !is.null(Km)
  if(is.null(gmeso))gmeso <- -999  # cannot pass NULL value to nls
  
  # Check data 
  if(nrow(data) == 0){
    Stop("No rows in data - check observations.")
  }
  
  # Make sure data is a proper dataframe, not some tibble or other nonsense.
  data <- as.data.frame(data)
  
  # Add PPFD and Tleaf to data, if needed (uses default values, or input values)
  data <- set_PPFD(varnames, data, PPFD, quiet)
  data <- set_Tleaf(varnames, data, Tleaf, quiet)
  
  # Set measured Rd if provided (or warn when provided but not used)
  Rd_meas <- set_Rdmeas(varnames, data, useRd, quiet)
  haveRd <- !is.na(Rd_meas)
  
  # Extract Ci and apply pressure correction
  data$Ci_original <- data[,varnames$Ci]
  data$Ci <- data[,varnames$Ci] * Patm/100
  
  # Extract measured net leaf photosynthesis
  data$ALEAF <- data[,varnames$ALEAF]
  
  # Calculate Km and GammaStar, if not input
  Km_v <- if(!kminput){
    TKm(data$Tleaf, Patm)
  } else {
    rep(Km, nrow(data))
  } 
  
  GammaStar_v <- if(!gstarinput){
    TGammaStar(data$Tleaf, Patm)
  } else {
    rep(GammaStar, nrow(data))
  }
  
  ### Fitting
  f <- do_fit_method_bilinear_bestcitrans(data, haveRd, fitTPU, alphag, Rd_meas, 
                                          Patm, Tcorrect, algorithm,
                                          alpha,theta,gmeso,EaV,EdVC,delsC,
                                          EaJ,EdVJ,delsJ,
                                          GammaStar_v, Km_v)
  
  # TPU. If it is not estimated, put something in because we need it later.
  if(!("TPU" %in% names(f)))f$TPU <- 1000
  
  # Only used to add 'Amodel' to the output
  acirun <- do_acirun(data,f,Patm,Tcorrect,  # Note several parameters are stored in 'f'
                      alpha=alpha,theta=theta,
                      gmeso=gmeso,EaV=EaV,
                      EdVC=EdVC,delsC=delsC,
                      EaJ=EaJ,EdVJ=EdVJ,
                      delsJ=delsJ,GammaStar=GammaStar_v, Km=Km_v)
  
  # If Ap is never actually limiting, set estimated TPU to 1000 (TPU not limiting)
  if(!any(acirun$Ap < acirun$Aj))f$TPU <- 1000
  
  ### Organize output
  l <- list()  
  runorder <- order(acirun$Ci)
  l$df <- acirun[runorder,]
  l$pars <- f$pars
  l$nlsfit <- f$fit
  l$Tcorrect <- Tcorrect
  
  formals(Photosyn)$Tleaf <- mean(data$Tleaf)
  formals(Photosyn)$Patm <- Patm
  formals(Photosyn)$PPFD <- mean(data$PPFD)
  formals(Photosyn)$Vcmax <- l$pars[1]
  formals(Photosyn)$Jmax <- l$pars[2]
  formals(Photosyn)$Rd <- l$pars[3]
  formals(Photosyn)$TPU <- 1000
  formals(Photosyn)$alphag <- alphag
  formals(Photosyn)$Tcorrect <- Tcorrect
  formals(Photosyn)$alpha <- alpha
  formals(Photosyn)$theta <- theta
  formals(Photosyn)$gmeso <- gmeso
  formals(Photosyn)$EaV <- EaV
  formals(Photosyn)$EdVC <- EdVC
  formals(Photosyn)$delsC <- delsC
  formals(Photosyn)$EaJ <- EaJ
  formals(Photosyn)$EdVJ <- EdVJ
  formals(Photosyn)$delsJ <- delsJ
  if(gstarinput)formals(Photosyn)$GammaStar <- GammaStar
  if(kminput)formals(Photosyn)$Km <- Km
  
  l$Photosyn <- Photosyn
  
  # The inverse - find Ci given a rate of photosyntheiss
  l$Ci <- function(ALEAF){
    
    O <- function(ci, photo){
      (l$Photosyn(Ci=ci)$ALEAF - photo)^2
    }
    o <- optimize(O, interval=c(1,10^5), photo =ALEAF)
    if(o$objective > 1e-02){
      return(NA)
    } else {
      return(o$minimum)
    }
  }
  
  # Transition points.
  trans <- findCiTransition(l$Photosyn)
  l$Ci_transition <- trans[1]
  l$Ci_transition2 <- trans[2]
  l$Rd_measured <- haveRd
  
  # Save GammaStar and Km (either evaluated at mean temperature, 
  # or input if provided)
  l$GammaStar <- mean(GammaStar_v)
  l$Km <- mean(Km_v)
  l$kminput <- kminput
  l$gstarinput <- gstarinput
  
  l$fitmethod <- fitmethod
  l$citransition <- ifelse(is.null(citransition), NA, citransition)
  l$gmeso <- ifelse(is.null(gmeso) || gmeso < 0, NA, gmeso)
  l$fitTPU <- fitTPU
  l$alphag <- alphag
  l$RMSE <- rmse_acifit(l)
  l$runorder <- runorder
  
  class(l) <- "acifit"
  
  if(plotit == TRUE){
    if(Tcorrect == T){
      # Recalculate parameter at Tleaf for plotting
      Vcmax_i <- l$pars[1] * TVcmax(mean(l$df$Tleaf),EaV, delsC, EdVC)
      Jmax_i <- l$pars[2] * TJmax(mean(l$df$Tleaf),EaJ, delsJ, EdVJ)
    } else {
      Vcmax_i <- l$pars[1]
      Jmax_i <- l$pars[2]
    }
    
    Cis <- seq(0,1800,50)
    Ac <- Vcmax_i*(Cis-TGammaStar(mean(l$df$Tleaf), Patm = Patm))/(Cis+TKm(mean(l$df$Tleaf), Patm = Patm))-l$pars[3]
    # J <- alpha*mean(l$df$PPFD)/sqrt(1+alpha^2*mean(l$df$PPFD)^2/Jmax_i^2)
    J <- Jfun(PPFD = mean(l$df$PPFD), alpha = alpha, Jmax = Jmax_i, theta = theta)
    Aj <- J*(Cis-TGammaStar(mean(l$df$Tleaf), Patm = Patm))/(4*Cis+8*TGammaStar(mean(l$df$Tleaf), Patm = Patm))-l$pars[3]
    
    if(add == FALSE){
      plot(-500, xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(0,1800), ylim = c(-5,max(l$df$Ameas)*1.1), bty = "L")
      points(l$df$Ameas ~ l$df$Ci, pch = 21, col = "black", bg = "gray85", cex = 1.4)
      axis(side = 1, font = 2)
      axis(side = 2, font = 2, las = 2)
      mtext(side = 1, line = 2.5, text = "Ci (ppm)", font = 2, cex = 1.2)
      mtext(side = 2, line = 2.5, text = expression(bold(paste("A"["net"], " (µmol m"^"-2", " s"^"-1", ")", sep = ""))), cex = 1.2)
    }
    points(Ac ~ Cis, type = "l", lty = 1, lwd = 2, col = icol)
    points(Aj ~ Cis, type = "l", lty = 2, lwd = 2, col = icol)
    if(Tcorrect == TRUE){
      print(paste0("At 25 °C: ", "Vcmax = ", round(l$pars[1],2), " | Jmax = ", round(l$pars[2],2), " | Rd = ", round(l$pars[3],2), " | RMSE = ", round(l$RMSE,3)))
    } else {
      print(paste0("At ", round(mean(l$df$Tleaf),1)," °C: ", "Vcmax = ", round(l$pars[1],2), " | Jmax = ", round(l$pars[2],2), " | Rd = ", round(l$pars[3],2), " | RMSE = ", round(l$RMSE,3)))
    }
  }
  return(l)
}
# End of function

### Fitting functions -------------------------------------------------------
do_acirun <- function(data,f,Patm,Tcorrect,...){
  
  acirun <- Photosyn(Ci=data$Ci, 
                     Vcmax=f$pars[1,1], Jmax=f$pars[2,1], 
                     Rd=f$pars[3,1], 
                     TPU=f$TPU,
                     PPFD=data$PPFD, 
                     Tleaf=data$Tleaf,
                     Patm=Patm,
                     Tcorrect=Tcorrect,...)
  
  acirun$Ameas <- data$ALEAF
  acirun$ELEAF <- NULL
  acirun$GS <- NULL
  acirun$Ca <- NULL
  acirun$Ci_original <- data$Ci_original
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  
  # shuffle
  avars <- match(c("Ci","Ameas","Amodel"),names(acirun))
  acirun <- acirun[,c(avars, setdiff(1:ncol(acirun), avars))]
  
  return(acirun)
}

do_fit_method_bilinear <- function(data, haveRd, alphag, Rd_meas, 
                                   Patm, citransition, citransition2=NULL,
                                   Tcorrect, algorithm,
                                   alpha,theta,gmeso,EaV,EdVC,delsC,
                                   EaJ,EdVJ,delsJ,
                                   GammaStar, Km, onepoint=FALSE){
  
  # Calculate T-dependent parameters
  ppar <- Photosyn(Tleaf=data$Tleaf, Patm=Patm, Tcorrect=Tcorrect,
                   alpha=alpha,theta=theta,
                   gmeso=gmeso,EaV=EaV,
                   EdVC=EdVC,delsC=delsC,
                   EaJ=EaJ,EdVJ=EdVJ,
                   delsJ=delsJ,
                   returnParsOnly=TRUE)
  if(!is.null(GammaStar))ppar$GammaStar <- GammaStar
  if(!is.null(Km))ppar$Km <- Km
  if(is.null(citransition2))citransition2 <- 10^6
  
  # Calculate Cc if gmeso included
  if(!is.null(gmeso) && gmeso > 0){
    Conc <- data$Ci - data$ALEAF/gmeso
  } else {
    Conc <- data$Ci
  }
  
  # Linearize
  data$vcmax_pred <- (Conc - ppar$GammaStar)/(Conc + ppar$Km)
  data$Jmax_pred <- (Conc - ppar$GammaStar)/(Conc + 2*ppar$GammaStar)
  data$TPU_part <- (Conc - ppar$GammaStar)/(Conc - (1+3*alphag)*ppar$GammaStar)
  
  # Fit Vcmax and Rd from linearized portion
  datv <- data[data$Ci < citransition & data$Ci < citransition2,]
  if(nrow(datv) == 0){
    return(list(pars=NA, fit=NA, TPU=NA, success=FALSE))
  }
  
  if(!haveRd){
    fitv <- lm(ALEAF ~ vcmax_pred, data=datv)
    Rd_fit <- coef(fitv)[[1]]
    Vcmax_fit <- coef(fitv)[[2]]
  } else {
    # If using measured Rd, add to Anet, and remove intercept from fit.
    datv$ALEAFg <- datv$ALEAF + Rd_meas
    fitv <- lm(ALEAFg ~ vcmax_pred-1, data=datv)
    Rd_fit <- -Rd_meas
    Vcmax_fit <- coef(fitv)[[1]]
  }
  
  # Fit Jmax from linearized portion
  datj <- data[data$Ci >= citransition & data$Ci < citransition2,]
  datp <- data[data$Ci >= citransition2,]
  
  # Manual fix: if only one point for TPU, and none for Jmax, abandon fit.
  # In this case it would be more defensible to use the single point for Jmax.
  if(nrow(datp) == 1 && nrow(datj) == 0){
    return(list(pars=NA, fit=NA, TPU=NA, success=FALSE))
  }
  
  
  if(nrow(datj) > 0){
    
    # Fit gross photo using fitted Rd
    datj$Agross <- datj$ALEAF - Rd_fit # Rd_fit is negative
    
    # One point, calculate directly
    if(nrow(datj) == 1){
      J_fit <- with(datj, 4 * Agross / Jmax_pred)
    } else {
      fitj <- lm(Agross ~ Jmax_pred-1, data=datj)
      J_fit <- 4 * coef(fitj)[[1]]
    }
    
    # And solve for Jmax from inverse non-rect. hyperbola
    Jmax_fit <- inverseJfun(mean(data$PPFD), alpha, J_fit, theta)
    
  } else {
    Jmax_fit <- 10^6  # not elegant but will do for now 
    # (avoids trouble elsewhere)
  }
  
  # TPU
  if(nrow(datp) > 0){
    datp$Agross <- datp$ALEAF - Rd_fit
    
    tpu_vals <- (1/3) * datp$Agross / datp$TPU_part
    
    TPU <- mean(tpu_vals)
  } else {
    TPU <- 1000  # same as default in Photosyn
  }
  
  # The above estimates are at the measured Tleaf.
  # Express at 25C?
  if(Tcorrect){
    Jmax_fit <- Jmax_fit / TJmax(mean(data$Tleaf), EaJ, delsJ, EdVJ)
    Vcmax_fit <- Vcmax_fit / TVcmax(mean(data$Tleaf),EaV, delsC, EdVC)
  }
  
  ses <- summary(fitv)$coefficients[,2]
  if(!haveRd){
    pars <- matrix(c(Vcmax_fit, Jmax_fit, -Rd_fit,
                     ses[2],NA,ses[1]), ncol=2)
  } else {
    pars <- matrix(c(Vcmax_fit, Jmax_fit, -Rd_fit,
                     ses[1],NA,NA), ncol=2)
  }
  
  
  rownames(pars) <- c("Vcmax","Jmax","Rd")
  colnames(pars) <- c("Estimate","Std. Error")
  
  
  return(list(pars=pars, fit=fitv, TPU=TPU, success=TRUE))
}


do_fit_method_bilinear_bestcitrans <- function(data, haveRd, fitTPU, alphag, 
                                               Rd_meas, Patm, Tcorrect,
                                               algorithm,alpha,theta,gmeso,EaV,
                                               EdVC,delsC,EaJ,EdVJ,
                                               delsJ,GammaStar, Km){
  
  # Possible Ci transitions
  ci <- data$Ci
  nci <- length(ci)
  citransitions <- diff(ci)/2 + ci[-nci] # Get mid points between Cis
  
  # at least two Ci values to estimate Vcmax and Rd, so delete first
  citransitions1 <- citransitions[-1]
  citransitions2 <- max(ci) + 1  # outside range, on purpose
  
  citransdf <- expand.grid(ci1=citransitions1, ci2=citransitions2)
  citransdf <- citransdf[citransdf$ci1 <= citransdf$ci2,]
  SS <- c()
  
  # Note that Tcorrect is set to FALSE inside the loop. If Tcorrect is needed, it is done
  # after the loop finishes (avoids a bug).
  for(i in seq_len(nrow(citransdf))){
    
    fit <- do_fit_method_bilinear(data, haveRd, alphag, Rd_meas, Patm, 
                                  citransdf$ci1[i], citransdf$ci2[i], 
                                  Tcorrect=FALSE, algorithm,
                                  alpha,theta,gmeso,EaV,EdVC,delsC,
                                  EaJ,EdVJ,delsJ,
                                  GammaStar, Km)
    
    if(fit$success && !any(is.na(fit$pars[,"Estimate"]))){
      run <- do_acirun(data,fit,Patm,Tcorrect=FALSE,
                       alpha=alpha,theta=theta,
                       gmeso=gmeso,EaV=EaV,
                       EdVC=EdVC,delsC=delsC,
                       EaJ=EaJ,EdVJ=EdVJ,
                       delsJ=delsJ,GammaStar=GammaStar,Km=Km)
      
      SS[i] <- sum((run$Ameas - run$Amodel)^2)  
    } else {
      SS[i] <- 10^6
    }
  }
  
  # Best Ci transitions
  bestcis <- citransdf[which.min(SS),]
  
  f <- do_fit_method_bilinear(data, haveRd, alphag, Rd_meas, Patm, 
                              bestcis$ci1, bestcis$ci2, Tcorrect, algorithm,
                              alpha,theta,gmeso,EaV,EdVC,delsC,EaJ,EdVJ,delsJ,
                              GammaStar, Km)
  
  if(f$pars["Jmax","Estimate"] < 0){
    Stop("Cannot invert light response curve to estimate Jmax - increase alpha or theta.")
  }
  
  
  return(f)  
}


rmse_acifit <- function(x)sqrt(sum((x$df$Ameas - x$df$Amodel)^2))


### Other sub-functions -----------------------------------------------------

.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}

TGammaStar <- function(Tleaf, Patm,
                       Egamma=37830.0, 
                       value25=42.75){  
  
  value25*arrh(Tleaf,Egamma)*Patm/100
}

# Find GammaStar at 25 degrees if measured at different temperature
GammaStar25 <- function(gstar, Tleaf, Patm, Egamma=37830.0){
  gstar / ((Patm / 100) * arrh(Tleaf,Egamma))
}

TKm <- function(Tleaf, Patm,
                Oi = 210,      # O2 concentration (mmol mol-1)
                Ec = 79430.0,  # activation energy for Kc 
                Eo = 36380.0,  # activation energy for Ko
                Kc25 = 404.9,  # Kc at 25C
                Ko25 = 278.4  # Ko at 25C
){
  
  Oi <- Oi * Patm / 100
  
  Ko <- Ko25*arrh(Tleaf, Eo)
  Kc <- Kc25*arrh(Tleaf, Ec)
  Km <- Kc * (1.0 + Oi / Ko)
  
  return(Km)
}

# Vcmax temperature response (Arrhenius)
TVcmax <- function(Tleaf, EaV, delsC, EdVC){
  
  if(EdVC > 0){
    V1 <- 1+exp((delsC*(25 + 273.15)-EdVC)/(.Rgas()*(25 + 273.15)))
    V2 <- 1+exp((delsC*(Tleaf+273.15)-EdVC)/(.Rgas()*(Tk(Tleaf))))
    f <- V1/V2
  } else f <- 1
  
  exp((Tleaf-25)*EaV/(.Rgas()*Tk(Tleaf)*Tk(25))) * f
}

# Jmax temperature response (Arrhenius)
TJmax <- function(Tleaf, EaJ, delsJ, EdVJ){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/.Rgas()/298.15)
  J2 <- 1+exp((Tk(Tleaf)*delsJ-EdVJ)/.Rgas()/Tk(Tleaf))
  exp(EaJ/.Rgas()*(1/298.15 - 1/Tk(Tleaf)))*J1/J2
}

# Non-rectangular hyperbola
Jfun <- function(PPFD, alpha, Jmax, theta){
  (alpha*PPFD + Jmax - 
     sqrt((alpha*PPFD + Jmax)^2 - 4*alpha*theta*PPFD*Jmax))/(2*theta)
}

# Given PPFD and J, what is Jmax? (inverse non-rectangular hyperbola)
inverseJfun <- function(PPFD, alpha, J, theta){
  J*(J*theta - alpha*PPFD)/(J - alpha*PPFD)
}

# Set PPFD
set_PPFD <- function(varnames, data, PPFD, quiet){
  
  if(!varnames$PPFD %in% names(data) & is.null(PPFD)){
    data$PPFD <- 1800
    if(!quiet)Warning("PARi not in dataset; assumed PARi = 1800.")
  } else {
    if(!is.null(PPFD))
      data$PPFD <- PPFD
    else
      data$PPFD <- data[,varnames$PPFD]
  }
  return(data)
}


# Set Tleaf
set_Tleaf <- function(varnames, data, Tleaf, quiet){
  # Check if Tleaf is provided
  if(!varnames$Tleaf %in% names(data) & is.null(Tleaf)){
    data$Tleaf <- 25
    if(!quiet)Warning("Tleaf not in dataset; assumed Tleaf = 25.")
  } else {
    if(!is.null(Tleaf))
      data$Tleaf <- Tleaf
    else
      data$Tleaf <- data[,varnames$Tleaf]
  }
  return(data)
}


# Set Rd
set_Rdmeas <- function(varnames, data, useRd, quiet){
  
  Rd_meas <- NA
  if(!is.null(varnames$Rd)){ 
    if(varnames$Rd %in% names(data) && useRd){
      
      # Has to be a single unique value for this dataset
      Rd_meas <- data[,varnames$Rd]
      Rd_meas <- unique(Rd_meas)
      if(length(Rd_meas) > 1){
        Stop("If Rd provided as measured, it must be a single",
             "\nunique value for an A-Ci curve.")
      }
      
      # Use positive value throughout.
      Rd_meas <- abs(Rd_meas)
      haveRd <- TRUE
    }
    if(varnames$Rd %in% names(data) && !useRd){
      if(!quiet)
        message(paste("Rd found in dataset but useRd set to FALSE.",
                      "Set to TRUE to use measured Rd."))
    }
  }
  return(Rd_meas)
}

### Sub-functions for gm fitting -----------------------------------------------------

Ac <- function(Ci, Vcmax, Rd, gi, Tleaf, Patm, gstar25_cc)
{
  if(gstar25_cc < 0){stop("gammastar negative")}
  a <- -1 / gi
  b <- (Vcmax - Rd) / gi + Ci + TKm(Tleaf, Patm, Ec = 80990, Eo = 23720, Kc25 = 272.38, Ko25 = 165.825)
  c <- Rd * (Ci + TKm(Tleaf, Patm, Ec = 80990, Eo = 23720, Kc25 = 272.38, Ko25 = 165.825)) - Vcmax * (Ci - TGammaStar(Tleaf, Patm, value25 = gstar25_cc, Egamma = 24460))
  return((-b + sqrt(b^2 - 4 * a * c)) / (2 * a))
}

Aj <- function(Ci, Jmax, Rd, gi, Tleaf, Patm, PPFD, Alpha, gstar25_cc)
{
  if(gstar25_cc < 0){stop("gammastar negative")}
  J <- Alpha * PPFD / sqrt(1 + Alpha^2 * PPFD^2 / Jmax^2)
  a <- -1 / gi
  b <- (J / 4 - Rd) / gi + Ci + 2 * TGammaStar(Tleaf, Patm, value25 = gstar25_cc, Egamma = 24460)
  c <- Rd * (Ci + 2 * TGammaStar(Tleaf, Patm, value25 = gstar25_cc, Egamma = 24460)) - J * (Ci - TGammaStar(Tleaf, Patm, value25 = gstar25_cc, Egamma = 24460)) / 4
  return((-b + sqrt(b^2 - 4 * a * c)) / (2 * a))
}

Afunc <- function(Ci, Vcmax, Jmax, Rd, gi, Ci_trans = transCi, Tleaf = mean(data$Tleaf), Patm = mean(data$Pa), PPFD = mean(data$Qin), Alpha = 0.24, gstar25 = 37.430001)
{
  Ci_c <- Ci[Ci < Ci_trans]
  Ci_j <- Ci[Ci > Ci_trans]
  
  # Negative values are for Cc rates, positives for Ci
  if(gstar25 > 0){
    gstar25_cc = gstar25 - (Rd / gi)
  } else {
    gstar25_cc = -gstar25
  }
  if(gstar25_cc < 0){stop("gammastar negative")}
  
  Acarb <- Ac(Ci_c, Vcmax, Rd, gi, Tleaf, Patm, gstar25_cc)
  Aelec <- Aj(Ci_j, Jmax,  Rd, gi, Tleaf, Patm, PPFD, Alpha, gstar25_cc)
  Anet <- c(Acarb, Aelec)
}

QUADP <- function(A,B,C){
  
  if((B^2 - 4*A*C) < 0){
    warning("IMAGINARY ROOTS IN QUADRATIC")
    return(0)
  }
  
  if(identical(A,0)){
    if(identical(B,0)){
      return(0)
    } else {
      return(-C/B)
    }
  } else {
    return((- B + sqrt(B^2 - 4*A*C)) / (2*A))
  }
  
}
