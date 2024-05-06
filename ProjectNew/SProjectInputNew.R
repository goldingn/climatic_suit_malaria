project<-function(i){  
  
  library(Rcpp)
  library(tempsuitcalc)
  
  load("SInputNew.RData")
  
  # Set up initial variables specifying month lengths and time steps
  spac = seq(from = -36,
             to = 398,
             by = (1 / 12))
  sinvec = sin((0:11) / 12 * 2 * pi)
  days = c(-31, 0, 30, 61, 91, 122, 153, 183, 214, 244, 275, 306, 334, 365, 395)
  outputweeks = (36 * 12) + seq(1, 4380, 12)
  
  riskf <- function(tempraw) {
    # Spline temp data to 2-hourly intervals
    means <- spline(days, tempraw[1:15], xout = spac)$y
    sins <- spline(days, tempraw[16:30], xout = spac)$y
    
    # Make teperature curve from mean curve and variance curve
    temp <- (means + sins * c(rep(sinvec, 434), 0))
 
    # Spline vapor pressure and rainfall data
    vapr <- spline(days,tempraw[31:45], xout=spac)$y
    drain <- spline(days,tempraw[46:60], xout=spac)$y
    
    # Calculate humitidy from water vapor pressure
    eqvapr <- (101.325/760)*10^(8.07131-1730.63/(233.426+temp))
    hum <- pmin(100*vapr/eqvapr,100) 
    
    # More Variables
    alt = tempraw[61]
    lat = abs(((3600-floor(tempraw[62]/8640))*150/3600)-60)
    logexp = log(hum/100,base=exp(1))+(7.625*temp/(243.04+temp))
    dewp = 243.04*logexp/(7.625-logexp)
    
    # Calculate a volume of water curve
    Vout = 0*temp
    .<- tempsuitcalc::WaterVolumeSim(temp, dewp, drain, alt, lat, Vout)
    
    # must used fixed lengths, as specific in Oli's code. He specifies them in hours
    # and then divides by 2:
    # simulation length = 434 days: 36d burn-in, 365d simulation, 33d post-simulation survival
    # L <- 10416 / 2
    # 33 days = maximum vector lifespan
    # L2 <- 792 / 2
    # 2 days = post-emergence non-biting period
    # MnB <- floor(47 / 2)
    
    # Convert volume to area
    Aout = (3*sqrt(pi)*Vout)^(2/3)
    
    # Shift area vector to account for delayed hatching
    Aoutmov = c(replicate(36,1000),Aout[1:5173])
   
    # degree day accumulation for pathogen
    dd <- pmax(0, temp - 16) / 1332
    
    # per-bin survival probabilities
    humsurv <- pmax(pmin(((hum-5)/37),1),0)
    
    surv <- exp(-1 / pmax((-52.8 + 15.72 * temp - 0.36 * temp ^ 2), 0))*(humsurv^(1/12))
    
    # replicate for all bin-ages (I.e. assume age has no effect on mortality)
    P <- t(replicate(396, surv))
    
    # vectors for results (these will be modified in place!)
    Zout <- dd * 0
    
    # run difference equations
    . <- tempsuitcalc::DLStempindexnew(dd, P, Zout, Aoutmov)
    
    return(Zout[outputweeks])
  }
  
  # Works for array spanning from 1 to 25
  # Crop subset of data
  finalind=min(246772*i,6169293)
  scrop <- ssrawrelval[(246772*i-246771):finalind,]
  rm(ssrawrelval)
  
  # Run main function
  output <- apply(scrop, 1, riskf)
  
  # Save output
  nam<-paste("output",i,sep="")
  assign(x=nam,value=output)
  filename<-paste(nam,".RData",sep="")
  save(list=nam,
       file=filename)
  
  print("Chunk Complete")
}

