project<-function(i){  

  library(Rcpp)
  library(tempsuitcalc)
  
  load("SInputNew.RData")
  
  spac = seq(from = -36,
             to = 398,
             by = (1 / 12))
  sinvec = sin((0:11) / 12 * 2 * pi)
  days = c(-31, 0, 30, 61, 91, 122, 153, 183, 214, 244, 275, 306, 334, 365, 395)
  outputweeks = (36 * 12) + seq(1, 4380, 12)
  
  riskf <- function(tempraw) {
    # Spline Temp Data
    
    means <- spline(days, tempraw[1:15], xout = spac)$y
    sins <- spline(days, tempraw[16:30], xout = spac)$y
    temp <- (means + sins * c(rep(sinvec, 434), 0))
    
    # must used fixed lengths, as specific in Oli's code. He specifies them in hours
    # and then divides by 2:
    # simulation length = 434 days: 36d burn-in, 365d simulation, 33d post-simulation survival
    # L <- 10416 / 2
    # 33 days = maximum vector lifespan
    # L2 <- 792 / 2
    # 2 days = post-emergence non-biting period
    # MnB <- floor(47 / 2)
    
    # degree day accumulation for pathogen
    dd <- pmax(0, temp - 16) / 1332
    
    # degree day accumulation for oviposition. We don't need this, but have to specify it anyway
    ovi <- dd
    
    # per-bin survival probabilities
    surv <-
      exp(-1 / pmax((-52.8 + 15.72 * temp - 0.36 * temp ^ 2), 0))
    
    # replicate for all bin-ages (I.e. assume age has no effect on mortality)
    P <- t(replicate(396, surv))
    
    # vectors for results (these will be modified in place!)
    Zout <- dd * 0
    Zoutovi <- dd * 0
    
    # run difference equations
    . <- tempsuitcalc::DLStempindex(dd, ovi, P, Zout, Zoutovi)
    
    return(Zout[outputweeks])
  }
  
  # Works for array spanning from 1 to 25
  # Crop subset of data
  finalind=min(246772*i,6169293)
  scrop <- ssrawrelval[(246772*i-246771):finalind,]
  rm(ssrawrelval)
  
  output <- apply(scrop, 1, riskf)
  nam<-paste("output",i,sep="")
  assign(x=nam,value=output)
  filename<-paste(nam,".RData",sep="")
  
  save(list=nam,
       file=filename)
  print("Chunk Complete")
}

