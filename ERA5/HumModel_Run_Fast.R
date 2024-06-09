project<-function(i){
  # Load libraries
  library(Rcpp)
  library(tempsuitcalc)
  
  # Load existing in progress files
  load(paste("TempDewpointMatrix",i,".RData",sep=""))
  
  # Define functions for calculating relative humidity
  mag_coeff_exp = function(temp) {
    return(exp((17.625 * temp) / (243.04 + temp)))
  }
  
  calc_hum = function(temp, dewpoint) {
    return(100 * (mag_coeff_exp(dewpoint) / mag_coeff_exp(temp)))
  }
  
  # Define the model
  riskf <- function(tempdew) {
    outputweeks = (36 * 12) + seq(1, 4380, 12)
    temp <- tempdew[1:5208]
    dewpoint = tempdew[5209:10416]
    # Convert dewpoint data to relative humidity
    hum <- mapply(calc_hum, temp, dewpoint)
    
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
    
    humsurv <- pmax(pmin(((hum-5)/37),1),0)
    
    surv <-
      exp(-1 / pmax((-52.8 + 15.72 * temp - 0.36 * temp ^ 2), 0))*(humsurv^(1/12))
    
    # replicate for all bin-ages (I.e. assume age has no effect on mortality)
    P <- t(replicate(396, surv))
    
    # vectors for results (these will be modified in place!)
    Zout <- dd * 0
    Zoutovi <- dd * 0
    
    # run difference equations
    . <- tempsuitcalc::DLStempindex(dd, ovi, P, Zout, Zoutovi)
    
    return(Zout[outputweeks])
  }
  
  # Apply the model to temp data
  output <- apply(temp_dewpoint_matrix, 1, riskf)
  
  # Save the output
  nam<-paste("HumOutput",i,sep="")
  assign(x=nam,value=output)
  filename<-paste(nam,".RData",sep="")
  save(list=nam, file=filename)
}
