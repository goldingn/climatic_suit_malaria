project<-function(i){
  # Load libraries
  library(sp)
  library(raster)
  library(Rcpp)
  library(tempsuitcalc)
  
  # Load existing in progress files
  load(paste("TempDewpointMatrix",i,".RData",sep=""))
  
  file_names <- c("AllData/2022_03.grib",
                  "AllData/2022_04.grib",
                  "AllData/2022_05.grib",
                  "AllData/2022_06.grib",
                  "AllData/2022_07.grib",
                  "AllData/2022_08.grib",
                  "AllData/2022_09.grib",
                  "AllData/2022_10.grib",
                  "AllData/2022_11.grib",
                  "AllData/2022_12.grib",
                  "AllData/2023_01.grib",
                  "AllData/2023_02.grib",
                  "AllData/2023_03.grib",
                  "AllData/2023_04.grib",
                  "AllData/2023_05.grib")
  
  # Extract the cells of this particular chunk
  valid_cells = na.omit(valid_cells_matrix[i,])
  rm(valid_cells_matrix)
  number_of_cells = length(valid_cells)

  # Extract temperatures from the data files
  # Running this direct without a for loop seems to cause memory issues
  extractPrecFromGrib = function(fileName) {
    temp_brick = brick(fileName)
    number_of_layers = nlayers(temp_brick) / 4
    return_matrix = matrix(0, number_of_cells, number_of_layers)
    for (layer in 1:number_of_layers) {
      print(paste(fileName, toString(layer), sep=" "))
      return_matrix[,layer] = temp_brick[[4 * layer]][valid_cells]
    }
    return(return_matrix)
  }
  prec_matrix = do.call(cbind, lapply(file_names, extractPrecFromGrib))[,1:5208]
  
  extractEvapFromGrib = function(fileName) {
    temp_brick = brick(fileName)
    number_of_layers = nlayers(temp_brick) / 4
    return_matrix = matrix(0, number_of_cells, number_of_layers)
    for (layer in 1:number_of_layers) {
      print(paste(fileName, toString(layer), sep=" "))
      return_matrix[,layer] = temp_brick[[4 * layer - 1]][valid_cells]
    }
    return(return_matrix)
  }
  evap_matrix = do.call(cbind, lapply(file_names, extractEvapFromGrib))[,1:5208]
  full_matrix = cbind(temp_dewpoint_matrix, prec_matrix, evap_matrix)
  
  # Save the output, so if the model needs to be re-run extraction doesn't need to be repeated
  save(valid_cells,
       full_matrix,
       file = paste("FullMatrix",i,".RData",sep=""))
  
  # Define functions for calculating relative humidity
  mag_coeff_exp = function(temp) {
    return(exp((17.625 * temp) / (243.04 + temp)))
  }
  
  calc_hum = function(temp, dewpoint) {
    return(100 * (mag_coeff_exp(dewpoint) / mag_coeff_exp(temp)))
  }
  
  # Define the model
  riskf <- function(fullraw) {
    outputweeks = (36 * 12) + seq(1, 4380, 12)
    # Convert Kelvin data to Celsius
    temp <- fullraw[1:5208] - 273.15
    dewpoint = fullraw[5209:10416] - 273.15
    # Convert dewpoint data to relative humidity
    hum <- mapply(calc_hum, temp, dewpoint)
    
    # Convert m data to mm
    prec = fullraw[10417:15624] * 1000
    evap = fullraw[15625:20832] * -1000
    
    # Calculate a volume of water curve
    Vout = 0*temp
    .<- tempsuitcalc::WaterVolumeSimWithEvap(temp, prec, evap, Vout)
    
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
  
  # Apply the model to temp data
  output <- apply(temp_dewpoint_matrix, 1, riskf)
  
  # Save the output
  nam<-paste("HumOutput",i,sep="")
  assign(x=nam,value=output)
  filename<-paste(nam,".RData",sep="")
  save(list=nam,
       file=filename)
}
