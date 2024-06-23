project<-function(i){
  # Load libraries
  library(Rcpp)
  library(tempsuitcalc)
  
  dataMatrixName = paste("HumData/TempDewpointMatrix",i,".RData",sep="")
  
  if (file.exists(dataMatrixName)) {
    # Load existing in progress files
    load(dataMatrixName)
  } else {
    library(sp)
    library(raster)
    
    # Load existing in progress files
    load(paste("TempData/TempMatrix",i,".RData",sep=""))
    
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
    
    number_of_cells = length(valid_cells)
  
    # Extract temperatures from the data files
    # Running this direct without a for loop seems to cause memory issues
    extractDewpointFromGrib = function(fileName) {
      temp_brick = brick(fileName)
      number_of_layers = nlayers(temp_brick) / 4
      return_matrix = matrix(0, number_of_cells, number_of_layers)
      for (layer in 1:number_of_layers) {
        print(paste(fileName, toString(layer), sep=" "))
        # Convert Kelvin data to Celsius
        return_matrix[,layer] = temp_brick[[4 * layer - 3]][valid_cells] - 273.15
      }
      return(return_matrix)
    }
    dewpoint_matrix = do.call(cbind, lapply(file_names, extractDewpointFromGrib))[,1:5208]
    temp_dewpoint_matrix = cbind(temp_matrix, dewpoint_matrix)
    
    # Save the output, so if the model needs to be re-run extraction doesn't need to be repeated
    save(valid_cells,
         temp_dewpoint_matrix,
         file = dataMatrixName)
  }
  
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
  filename<-paste("HumOutput/",nam,".RData",sep="")
  save(list=nam, file=filename)
}
