project<-function(i){
  # Load libraries
  library(Rcpp)
  library(tempsuitcalc)
  
  dataMatrixName = paste("TempMatrix",i,".RData",sep="")
  
  if (file.exists(dataMatrixName)) {
    # Load existing in progress files
    load(dataMatrixName)
  } else {
    library(sp)
    library(raster)
    
    # Load existing in progress files
    load("ValidCells.RData")
    
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
    extractTempFromGrib = function(fileName) {
      temp_brick = brick(fileName)
      number_of_layers = nlayers(temp_brick) / 4
      return_matrix = matrix(0, number_of_cells, number_of_layers)
      for (layer in 1:number_of_layers) {
        print(paste(fileName, toString(layer), sep=" "))
        # Convert Kelvin data to Celsius
        return_matrix[,layer] = temp_brick[[4 * layer - 2]][valid_cells] - 273.15
      }
      return(return_matrix)
    }
    temp_matrix = do.call(cbind, lapply(file_names, extractTempFromGrib))[,1:5208]
    
    # Save the output, so if the model needs to be re-run extraction doesn't need to be repeated
    save(valid_cells,
         temp_matrix,
         file = dataMatrixName)
  }
  
  # Define the model
  riskf <- function(temp) {
    outputweeks = (36 * 12) + seq(1, 4380, 12)

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
  
  # Apply the model to temp data
  output <- apply(temp_matrix, 1, riskf)
  
  # Save the output
  nam<-paste("TempOutput",i,sep="")
  assign(x=nam,value=output)
  filename<-paste(nam,".RData",sep="")
  save(list=nam, file=filename)
}
