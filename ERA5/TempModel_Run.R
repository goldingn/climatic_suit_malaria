project<-function(i){
  # Load libraries
  library(sp)
  library(raster)
  library(Rcpp)
  library(tempsuitcalc)
  
  # Load existing in progress files
  load("ValidCells.RData")
  
  file_names <- c("TempData/2m_temp_2022_03.grib",
                  "TempData/2m_temp_2022_04.grib",
                  "TempData/2m_temp_2022_05.grib",
                  "TempData/2m_temp_2022_06.grib",
                  "TempData/2m_temp_2022_07.grib",
                  "TempData/2m_temp_2022_08.grib",
                  "TempData/2m_temp_2022_09.grib",
                  "TempData/2m_temp_2022_10.grib",
                  "TempData/2m_temp_2022_11.grib",
                  "TempData/2m_temp_2022_12.grib",
                  "TempData/2m_temp_2023_01.grib",
                  "TempData/2m_temp_2023_02.grib",
                  "TempData/2m_temp_2023_03.grib",
                  "TempData/2m_temp_2023_04.grib",
                  "TempData/2m_temp_2023_05.grib")
  
  # Extract the cells of this particular chunk
  valid_cells = na.omit(valid_cells_matrix[i,])
  rm(valid_cells_matrix)
  number_of_cells = length(valid_cells)

  # Extract temperatures from the data files
  extractFromGrib = function(fileName) {
    return(brick(fileName)[valid_cells])
  }
  temp_matrix = do.call(cbind, lapply(file_names, extractFromGrib))[,1:5208]
  
  # Save the output, so if the model needs to be re-run extraction doesn't need to be repeated
  save(valid_cells,
       temp_matrix,
       file = paste("TempMatrix",i,".RData",sep=""))
  
  # Define the model
  riskf <- function(tempraw) {
    # Convert Kelvin data to Celsius
    temp <- tempraw - 273.15
    
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
  nam<-paste("output",i,sep="")
  assign(x=nam,value=output)
  filename<-paste(nam,".RData",sep="")
  save(list=nam,
       file=filename)
}
