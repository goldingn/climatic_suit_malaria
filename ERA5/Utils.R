# Take a rasterlayer cell number and convert it to a lat/long.
# This function assumes the rasterlayer is of the entire available region
convertCellNumberIntoCoords = function(cell) {
  offsetCell = cell - 1
  xcoord = 90 - (floor(offsetCell / 3600) / 10)
  ycoord = ((offsetCell %% 3600) / 10) - 180
  return (c(xcoord, ycoord))
}

# Take a lat/long and convert it to a rasterlayer cell
# This function assumes the rasterlayer is of the entire available region
convertCoordsIntoCellNumber = function(coords) {
  numberCellsDown = (90 - coords[1]) * 10
  numberCellsAcross = (coords[2] + 180) * 10
  return ((numberCellsDown * 3600) + numberCellsAcross)
}

# Take the outputs of a model (TempOutput, HumOutput and RainfallOutput), and plot some pngs
# maxz controls the upper limit of the plot. For temp and hum models, this can be set to 2500. Rainfall model should set this to 4000.
plotRasterize = function(model, maxz) {
  library(raster)
  library(RColorBrewer)
  # First day of the model is 6th April 2022
  # These days and titles will provide snapshots on the 15th of each month. This can be changed to produce a different set of outputs.
  days<-c(285, 316, 344, 10, 40, 71, 101, 132, 163, 193, 224, 254)
  titles=c("January", "February", "March", "April", "May", "June",
           "July", "August", "September", "October", "November", "December")
  mont<-matrix(nrow=length(days),ncol=0)
  
  load('ValidCells.RData')
  valid_cells = as.vector(na.omit(as.vector(t(valid_cells_matrix))))
  rm(valid_cells_matrix)
  
  colours=brewer.pal(9,"YlOrRd")[2:9]
  palette=colorRampPalette(colours)
  
  for(i in 1:25) {
    load(paste(model, '/', model, i,'.RData',sep=''))
    temp<-get(paste(model,i,sep=''))
    mont<-cbind(mont,temp[days,])
    rm(list=paste(model,i,sep=''))
  }
  for(index in 1:length(days)) {
    blankRaster = raster('blankNew.gri')
    blankRaster[valid_cells] = mont[index,]
    png(filename = paste(model,'/',model,'Image',index,'.png',sep=''), height=720, width=949)
    plot(blankRaster,col=palette(100),axes=FALSE,zlim=c(0,maxz))
    title(titles[index], cex.main=2.5)
    dev.off()
  }
}

# Take the outputs of a model (TempOutput, HumOutput and RainfallOutput), aggregates the results over the 365 days, and produces an output image
# The output raster is scaled relative to the model's maxValue
plotAggregate = function(model) {
  library(raster)
  library(RColorBrewer)
  
  agg = vector()
  
  load('ValidCells.RData')
  valid_cells = as.vector(na.omit(as.vector(t(valid_cells_matrix))))
  rm(valid_cells_matrix)
  
  colours=brewer.pal(9,"YlOrRd")[2:9]
  palette=colorRampPalette(colours)
  
  for(i in 1:25) {
    load(paste(model, '/', model, i,'.RData',sep=''))
    temp<-get(paste(model,i,sep=''))
    agg = c(agg, colMeans(temp))
    rm(list=paste(model,i,sep=''))
  }
  outputRaster = raster('blankNew.gri')
  outputRaster[valid_cells] = agg
  png(filename = paste(model,'Agg.png',sep=''), height=720, width=949)
  zmax = maxValue(outputRaster)
  outputRaster = outputRaster / zmax
  plot(outputRaster,col=palette(100),axes=FALSE,zlim=c(0,1))
  dev.off()
  writeRaster(outputRaster, filename=paste(model,'Agg.grd',sep=''),overwrite=TRUE)
}

# Take the outputs of a model (TempOutput, HumOutput and RainfallOutput), aggregates the results over the seasons, and produces 4 output images
# maxz controls the upper limit of the plot. For temp and hum models, this can be set to 2500. Rainfall model should set this to 4000.
plotSeasonalAverages = function(model,maxz) {
  library(raster)
  library(RColorBrewer)
  
  agg = matrix(nrow=4,ncol=0)
  titles=c("Mar - May", "Jun - Aug", "Sep - Nov", "Dec - Feb")
  load('ValidCells.RData')
  valid_cells = as.vector(na.omit(as.vector(t(valid_cells_matrix))))
  rm(valid_cells_matrix)
  
  colours=brewer.pal(9,"YlOrRd")[2:9]
  palette=colorRampPalette(colours)
  
  for(i in 1:25) {
    load(paste(model, '/', model, i,'.RData',sep=''))
    temp<-get(paste(model,i,sep=''))
    rm(list=paste(model,i,sep=''))
    means=rbind(
      # Spring
      colMeans(rbind(temp[330:365,], temp[1:56,])),
      # Summer
      colMeans(temp[57:148,]),
      # Autumn
      colMeans(temp[149:239,]),
      # Winter
      colMeans( temp[240:329,])
    )
    agg = cbind(agg, means)
  }
  for(index in 1:4) {
    blankRaster = raster('blankNew.gri')
    blankRaster[valid_cells] = agg[index,]
    png(filename = paste(model,'Agg',index,'.png',sep=''), height=720, width=949)
    plot(blankRaster,col=palette(100),axes=FALSE,zlim=c(0,maxz))
    title(titles[index], cex.main=2.5)
    dev.off()
  }
}
