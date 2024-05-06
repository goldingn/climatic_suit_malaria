library(sp)
library(raster)
library(mgcv)
library(rgdal)
library(Rcpp)
library(tempsuitcalc)
library(parallel)

# load("ProjectNew.RData")
# load("InputNew.RData")

tmin_short <- raster::stack("Air_Temp_Min_5km.grd")
tmin <- tmin_short[[c(3:12, 1:5)]]
rm(tmin_short)

tmax_short <- raster::stack("Air_Temp_Max_5km.grd")
tmax <- tmax_short[[c(3:12, 1:5)]]

meanTemp = overlay(tmax, tmin, fun = function(x,y){ (x+y)/2 } )
writeRaster(meanTemp, filename="meanTemp.grd", overwrite=TRUE)
sinTemp = overlay(tmax, tmin, fun = function(x,y){ (x-y)/2 } )
writeRaster(sinTemp, filename="sinTemp.grd", overwrite=TRUE)

naRaster <- meanTemp[[1]]
for (i in 2:12) {
  naRaster <- raster::mask(naRaster, meanTemp[[i]])
}
blank<-naRaster
blank[!is.na(blank[])] <- 0
writeRaster(blank, filename="blankNew.grd",overwrite=TRUE)

basena = extract(calc(naRaster, function(x)(is.na(x))), 1:ncell(naRaster))
na = which(basena!=1)

maxN <- function(x, N=3){
  len <- length(x)
  return(sort(x,partial=len-N+1)[len-N+1])
}

findmax <- extract(tmax_short, na)
thirdhigh <- apply(findmax, 1, maxN)
rel=na[which(thirdhigh>=16)]

relcoordy=as.integer(floor(rel/8640))
relcoordraw=cbind(as.integer(rel-relcoordy*8640),relcoordy)
relcoord=matrix(nrow=nrow(relcoordraw),ncol=2)
relcoord[,1]<-(relcoordraw[,1]/24)-180-(1/144)
relcoord[,2]<-((3600-relcoordraw[,2])/24)-60-(1/144)

save(na,
     rel,
     relcoord,
     blank,
     meanTemp,
     sinTemp,
     
     file = "ProjectNew.RData")

smeanval<-round(getValues(meanTemp),2)
smeanval<-smeanval[rel,]
ssinval<-round(getValues(sinTemp),2)
ssinval<-ssinval[rel,]
srawrelval<-cbind(smeanval,ssinval)

save(srawrelval, file = "InputNew.RData")
