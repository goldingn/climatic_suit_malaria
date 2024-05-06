library(sp)
library(raster)
library(rgdal)

load("ProjectNew.RData")

vapr=stack("wc2-5/vapr_03.tif","wc2-5/vapr_04.tif","wc2-5/vapr_05.tif","wc2-5/vapr_06.tif","wc2-5/vapr_07.tif",
           "wc2-5/vapr_08.tif","wc2-5/vapr_09.tif","wc2-5/vapr_10.tif","wc2-5/vapr_11.tif","wc2-5/vapr_12.tif",
           "wc2-5/vapr_01.tif","wc2-5/vapr_02.tif","wc2-5/vapr_03.tif","wc2-5/vapr_04.tif","wc2-5/vapr_05.tif")
vapr<-crop(vapr,extent(vapr,0,3600,0,8640))

naRaster <- blank
for (i in 1:12) {
  naRaster <- raster::mask(naRaster, vapr[[i]])
}
hblank<-naRaster
hblank[!is.na(hblank[])] <- 0
writeRaster(hblank, filename="hblankNew.grd",overwrite=TRUE)

basena = extract(calc(naRaster, function(x)(is.na(x))), 1:ncell(naRaster))
hna = which(basena!=1)

hrel=intersect(hna,rel)

hrelcoordy=as.integer(floor(hrel/8640))
hrelcoordraw=cbind(as.integer(hrel-hrelcoordy*8640),hrelcoordy)
hrelcoord=matrix(nrow=nrow(hrelcoordraw),ncol=2)
hrelcoord[,1]<-(hrelcoordraw[,1]/24)-180-(1/144)
hrelcoord[,2]<-((3600-hrelcoordraw[,2])/24)-60-(1/144)

save(hna,
     hrel,
     hrelcoord,
     hblank,
     
     file = "HProjectNew.RData")


shumval<-getValues(vapr)
rm(vapr)
shumval<-shumval[hrel,]
meanTemp<-raster::stack("meanTemp.grd")
sinTemp<-raster::stack("sinTemp.grd")
smeanval<-round(getValues(meanTemp),2)
smeanval<-smeanval[hrel,]
ssinval<-round(getValues(sinTemp),2)
ssinval<-ssinval[hrel,]
shrawrelval<-cbind(smeanval,ssinval,shumval)

save(shrawrelval,file="HInputNew.RData")
