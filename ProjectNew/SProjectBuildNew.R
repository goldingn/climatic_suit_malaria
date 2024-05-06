library(sp)
library(raster)
library(rgdal)

load("HProjectNew.RData")

rainf=stack("wc2.0_2.5m_prec/wc2.0_2.5m_prec_03.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_04.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_05.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_06.tif",
            "wc2.0_2.5m_prec/wc2.0_2.5m_prec_07.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_08.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_09.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_10.tif",
            "wc2.0_2.5m_prec/wc2.0_2.5m_prec_11.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_12.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_01.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_02.tif",
            "wc2.0_2.5m_prec/wc2.0_2.5m_prec_03.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_04.tif","wc2.0_2.5m_prec/wc2.0_2.5m_prec_05.tif")
rainf<-crop(rainf,extent(rainf,0,3600,0,8640))
alt=raster('alt2-5/alt.bil')

srainval<-getValues(rainf)
rm(rainf)
srainval<-srainval[hrel,]
srainval<-t((t(srainval/12))/c(31,30,31,30,31,31,30,31,30,31,31,28,31,30,31))

altval<-getValues(alt)
altval=altval[hrel]

load("HInput.RData")
ssrawrelval<-cbind(shrawrelval,srainval,altval,hrel)

save(ssrawrelval,file="SInputNew.RData")
