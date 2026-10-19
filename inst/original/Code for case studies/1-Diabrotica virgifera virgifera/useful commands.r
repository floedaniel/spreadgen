library(sp)
library(raster)
library(rgdal)
climexcsv=F
elevmax=F
habitatfile=T
habitatformat=".tif"
# CLIMEX results
plotRA() # risk area, where EI>0 and host is present 
plotGI()
plothabitat()
# model A
N0 = 1/3104*100
res=lgecon(N0=N0,r=0.33,t=5,econraster=T,econformat=".tif",mult=1,sim=10)
# model B
res = radial(RR=80,t=5,coord=c(20.30,44.82)) 
# hybrid model:
res = radialrand(RR=80,t=5,N0=N0,r=0.33,coord=c(20.30,44.82))
# model C
res = slg(N0=100*100/6.3E10,lmax=40,t=15)
# model D
res = dispk(N0=NULL, t=5, lmax=40, p=5, u=80, presencefile=T)

