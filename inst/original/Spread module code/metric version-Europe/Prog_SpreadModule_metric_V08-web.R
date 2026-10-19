
#-------------------------------------------------------------------------------------------
# Description: Spread module (task 2.5.3) for the EU project PRATIQUE
# Authors: Christelle Robinet (INRA, France) and Wopke van der Werf (WU, The Netherlands)
# Start: 19 March 2010
# Last modification: 30 May 2011 
#-------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LIST OF THE MAIN VARIABLES USED IN THIS CODE
#
# EI		: ecoclimatic index from CLIMEX
# GI		: growth index from CLIMEX
# xygrid	: y and x coordinates in decimal degree of the points from CLIMEX
# nxy		: the number of CLIMEX points
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("_________________________________________________________________________________________________________________________________")
print("              Version 8, metric version for Europe, May 30th 2011                                                                ")
print("              LAEA projection system at 10 km resolution                                                                         ")
print("                                                                                                                                 ")
print("                                                                                                                                 ")
print(" CITATIONS:                                                                                                                      ")
print(" Robinet et al. 2012, A suite of models to support the quantitative assessment of spread in pest risk analysis, PLoS ONE         ")
print(" Kehlenbeck et al. 2012, Modelling and mapping spread in pest risk analysis: a generic approach, Bull OEPP/EPPO Bull 42: 74-80   ")
print("                                                                                                                                 ")
print("_________________________________________________________________________________________________________________________________")
print("                                                                                                                                 ")


# points delimiting Europe
europeL <- read.table("EuropeLAEA.txt",header=F)

#Create projections
projWGS84="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
projETRS89="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
prj=as.character(expression(PROJCS["ETRS_1989_LAEA",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],UNIT["Meter",1.0]]))

#Extent raster for the Climex model
extrast = raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs=projWGS84)
res(extrast)=0.5

#Create a mask for the ouput file
mask=raster(nrows=500, ncols=450, xmn=2500000, xmx=7500000, ymn=1000000, ymx=5500000, crs=projETRS89)
res(mask)=10000

# grid for the models
xvect <- seq(2500000,7500000,by=10000)/1000
yvect <- seq(1000000,5500000,by=10000)/1000
xygrid <- expand.grid(xvect,yvect)
nxy <- length(xygrid[,1])


############# PROJECTION OF THE CLIMEX FILE ######################

if (climexcsv)
{
r = read.table("ClimexOutput.csv",sep=",",header=T,as.is=T,skip=7)

xyz = data.frame(x=r[,"Longitude"], y=r[,"Latitude"], z=r[,"EI"])
coordinates(xyz) = ~x+y
rEI=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
EI=projectRaster(rEI,mask)
EI <- xyValues(EI,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
exprNAEI <- is.na(EI)
EI[exprNAEI] <- 0


xyz = data.frame(x=r[,"Longitude"], y=r[,"Latitude"], z=r[,"GI"])
coordinates(xyz) = ~x+y
rGI=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
GI=projectRaster(rGI,mask)
GI <- xyValues(GI,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
GI <- GI*100/max(GI,na.rm=T) # we rescale GI values
GI[is.na(GI)] <- 0
}

if(!climexcsv)
{
# Re-worked file exported from CLIMEX
# Latitude, Longitude, EI, GI
tab <- read.table("ClimexOutput.txt",header=F)

xyz = data.frame(x=tab[,2], y=tab[,1], z=tab[,3])
coordinates(xyz) = ~x+y
rEI=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
EI=projectRaster(rEI,mask)
EI <- xyValues(EI,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
exprNAEI <- is.na(EI)
EI[exprNAEI] <- 0


xyz = data.frame(x=tab[,2], y=tab[,1], z=tab[,4])
coordinates(xyz) = ~x+y
rGI=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
GI=projectRaster(rGI,mask)
GI <- xyValues(GI,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
GI <- GI*100/max(GI,na.rm=T) # we rescale GI values
GI[is.na(GI)] <- 0
}



############# ELEVATION AND HOST RASTERS ######################

if(is.numeric(elevmax))
{
elevminr  <- raster("ZMIN.tif")
elevgridmin <- xyValues(elevminr,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
elevgridmin[is.na(elevgridmin)]<-0
elevmaxr  <- raster("ZMAX.tif")
elevgridmax <- xyValues(elevmaxr,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
elevgridmax[is.na(elevgridmax)]<-0
elevgrid <- (elevgridmin+elevgridmax)/2

# correct EI and GI according to the maximum elevation limit entered by the end user
EI <- EI*(elevgrid<elevmax)
GI <- GI*(elevgrid<elevmax)
}


if (habitatfile)
{
host <- raster(paste("habitat",habitatformat,sep=""))
hostgrid <- xyValues(host,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
exprNAhost <- is.na(hostgrid)
exprNATOT <- (exprNAEI)|(exprNAhost) 
hostgrid[is.na(hostgrid)]<-0

# correct EI according to host distribution
EI <- EI*(hostgrid>0)
GI <- GI*(hostgrid>0)
}
if (!habitatfile)
{
exprNATOT <- (exprNAEI)
} 

tab <- cbind(xygrid,EI,GI)  # be careful: longitude, latitude, EI and GI (lat long EI GI in non metric models)

#############################################################

plotRA <- function()
{
	par(bg = "lightgrey")
	mat <- matrix(EI,ncol=length(yvect),byrow=F)
	matNA <- matrix(exprNATOT,ncol=length(yvect),byrow=F)
	mat[matNA] <- (-9999)
	image(xvect,yvect,mat,breaks=c(-10000,-1,0,101),col=c("lightgrey","white","orange"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Area of potential establishment",asp=1)
	points(europeL,pch=".")
	par(bg = "transparent")
	print(" LEGEND:                                            ")
	print(" orange dots: EI>0, white dots: EI=0, grey: no data ")
}


plotGI <- function()
{
	par(bg = "lightgrey")
	mat <- matrix(GI,ncol=length(yvect),byrow=F)
	matNA <- matrix(exprNATOT,ncol=length(yvect),byrow=F)
	mat[matNA] <- (-9999)
	coulpal <- colorRampPalette(c("white","yellow","orange","red"),space="rgb")
	coul <- coulpal(11)
	image(xvect,yvect,mat,breaks=c(-10000,-1,0,10,20,30,40,50,60,70,80,90,100),col=c("lightgrey","white",coul[2],coul[3],coul[4],coul[5],coul[6],coul[7],coul[8],coul[9],coul[10],coul[11]),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Growth potential",asp=1)
	points(europeL,pch=".")
	par(bg = "transparent")
	print(" LEGEND:                                                                                          ")
	print(" graduated colors from white (GI=0) to red (GI=100) for the rescaled GI. Grey means no data")
}

plothabitat<- function()
{
	if(length(unique(hostgrid))<3)
	{
	par(bg = "lightgrey")
	mat <- matrix(hostgrid>0,ncol=length(yvect),byrow=F)
	matNA <- matrix(exprNATOT,ncol=length(yvect),byrow=F)
	mat[matNA] <- (-9999)
	image(xvect,yvect,mat,breaks=c(-10000,-1,0,101),col=c("lightgrey","white","green"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Habitat distribution",asp=1)
	points(europeL,pch=".")
	print(" LEGEND:                                                          ")
	print(" green dots: habitat present, white dots: habitat absent, grey: no data ")
	par(bg = "transparent")
	}
	if(length(unique(hostgrid))>2)
	{
	par(bg = "lightgrey")
	mat <- matrix(hostgrid,ncol=length(yvect),byrow=F)
	matNA <- matrix(exprNATOT,ncol=length(yvect),byrow=F)
	mat[matNA] <- (-9999)
	coulpal <- colorRampPalette(c("white","lightgreen","darkgreen"),space="rgb")
	coul <- coulpal(11)
	seqh <- seq(from=min(hostgrid,na.rm=T),to=max(hostgrid,na.rm=T),length.out=10)
	if(seqh[1]>0)
	{
	image(xvect,yvect,mat,breaks=c(-10000,-1,seqh[1],seqh[2],seqh[3],seqh[4],seqh[5],seqh[6],seqh[7],seqh[8],seqh[9],seqh[10]),col=c("lightgrey",coul[2],coul[3],coul[4],coul[5],coul[6],coul[7],coul[8],coul[9],coul[10],coul[11]),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Habitat distribution",asp=1)
	}
	if(seqh[1]==0)
	{
	image(xvect,yvect,mat,breaks=c(-10000,-1,0,seqh[2],seqh[3],seqh[4],seqh[5],seqh[6],seqh[7],seqh[8],seqh[9],seqh[10]),col=c("lightgrey","white",coul[3],coul[4],coul[5],coul[6],coul[7],coul[8],coul[9],coul[10],coul[11]),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Habitat distribution",asp=1)
	}
	points(europeL,pch=".")
	par(bg = "transparent")
	print("LEGEND:                                                                                ")
	nom <- paste("graduated colors from white (no habitat) to dark green (max=",round(max(hostgrid,na.rm=T),2),"). Grey means no data",sep="")
	print(nom)
	}
}


############################################################

# values needed to parametrize the models:

printinfo <- function()
{
print("                                                          ")
print("-------------  For this specific case study: -------------")

resol <- paste("The grid resolution is 10 km")
print(resol)

calarea <- paste("The area of a cell is 100 km²")
print(calarea)

calnRA <- paste("The number of suitable cells is: ",sum(EI>0),sep="")
print(calnRA)

calntot <- paste("The total number of cells is: ",length(EI),sep="")
print(calntot)

calNAclimex <- paste("The number of missing values in the CLIMEX dataset is: ",sum(exprNAEI),sep="")
print(calNAclimex)

if(habitatfile)
{
calNAhost <- paste("The number of missing values in the habitat distribution is: ",sum(exprNAhost),sep="")
print(calNAhost)

calNATOT <- paste("The total number of missing values (CLIMEX or host) is: ",sum(exprNATOT),sep="")
print(calNATOT)
}

print("                                                          ")
print("-----  Type the commands below to call the functions -----")
print("                                                          ")
}

printinfo()

##############################################
# SIMPLE LOGISTIC GROWTH: model C				
##############################################

slg <- function(N0=1,lmax=7,movie=F,t=1,export=F,name="SLGoutput")
{
	growth <- GI*log(lmax)/100
	val <- 100*N0*exp(growth*t)/(100+N0*(exp(growth*t)-1)) 
	expr0 <- (EI==0)
	res<- (!expr0)*val # percentage of the cell invaded

	if(!movie)
	{
	nom <- paste("SLG model t = ",t,sep="")
	mat <- matrix(res,ncol=length(yvect),byrow=F)
	matNA <- matrix(exprNATOT,ncol=length(yvect),byrow=F)
	mat[matNA] <- (-9999)
	image(xvect,yvect,mat,breaks=c(-10000,-1,0,25,50,75,100),col=c("white","lightgrey","blue","green","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main=nom,asp=1)
	points(europeL,pch=".")
	}

	if(movie)
	{
	for (k in 1:t)
	{
	if(k<t){nom <- paste("SLG model t = ",k," (CLICK ON THE MAP)",sep="")}
	if(k==t){nom <- paste("SLG model t = ",k,sep="")}
	val <- 100*N0*exp(growth*k)/(100+N0*(exp(growth*k)-1)) 
	res<- (!expr0)*val # percentage of the cell invaded
	mat <- matrix(res,ncol=length(yvect),byrow=F)
	matNA <- matrix(exprNATOT,ncol=length(yvect),byrow=F)
	mat[matNA] <- (-9999)
	image(xvect,yvect,mat,breaks=c(-10000,-1,0,25,50,75,100),col=c("white","lightgrey","blue","green","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main=nom,asp=1)
	points(europeL,pch=".")
	if(k<t){locator(1)}
	}

	}

	if(export)
	{
	res[exprNATOT] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=res)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}

print("LEGEND:                                               ")
print("Nt: population abundance (%)                          ")
print("grey: Nt=0, blue:  0 < Nt < = 25, green: 25 < Nt < = 50, orange: 50 < Nt < = 75, red: 75 < Nt < = 100, white: no data")
return(list(slg=res,ntot=length(EI),nRA=sum(!expr0),sum0=sum(mat==0,na.rm=T),sum25=sum((mat>0)&(mat<=25),na.rm=T),sum50=sum((mat>25)&(mat<=50),na.rm=T),sum75=sum((mat>50)&(mat<=75),na.rm=T),sum100=sum((mat>75)&(mat<=100),na.rm=T)))
}


####################################################
# LOGISTIC GROWTH BASED ON ECONOMIC VALUES: model A		
####################################################

lgecon <- function(N0=1,r=1,t=1,econraster=F,mult=1,sim=1,habitatraster=F,valperhost=NULL,econformat=NULL,export=F,name="LGECONoutput")
{

if(econraster & habitatraster){print("BE CAREFUL: econraster=T and habitatraster=T. The code used the econraster file for the calculations")}

if (econraster)
{
	econ <- raster(paste("econ",econformat,sep=""))
	econgrid <- xyValues(econ,cbind(xygrid[,1]*1000,xygrid[,2]*1000))
	exprEconNA <- is.na(econgrid)
	if(sum(exprEconNA)>0){print("CAUTION: the model excludes areas where economic values are missing")}
	econgrid[is.na(econgrid)] <- 0     # where there is missing data, we assume that the value is 0 
	econgrid <- econgrid*mult	   # to convert in the right unit
	breakvalues <- seq(0,max(econgrid),length.out=6)
	tmp <- paste("The number of missing values in the econ file is: ",sum(exprEconNA),sep="")
	print(tmp)
}

if (habitatraster & (!econraster))
{
	if(is.null(valperhost)){print("ERROR: valperhost, value per km² if habitat everywhere, is required in this case")}
	econgrid <- hostgrid*valperhost
	econgrid <- econgrid*mult*(res(mask)[1]/1000)*(res(mask)[2]/1000)   
	breakvalues <- seq(0,max(econgrid),length.out=6)
	exprEconNA <- is.na(econgrid)
}

if(!(econraster|habitatraster)) 
{
	econgrid <- runif(n=length(EI),min=1,max=100)
	econgrid <- econgrid*mult*(res(mask)[1]/1000)*(res(mask)[2]/1000)  
	breakvalues <- seq(0,max(econgrid),length.out=6)
	print("CAUTION: no economic file was provided => economic values are taken at random")
	exprEconNA <- is.na(econgrid)
}
	
nRA <- sum((EI>0)&(!exprEconNA)) # number of pixels included in the risk area (where EI>0)
val <- 100*N0*exp(r*t)/(100+N0*(exp(r*t)-1)) # is the percentage of assets invaded
ninv <- round(val*nRA/100) # number of pixels invaded 

	vectbest <- (EI<0); vectbest[EI==0] <- (-1); vectbest[exprNATOT|exprEconNA] <- (-9999) 	# 0 -non invaded- everywhere at the beginning,  -1 for non suitable (EI=0) and then -9999 for missing data
	vectworst <- vectbest
	vectrand <- vectbest 

	# the code is going to put 1 where the cell is invaded. 
	nID <- 1:nxy    	
	nIDselect <- nID[(EI>0)&(!exprEconNA)]

if(ninv>0)
{
	par(mfrow=c(3,2))
	# selection of the points with the risk area and where we have economic data:
	ecovalTOT <- econgrid[(EI>0)&(!exprEconNA)]

# Best case scenario: pixels with low values are invaded first
#--------------------------------------------------------------
orderbest <- order(ecovalTOT,decreasing=FALSE)
ecobest <- ecovalTOT[t(orderbest)]; resbest <- sum(ecobest[1:ninv])	# sum of economic values in invaded cells
nIDbest <- nIDselect[t(orderbest)]; nIDbest <- nIDbest[1:ninv]		# ID of the lines of invaded cells
vectbest[nIDbest] <- 1							# update the file for invaded cells
	mat <- matrix(vectbest,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-100,-1,0,1),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Best case",asp=1)
	points(europeL,pch=".")
	hist(ecobest[1:ninv],main="Best case",xlab="Economic value (euros/cell)",breaks=breakvalues)

# Medium case scenario: random selection of the pixels
#-----------------------------------------------------

# we make some replicate simulations:
memreal <- NULL
mem <- NULL
for (i in 1:sim)
{
ni <- sample(1:length(ecovalTOT),size=ninv,replace=F)
resmed <- sum(ecovalTOT[ni])
nIDrand <- nIDselect[ni]
vectrand[nIDrand] <- 1	
memreal <- rbind(memreal,cbind(i,resmed)) # keep the total value for each simulation
mem <- rbind(mem,rbind(t(t(ecovalTOT[ni])))) # keep all the values 
}

	mat <- matrix(vectrand,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-100,-1,0,1),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Random case",asp=1)
	points(europeL,pch=".")
	nom <- paste("Random case - ",sim," replicate simulations",sep="")
	hist(mem,main=nom,xlab="Economic value (euros/cell)",breaks=breakvalues)
	resrand <- summary(memreal[,2])

# Worst case scenario: pixels with high values are invaded first
#-----------------------------------------------------------------
orderworst <- order(ecovalTOT,decreasing=TRUE)
ecoworst <- ecovalTOT[t(orderworst)]; resworst <- sum(ecoworst[1:ninv])		# sum of economic values in invaded cells
nIDworst <- nIDselect[t(orderworst)]; nIDworst <- nIDworst[1:ninv]		# ID of the lines of invaded cells
vectworst[nIDworst] <- 1							# update the file for invaded cells
	mat <- matrix(vectworst,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-100,-1,0,1),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Worst case",asp=1)
	points(europeL,pch=".")
	hist(ecoworst[1:ninv],main="Worst case",xlab="Economic value (euros/cell)",breaks=breakvalues)

par(mfrow=c(1,1))
}

if(ninv==0)
{ 
	mat <- matrix(vectbest,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-100,-1,0,1),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="LG-Econ",asp=1)
	points(europeL,pch=".")
	resworst=0
	resbest=0
	resrand=0
}
	if(export)
	{
	vectbest[exprNATOT|exprEconNA] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=vectbest)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,"_best.tiff",sep=""),format="GTiff")

	vectworst[exprNATOT|exprEconNA] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=vectworst)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,"_worst.tiff",sep=""),format="GTiff")

	vectrand[exprNATOT|exprEconNA] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=vectrand)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,"_rand.tiff",sep=""),format="GTiff")

	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}

print("LEGEND:                                          ")
print("red: invaded cells, orange: non invaded but suitable cells, grey: EI=0, white: no data")
return(list(nRA=nRA,ninv=ninv,ntot=dim(xygrid)[1],pniche=val,worst=resworst,best=resbest,rand=resrand,econvalues=range(econgrid)))
}



##############################################
# RADIAL RANGE EXPANSION: model B   				
##############################################


radial <- function(RR=20,t=1,coord=NULL,export=F,name="radialOutput")
{
	expr0 <- (EI==0)

# ENTRY POINT:
#-------------
# if random selection of entry point:
if(is.null(coord))
{
entry <- sample(1:sum(!expr0),size=1)
tmp <- xygrid[!expr0,]
X <- tmp[entry,1]; Xtot <- X
Y <- tmp[entry,2]; Ytot <- Y
exprentry <- (xygrid[,2]==Y)&(xygrid[,1]==X)
ientry <- (1:nxy)[exprentry]
}

# if entry point provided:
if(!is.null(coord))
{
Xtot <- NULL
Ytot <- NULL
nbpoint <- length(coord)/2
ientry <- NULL
for (i in 1:nbpoint)
{
x <- coord[2*i-1]	# longitude, not necessarily from the grid cell=> will find the closest grid cell
y <- coord[2*i]		# latitude,  not necessarily from the grid cell=> will find the closest grid cell
tmp <- xygrid[!expr0,]
expr <- (min((tmp[,1]-x)^2+(tmp[,2]-y)^2)==(tmp[,1]-x)^2+(tmp[,2]-y)^2)
X <- tmp[expr,1]
Y <- tmp[expr,2]
if(sum(expr)>1){X <- X[1]; Y <- Y[1]}
exprentry <- (xygrid[,2]==Y)&(xygrid[,1]==X)
tmpientry <- (1:nxy)[exprentry] 
ientry <- rbind(ientry,tmpientry)
Xtot <- rbind(Xtot,X)
Ytot <- rbind(Ytot,Y)
}
}

# calculation for time t
#------------------------
	vect <- (EI<0); vect[EI==0] <- (-1); vect[exprNATOT] <- (-9999) 	# 0 -non invaded- everywhere at the beginning,  -1 for non suitable (EI=0) and then -9999 for missing data
	vect[ientry] <- 1	# the code is going to put 1 where the cell is invaded. 
	nID <- 1:nxy    	# ID of the line 
	ninv <- nID[vect==1] 
	nnon <- nID[vect==0]	
	vect2 <- vect
	for (i in ninv) # for all cells invaded
	{
	R <- sqrt( (xygrid[i,1]-xygrid[nnon,1])^2 + (xygrid[i,2]-xygrid[nnon,2])^2 )
	expr <- (R <= RR*t)
	ninvaded <- nnon[expr]
	vect2[ninvaded] <- 1
	}
	mat <- matrix(vect2,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-100,-1,0,1),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Radial expansion",asp=1)
	points(europeL,pch=".")
	points(Xtot,Ytot,pch=19)

	if(export)
	{
	vect2[exprNATOT] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=vect2)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}

pniche=100*sum(mat==1)/sum(!expr0)
print("LEGEND:                                                   ")
print("red: invaded cells, orange: non invaded but suitable cells, grey: EI=0, white: no data")
return(list(radial=vect2,ninv=sum(mat==1),nRA=sum(!expr0),ntot=length(EI),pniche=pniche))
}



##############################################
# COMBINED MODEL: model A random + model B			
##############################################

radialrand <- function(N0=1,r=1,t=1,RR=20,coord=NULL,export=F,name="hybridOutput")
{
	expr0 <- (EI==0)

# ENTRY POINT:
#-------------
# if random selection of entry point:
if(is.null(coord))
{
entry <- sample(1:sum(!expr0),size=1)
tmp <- xygrid[!expr0,]
X <- tmp[entry,1]
Y <- tmp[entry,2]
exprentry <- (xygrid[,2]==Y)&(xygrid[,1]==X)
ientry <- (1:nxy)[exprentry]
Xtot <- X
Ytot <- Y
}

# if entry point provided:
if(!is.null(coord))
{
Xtot <- NULL
Ytot <- NULL
nbpoint <- length(coord)/2
ientry <- NULL
for (i in 1:nbpoint)
{
x <- coord[2*i-1]	# longitude, not necessarily from the grid cell=> will find the closest grid cell
y <- coord[2*i]		# latitude,  not necessarily from the grid cell=> will find the closest grid cell
tmp <- xygrid[!expr0,]
expr <- (min((tmp[,1]-x)^2+(tmp[,2]-y)^2)==(tmp[,1]-x)^2+(tmp[,2]-y)^2)
X <- tmp[expr,1]
Y <- tmp[expr,2]
if(sum(expr)>1){X <- X[1]; Y <- Y[1]}
exprentry <- (xygrid[,2]==Y)&(xygrid[,1]==X)
tmpientry <- (1:nxy)[exprentry] 
ientry <- rbind(ientry,tmpientry)
Xtot <- rbind(Xtot,X)
Ytot <- rbind(Ytot,Y)
}
}

# calculation for time t
#------------------------
	vect <- (EI<0); vect[EI==0] <- (-1); vect[exprNATOT] <- (-9999) 	# 0 -non invaded- everywhere at the beginning,  -1 for non suitable (EI=0) and then -9999 for missing data
	vect[ientry] <- 1	# the code is going to put 1 where the cell is invaded. 
	nID <- 1:nxy    	# ID of the line 
	ninv <- nID[vect==1] 
	nnon <- nID[vect==0]	
	vect2 <- vect
	for (i in ninv) # for all cells invaded
	{
	R <- sqrt( (xygrid[i,1]-xygrid[nnon,1])^2 + (xygrid[i,2]-xygrid[nnon,2])^2 )
	expr <- (R <= RR*t)
	ninvaded <- nnon[expr]
	vect2[ninvaded] <- 1
	}

# random selection of the pixels
#-----------------------------------------------------
	nRA <- sum(EI>0) # number of pixels included in the risk area (where EI>0)
	val <- 100*N0*exp(r*t)/(100+N0*(exp(r*t)-1)) # is the percentage of assets invaded
	ninv <- round(val*nRA/100) # number of pixels invaded 
	if (ninv>sum(vect2==1)){vect3 <- vect2}
	if (ninv<=sum(vect2==1))
	{
	ni <- sample(1:sum(vect2==1),size=ninv,replace=F)
	nvect2 <- nID[vect2==1]
	vect3 <- vect
	vect3[nvect2[ni]] <- 1
	}
	
	mat <- matrix(vect3,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-100,-1,0,1),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Radial random expansion",asp=1)
	points(europeL,pch=".")
	points(Xtot,Ytot,pch=19)

	if(export)
	{
	vect3[exprNATOT] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=vect3)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}

pniche=100*sum(vect3==1)/sum(!expr0)
print("LEGEND:                                                   ")
print("red: invaded cells, orange: non invaded but suitable cells, grey: EI=0, white: no data")
return(list(radialrand=vect3,ninv=sum(vect3==1),nRA=sum(!expr0),ntot=length(EI),pniche=pniche))
}


####################################################
# DISPERSAL KERNEL: DETERMINISTIC VERSION - model D
####################################################

plotkernel <- function(vectpop0,presence=NULL,t=NULL,threshold=NULL,legend=T)
{
if(is.null(threshold))
{
	vect <- vectpop0
	vect[EI==0] <- -1
	vect[exprNATOT] <- (-9999)
	mat <- matrix(vect,ncol=length(yvect),byrow=F)
	coulpal <- colorRampPalette(c("yellow","orange","red"),space="rgb")
	coul <- coulpal(9)
	if(!is.null(t)){nom <- paste("Dispersal kernel for time t = ",t,sep="")}
	if(is.null(t)){nom <- paste("Dispersal kernel")}
	image(xvect,yvect,mat,breaks=c(-10000,-10,-1,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100),col=c("white","lightgrey",coul[1],coul[2],coul[3],coul[4],coul[5],coul[6],coul[7],coul[8],coul[9]),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main=nom,asp=1)
	points(europeL,pch=".")
	points(xygrid[presence>0,1],xygrid[presence>0,2],col="black",pch=19)

if(legend)
{
print("LEGEND:                                                                ")
print("Nt: population abundance (%)                                           ")
print(" graduated colors from white (Nt < 10^(-6) %) to yellow, orange and red (Nt > = 10%). Grey means no data")
}
}

if(!is.null(threshold))
{
	nom <- paste("Dispersal kernel for a threshold Nt = ",threshold,sep="")
	vect <- vectpop0
	vect[EI==0] <- -1
	vect[exprNATOT] <- (-9999)
	mat <- matrix(vect,ncol=length(yvect),byrow=F)
	image(xvect,yvect,mat,breaks=c(-10000,-10,-1,threshold,100),col=c("white","lightgrey","orange","red"),xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main=nom,asp=1)
	points(europeL,pch=".")
	if(!is.null(presence)){points(xygrid[presence>0,1],xygrid[presence>0,2],col="black",pch=19)}

if(legend)
{
print("LEGEND:                                                                ")
print("Nt: population abundance (%)                                           ")
print("red dots: cells where Nt > threshold, orange: cells where Nt<= threshold, grey: unsuitable cells (EI=0), white: no data")
}
}
}

exportkernel <- function(vectpop0,name="kernelOutput")
{
	vectpop0[exprNATOT] <- NA
	xyz = data.frame(x=1000*xygrid[,1], y=1000*xygrid[,2], z=vectpop0)
	coordinates(xyz) = ~x+y
	r=pointsToRaster(mask,xyz,values=xyz$z,fun=mean)
	writeRaster(r,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
}





dispk <- function(N0=NULL, t=1,lmax=7, P=1, p=50, u=100,presencefile=T,nentry=NULL)  
{

#--------------------------------------------
if(!presencefile)
{# case 3
entry <- sample(1:sum(!(EI==0)),size=nentry)
tmp <- xygrid[!(EI==0),]
xentry <- tmp[entry,1]
yentry <- tmp[entry,2]
presence <- cbind(xentry,yentry)
presencegrid <- rep(F,length(xygrid[,1]))
for (i in 1:nentry)
{
	expr <- (min( ((presence[i,1]-xygrid[,1])^2)+((presence[i,2]-xygrid[,2])^2) ) == ( ((presence[i,1]-xygrid[,1])^2)+((presence[i,2]-xygrid[,2])^2) ))
	presencegrid[expr] <- T
}
presence <- presencegrid
# vector of TRUE/FALSE
}

if(presencefile)
{
presence <- read.table("presence.txt",header=F)

ncolumn <- dim(presence)[2]

if(ncolumn==2)
{ # case 2
presencegrid <- rep(F,length(xygrid[,1]))
for (i in 1:dim(presence)[1])
{
	expr <- (min( ((presence[i,1]-xygrid[,1])^2)+((presence[i,2]-xygrid[,2])^2) ) == ( ((presence[i,1]-xygrid[,1])^2)+((presence[i,2]-xygrid[,2])^2) ))
	presencegrid[expr] <- T
}
presence <- presencegrid
# vector of TRUE/FALSE
}

if(ncolumn==3)
{ # case 2
presencegrid <- rep(0,length(xygrid[,1]))
count  <- rep(0,length(xygrid[,1]))
for (i in 1:dim(presence)[1])
{
	expr <- (min( ((presence[i,1]-xygrid[,1])^2)+((presence[i,2]-xygrid[,2])^2) ) == ( ((presence[i,1]-xygrid[,1])^2)+((presence[i,2]-xygrid[,2])^2) ))
	presencegrid[expr] <- presencegrid[expr]+presence[i,3]
	count[expr] <- count[expr]+1
}
presence <- presencegrid
presence[!presence==0] <- presence[!presence==0]/count[!presence==0] # we calculate the mean population density in case there is more than one presence point close to the grid point 
# vector of population density
}

}
#------------------------------

if(is.logical(presence)&(is.null(N0))){print("You should enter a value N0")}
if((!is.logical(presence))&(!is.null(N0))){print("The model used the N0 values provided in the file presence.txt and not the value given in the function")}


if(is.logical(presence)){presence <- as.numeric(presence)*N0}
vectpop0 <- presence



for (k in 1:t)
{
print(c("t = ",k))
popdisp <- rep(0,length(xygrid[,1]))



print(paste("Total number of iterations= ",sum(!(EI==0)),sep=""))
for (i in (1:length(xygrid[,1]))[(!(EI==0))|(presence>0)]) # we assume that a source point is necessarilly where the EI>0 (to reduce the potential source points and calculations)
{
 
	x0 <- xygrid[i,1]
	y0 <- xygrid[i,2]

	R <- sqrt( (x0-xygrid[,1])^2 + (y0-xygrid[,2])^2 )	
	f <- (1/((u^2)*p*pi)) * gamma((p+1)/2)/( gamma((p-1)/2)*(1+(R^2)/(p*(u^2)))^((p+1)/2) )
	probaf <- f*((res(mask)[1]/1000))*((res(mask)[2]/1000))

	popdisp <- popdisp + probaf*P*vectpop0[i] # we consider only individuals engaged in dispersal
	popdisp[EI==0] <- 0 # to prevent dispersal in cells where EI=0
}


popdisp <- popdisp + (1-P)*vectpop0 # we add individuals staying at the same place

popdisp[EI==0] <- 0 # we assume that individuals die where the EI=0

growth <- GI*log(lmax)/100
vectpop0 <- (100*popdisp*exp(growth)/(100+popdisp*(exp(growth)-1)))

plotkernel(vectpop0=vectpop0,presence=presence,t=k,legend=F)
}

print("LEGEND:                                                                ")
print("Nt: population abundance (%)                                           ")
print(" graduated colors from white (Nt < 10^(-6) %) to yellow, orange and red (Nt > = 10%). Grey means no data")

# percentage of the niche invaded:
pniche <- (100*sum(vectpop0>1)/sum(!(EI==0)))

return(list(dispk=vectpop0,pniche=pniche,ntot=length(vectpop0),nRA=sum(!(EI==0)),ninv= sum(vectpop0>1),presence=presence))
}




plotkp <- function(u,p,xmax=1000,add=F,color="red")
{
	R <- seq(0,xmax,by=1)
	f <- (1/((u^2)*p*pi)) * gamma((p+1)/2)/( gamma((p-1)/2)*(1+(R^2)/(p*(u^2)))^((p+1)/2) )
	if(!add)
	{
	plot(R,f,type="l",lwd=2,xlab="Distance (km)",ylab="Probability",cex.axis=1.5,cex.lab=1.5,col=color)
	}
	if(add)
	{
	lines(R,f,lwd=2,col=color)
	}
	
}









##########################################################################
# supplementary function(s)



proj <- function(x,y)  # function to transform the coordinates in decimal degrees into the LAEA projection in km
{
xy=matrix(c(x,y),1)
return(project(xy,projETRS89)/1000)
}

projinv <- function(x,y)  # function to transform the coordinates in LAEA projection in km into decimal degrees
{
xy=matrix(c(1000*x,1000*y),1)
return(project(xy,projETRS89,inv=T))
}



pointtocell <- function(name)
{
listpoints <- read.table(name,header=F)
Xtot <- NULL
Ytot <- NULL
nbpoint <- dim(listpoints)[1]
for (i in 1:nbpoint)
{
x <- listpoints[i,1]	
y <- listpoints[i,2]
tmp <- xygrid[EI>0,]
expr <- (min((tmp[,1]-x)^2+(tmp[,2]-y)^2)==(tmp[,1]-x)^2+(tmp[,2]-y)^2)
X <- tmp[expr,1]
Y <- tmp[expr,2]
if(sum(expr)>1){X <- X[1]; Y <- Y[1]}
if(is.null(Xtot))
{
	Xtot <- rbind(Xtot,X)
	Ytot <- rbind(Ytot,Y)
}
if(!is.null(Xtot))
{	if(min( (Xtot-X)^2+(Ytot-Y)^2) >0)
	{
	Xtot <- rbind(Xtot,X)
	Ytot <- rbind(Ytot,Y)
	}
}
}
return(length(Xtot))
}
