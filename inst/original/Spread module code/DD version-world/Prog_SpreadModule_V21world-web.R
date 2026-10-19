
#-------------------------------------------------------------------------------------------
# Description: Spread module (task 2.5.3) for the EU project PRATIQUE
# Authors: Christelle Robinet (INRA, France) and Wopke van der Werf (WU, The Netherlands)
# Start: 19 March 2010
# Last modification:  26 Aug 2011 
#-------------------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LIST OF THE MAIN VARIABLES USED IN THIS CODE
#
# tab		: output from CLIMEX (col1: latitude in DD, col2: longitude in DD, col3: EI, col4: GI
# EI		: ecoclimatic index from "tab"
# GI		: growth index from "tab"
# xygrid	: y and x coordinates in DD of the points from "tab"
# nxy		: the number of CLIMEX points

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


print("________________________________________________________________________________________________________________________________")
print("              Version 21, world version (in decimal degree), Aug 26th 2011                                                      ")
print("                                                                                                                                ")
print("                                                                                                                                ")
print(" CITATIONS:                                                                                                                     ")
print(" Robinet et al. 2012, A suite of models to support the quantitative assessment of spread in pest risk analysis, PLoS ONE        ")
print(" Kehlenbeck et al. 2012, Modelling and mapping spread in pest risk analysis: a generic approach, Bull OEPP/EPPO Bull 42: 74-80  ")
print("                                                                                                                                ")
print("________________________________________________________________________________________________________________________________")
print("                                                                                                                                ")



############# DOWNLOAD THE DATA ######################

# points delimiting the countries
folder="shpfile"
lyr="cntry06"
shp <- readOGR(dsn=folder,layer=lyr)



if(climexcsv==T)
{ # direct output of CLIMEX
r = read.table("ClimexOutput.csv",sep=",",header=T,as.is=T,skip=7)
tab <- cbind(r$Latitude, r$Longitude, r$EI,r$GI)
EI <- tab[,3]
exprNAEI <- is.na(EI)
EI[exprNAEI] <- 0
tab[exprNAEI,3] <- 0
GI <- tab[,4]
GI <- GI*100/max(GI,na.rm=T) # we rescale GI values
exprNAGI <- is.na(GI)
GI[exprNAGI] <- 0
tab[exprNAGI,4] <- 0
exprNATOT <- exprNAEI|exprNAGI
}

if(climexcsv==F)
{
# Re-worked file exported from CLIMEX
# Latitude, Longitude, EI, GI
tab <- read.table("ClimexOutput.txt",header=F)
EI <- tab[,3]
exprNAEI <- is.na(EI)
EI[exprNAEI] <- 0
tab[exprNAEI,3] <- 0
GI <- tab[,4]*100/max(tab[,4],na.rm=T) # we rescale GI values
exprNAGI <- is.na(GI)
GI[exprNAGI] <- 0
tab[exprNAGI,4] <- 0
exprNATOT <- exprNAEI|exprNAGI
}

# create the CLIMEX grid in decimal degrees:
xygrid <- tab[,1:2]
nxy <- length(xygrid[,1])

if(is.numeric(elevmax))
{
elev  <- raster("elev.asc")
projWGS84="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
projection(elev)=projWGS84
elevgrid <- xyValues(elev,cbind(xygrid[,2],xygrid[,1]))
# correct EI and GI according to the maximum elevation limit entered by the end user
EI <- EI*(elevgrid<elevmax); EI[is.na(EI)] <- 0
tab[,3] <- EI
GI <- GI*(elevgrid<elevmax); GI[is.na(GI)] <- 0
tab[,4] <- GI
exprNATOT <- exprNATOT|is.na(elevgrid)
}


if (habitatfile)
{
host <- raster(paste("habitat",habitatformat,sep=""))
hostgrid <- xyValues(host,cbind(xygrid[,2],xygrid[,1]))
# correct EI according to host distribution
hostgrid[is.na(hostgrid)] <- 0
EI <- EI*(hostgrid>0)
tab[,3] <- EI
GI <- GI*(hostgrid>0)
tab[,4] <- GI
exprNATOT <- exprNATOT|is.na(hostgrid)
}

EI[is.na(EI)] <- 0; tab[,3] <- EI
GI[is.na(GI)] <- 0; tab[,4] <- GI


############################################################
# calculations in metric system

tmpx <- sort(unique(tab[,2]))
tmpy <- sort(unique(tab[,1]))
dx <- min(tmpx[-1]-tmpx[-length(tmpx)]) 
dy <- min(tmpy[-1]-tmpy[-length(tmpy)]) 
dx <- round(dx,1)
dy <- round(dy,1)

if((min(tmpx)<0)&(max(tmpx)>0))
{
seqx1 <- seq(from=min(tmpx),to=max(tmpx[tmpx<0]),by=dx)
seqx2 <- seq(from=min(tmpx[tmpx>0]),to=max(tmpx),by=dx)
seqx <- rbind(t(t(seqx1)),t(t(seqx2)))
}
if(!((min(tmpx)<0)&(max(tmpx)>0)))
{
seqx <- seq(from=min(tmpx),to=max(tmpx),by=dx)
}
seqy <- seq(from=min(tmpy),to=max(tmpy),by=dy)

vectxy <- expand.grid(seqx,seqy)
vectxy[,1] <- round(vectxy[,1],3); vectxy[,2] <- round(vectxy[,2],3)   # Definition of the grid




# calculation of the real cell size dX in degrees:
dX <- NULL # vector of dx
dY <- dy
for (i in 1:length(seqx))
{
	if ((i==1)|(i==length(seqx)))
	{
	dX <- rbind(dX,dx)
	}
	if ((!(i==1))&(!(i==length(seqx))))
	{
	dX <- rbind(dX,((seqx[i+1]+seqx[i])/2)-((seqx[i]+seqx[i-1])/2))
	}
}

dXexpand <- rep(dx,length(vectxy[,1]))
for (i in 1:length(seqx))
{
	expr <- (vectxy[,1]==seqx[i])
	dXexpand[expr] <- dX[i]
		
}
dXexpand <- round(dXexpand,2)


# calculation of the area:
#--------------------------
longref <- 20
X <- (40000/360)*(vectxy[,1]-longref)*cos(vectxy[,2]*pi/180) 	# X-coordinates in metric system
Y <- (40000/360)*vectxy[,2]					# Y-coordinates in metric system
			
dxkm <- (40000/360)*dXexpand*cos(vectxy[,2]*pi/180)
dykm <- (40000/360)*dy

area <- abs(dxkm*dykm)


#**************************************
# Define projection WGS 194 and the resolution (given by CLIMEX)
projWGS84="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#Extent raster for the Climex model
extrast = raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs=projWGS84)
res(extrast)[1]=dx
res(extrast)[2]=dy

#****************************************


# transformation of EI and GI in vector associated with vectxy (expandgrid)
vectEI <- rep(0,length(vectxy[,1]))
vectGI <- rep(0,length(vectxy[,1]))
vectsea <- rep(T,length(vectxy[,1]))
for(i in 1:length(vectxy[,1]))
{
	xi <- vectxy[i,1]
	yi <- vectxy[i,2]
	expr <- ( (tab[,1]==yi) & (tab[,2]==xi) )
	if(sum(expr)==1) {vectEI[i] <- EI[expr]; vectGI[i] <- GI[expr];vectsea[i] <- F}
}

vectsea <- as.logical(vectsea)

vectexpr0 <- (vectEI==0)


plotRA <- function()
{
	par(bg = "lightgrey")
	plot(tab[,2],tab[,1],pch=19,cex=0.7,col="white",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Area of potential establishment",asp=1.3)
	points(tab[EI>0,2],tab[EI>0,1],pch=19,cex=0.7,col="orange")
	plot(shp,pch=".",add=T)
	par(bg = "transparent")
	print(" LEGEND:                                            ")
	print(" orange dots: EI>0, white dots: EI=0, grey: no data ")
}

plotGI <- function()
{
	par(bg = "lightgrey")
	plot(tab[,2],tab[,1],pch=19,cex=0.7,col="grey",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Growth potential",asp=1.3,bg="lightgrey")
	expr00 <- (GI==0)
	expr0 <- (GI>0)&(GI<=10)
	expr1 <- (GI>10)&(GI<=20)
	expr2 <- (GI>20)&(GI<=30)
	expr3 <- (GI>30)&(GI<=40)
	expr4 <- (GI>40)&(GI<=50)
	expr5 <- (GI>50)&(GI<=60)
	expr6 <- (GI>60)&(GI<=70)
	expr7 <- (GI>70)&(GI<=80)
	expr8 <- (GI>80)&(GI<=90)
	expr9 <- (GI>90)&(GI<=100)
	coulpal <- colorRampPalette(c("white","yellow","orange","red"),space="rgb")
	coul <- coulpal(11)
	points(tab[expr00,2],tab[expr00,1],col="white",pch=19)
	points(tab[expr0,2],tab[expr0,1],col=coul[2],pch=19)
	points(tab[expr1,2],tab[expr1,1],col=coul[3],pch=19)
	points(tab[expr2,2],tab[expr2,1],col=coul[4],pch=19)
	points(tab[expr3,2],tab[expr3,1],col=coul[5],pch=19)
	points(tab[expr4,2],tab[expr4,1],col=coul[6],pch=19)
	points(tab[expr5,2],tab[expr5,1],col=coul[7],pch=19)
	points(tab[expr6,2],tab[expr6,1],col=coul[8],pch=19)
	points(tab[expr7,2],tab[expr7,1],col=coul[9],pch=19)
	points(tab[expr8,2],tab[expr8,1],col=coul[10],pch=19)
	points(tab[expr9,2],tab[expr9,1],col=coul[11],pch=19)
	plot(shp,pch=".",add=T)
	par(bg = "transparent")
	print(" LEGEND:                                                                                          ")
	print(" graduated colors from white (GI=0) to red (GI=100) for the rescaled GI. Grey means no data")
}

plothabitat<- function()
{
	if(length(unique(hostgrid))<3)
	{
	par(bg = "lightgrey")
	plot(tab[,2],tab[,1],pch=19,cex=0.7,col="white",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Habitat distribution",asp=1.3)
	points(tab[hostgrid>0,2],tab[hostgrid>0,1],pch=19,cex=0.7,col="green")
	plot(shp,pch=".",add=T)
	print(" LEGEND:                                                          ")
	print(" green dots: habitat present, white dots: habitat absent, grey: no data ")
	par(bg = "transparent")
	}
	if(length(unique(hostgrid))>2)
	{
	par(bg = "lightgrey")
	seqh <- seq(from=min(hostgrid,na.rm=T),to=max(hostgrid,na.rm=T),length.out=10)
	expr0 <- (hostgrid==seqh[1])
	expr1 <- (hostgrid>seqh[1])&(hostgrid<=seqh[2])
	expr2 <- (hostgrid>seqh[2])&(hostgrid<=seqh[3])
	expr3 <- (hostgrid>seqh[3])&(hostgrid<=seqh[4])
	expr4 <- (hostgrid>seqh[4])&(hostgrid<=seqh[5])
	expr5 <- (hostgrid>seqh[5])&(hostgrid<=seqh[6])
	expr6 <- (hostgrid>seqh[6])&(hostgrid<=seqh[7])
	expr7 <- (hostgrid>seqh[7])&(hostgrid<=seqh[8])
	expr8 <- (hostgrid>seqh[8])&(hostgrid<=seqh[9])
	expr9 <- (hostgrid>seqh[9])&(hostgrid<=seqh[10])
	coulpal <- colorRampPalette(c("white","lightgreen","darkgreen"),space="rgb")
	coul <- coulpal(11)
	plot(tab[,1],tab[,2],pch=".",col="white",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Habitat distribution",asp=1.3)
	if(seqh[1]==0){points(tab[expr0,2],tab[expr0,1],col="white",pch=19)}
	if(seqh[1]>0){points(tab[expr0,2],tab[expr0,1],col=coul[2],pch=19)}
	points(tab[expr1,2],tab[expr1,1],col=coul[3],pch=19)
	points(tab[expr2,2],tab[expr2,1],col=coul[4],pch=19)
	points(tab[expr3,2],tab[expr3,1],col=coul[5],pch=19)
	points(tab[expr4,2],tab[expr4,1],col=coul[6],pch=19)
	points(tab[expr5,2],tab[expr5,1],col=coul[7],pch=19)
	points(tab[expr6,2],tab[expr6,1],col=coul[8],pch=19)
	points(tab[expr7,2],tab[expr7,1],col=coul[9],pch=19)
	points(tab[expr8,2],tab[expr8,1],col=coul[10],pch=19)
	points(tab[expr9,2],tab[expr9,1],col=coul[11],pch=19)
	plot(shp,pch=".",add=T)
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

resol <- paste("The grid resolution is: ",dx," decimal degrees")
print(resol)

resoldxkm <- paste("The mean cell size: ",round(mean(dxkm))," km longitude and ",round(mean(dykm))," km latitude")
print(resoldxkm)

calarea <- paste("The mean area of a cell is: ",round(mean(area))," km²",sep="")
print(calarea)

calnRA <- paste("The number of cells in the area of potential establishment is: ",sum(EI>0),sep="")
print(calnRA)

calntot1 <- paste("The total number of cells in SLG, LG-Econ, radial, rand-radial is: ",length(EI),sep="")
print(calntot1)

calntot2 <- paste("The total number of cells in the grid used by the dispersal kernel model is: ",length(vectEI),sep="")
print(calntot2)

print("                                                          ")
print("-----  Type the commands below to call the functions -----")
print("                                                          ")
}

printinfo()

##############################################
# SIMPLE LOGISTIC GROWTH: Model C
##############################################


slg <- function(N0=1,lmax=7,movie=F,t=1,export=F,name="SLG-DD-output")
{
	growth <- GI*log(lmax)/100
	val <- 100*N0*exp(growth*t)/(100+N0*(exp(growth*t)-1)) 
	expr0 <- (EI==0)
	res<- (!expr0)*val # percentage of the cell invaded

	if(!movie)
	{
	nom <- paste("SLG model t = ",t,sep="")
	plot(xygrid[,2],xygrid[,1],pch=".",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main=nom,asp=1.3,col="white")
	points(points(xygrid[expr0,2],xygrid[expr0,1],pch=19,cex=0.7,col="grey"))
	expr25 <- ((res>0) & (res<=25))
	expr50 <- ((res>=25) & (res<50))
	expr75 <- ((res>=50) & (res<75))
	expr100 <- (res>=75)
	points(xygrid[expr25,2],xygrid[expr25,1],pch=19,cex=0.7,col="blue")
	points(xygrid[expr50,2],xygrid[expr50,1],pch=19,cex=0.7,col="green")
	points(xygrid[expr75,2],xygrid[expr75,1],pch=19,cex=0.7,col="orange")
	points(xygrid[expr100,2],xygrid[expr100,1],pch=19,cex=0.7,col="red")
	plot(shp,pch=".",add=T)
	}

	if(movie)
	{
	for (k in 1:t)
	{
	if(k<t){nom <- paste("SLG model t = ",k," (CLICK ON THE MAP)",sep="")}
	if(k==t){nom <- paste("SLG model t = ",k,sep="")}
	plot(xygrid[,2],xygrid[,1],pch=".",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main=nom,asp=1.3,col="white")
	val <- 100*N0*exp(growth*k)/(100+N0*(exp(growth*k)-1)) 
	res<- (!expr0)*val # percentage of the cell invaded
	points(points(xygrid[expr0,2],xygrid[expr0,1],pch=19,cex=0.7,col="grey"))
	expr25 <- ((res>0) & (res<=25))
	expr50 <- ((res>=25) & (res<50))
	expr75 <- ((res>=50) & (res<75))
	expr100 <- (res>=75)
	points(xygrid[expr25,2],xygrid[expr25,1],pch=19,cex=0.7,col="blue")
	points(xygrid[expr50,2],xygrid[expr50,1],pch=19,cex=0.7,col="green")
	points(xygrid[expr75,2],xygrid[expr75,1],pch=19,cex=0.7,col="orange")
	points(xygrid[expr100,2],xygrid[expr100,1],pch=19,cex=0.7,col="red")
	plot(shp,pch=".",add=T)
	if(k<t){locator(1)}
	}

	}


	if(export)
	{
	# to export in the LAEA projection (for combination in MCAS)
	res[exprNATOT] <- NA
	xyz = data.frame(x=xygrid[,2], y=xygrid[,1], z=res)
	coordinates(xyz) = ~x+y
	rasterSLG=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterSLG,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}


print("LEGEND:                                               ")
print("Nt: population abundance (%)                          ")
print("blue dots:  0 < Nt < = 25, green dots: 25 < Nt < = 50, orange dots: 50 < Nt < = 75, red dots: 75 < Nt < = 100, grey dots: Nt=0, white: no data")
return(list(slg=res,ntot=dim(tab)[1],nRA=sum(!expr0),sum0=sum(res==0),sum25=sum(expr25),sum50=sum(expr50),sum75=sum(expr75),sum100=sum(expr100)))
}


####################################################
# LOGISTIC GROWTH BASED ON ECONOMIC VALUES: model A
####################################################

lgecon <- function(N0=1,r=1,t=1,econraster=F,mult=1,sim=1,habitatraster=F,valperhost=NULL,econformat=NULL,export=F,name="LGECON-DD-Output")
{

if(econraster & habitatraster){print("You should choose between econraster and habitatraster")}

if (econraster)
{
	# read economic file and get values
	econ <- raster(paste("econ",econformat,sep=""))
	econvalues <- xyValues(econ,cbind(xygrid[,2],xygrid[,1]))
	econmem <- econvalues
}
if (habitatraster & (!econraster))
{
	if(is.null(valperhost)){print("ERROR: valperhost is required in this case")}
	econvalues <- hostgrid*valperhost
	econmem <- econvalues
}
if (econraster|habitatraster)
{
exprEconNA <- is.na(econvalues)
if(sum(exprEconNA)>0){print("CAUTION: the model excludes areas where economic values are missing")}
econvalues[is.na(econvalues)] <- 0 # where there is missing data, we assume that the value is 0 
econvalues[econvalues<0] <- 0  
econvalues <- econvalues*mult	   # to convert in the right unit
# application to the area of the cell
econvaluereal <- econvalues
for (i in (1:length(econvalues))[econvalues>0])
{
expr <- ((tab[i,2]-vectxy[,1])^2+(tab[i,1]-vectxy[,2])^2) == min((tab[i,2]-vectxy[,1])^2+(tab[i,1]-vectxy[,2])^2)
econvaluereal[i] <- econvalues[i]*area[expr] 
} 
breakvalues <- seq(0,max(econvalues),length.out=6)
}


if(!(econraster|habitatraster)) 
{
econvalues <- runif(n=length(EI),min=1,max=100)
exprEconNA <- is.na(econvalues)
econmem <- econvalues
# application to the area of the cell
econvaluereal <- econvalues
for (i in (1:length(econvalues))[econvalues>0])
{
expr <- ((tab[i,2]-vectxy[,1])^2+(tab[i,1]-vectxy[,2])^2) == min((tab[i,2]-vectxy[,1])^2+(tab[i,1]-vectxy[,2])^2)
econvaluereal[i] <- econvalues[i]*area[expr]
} 
breakvalues <- c(0,20,40,60,80,100)
print("CAUTION: no economic file was provided => economic values are taken at random")
}
	
	nRA <- sum(EI>0) # number of pixels included in the risk area (where EI>0)
	val <- 100*N0*exp(r*t)/(100+N0*(exp(r*t)-1)) # is the percentage of assets invaded
	ninv <- round(val*nRA/100) # number of pixels invaded 

	vectbest <- (EI<0); vectbest[EI==0] <- (-1); vectbest[exprNATOT|exprEconNA] <- (-9999) 	# 0 -non invaded- everywhere at the beginning,  -1 for non suitable (EI=0) and then -9999 for missing data
	vectworst <- vectbest
	vectrand <- vectbest 

	# the code is going to put 1 where the cell is invaded. 
	nID <- 1:nxy    	# ID of the line 
	nIDselect <- nID[(EI>0)&(!exprEconNA)]


if(ninv>0)
{
par(mfrow=c(3,2))
# selection of the points with the risk area and where we have economic data:
tabRA <- tab[(EI>0)&(!exprEconNA),]
ecovalTOT <- econvalues[(EI>0)&(!exprEconNA)]
ecovalTOTreal <- econvaluereal[(EI>0)&(!exprEconNA)]


# Best case scenario: pixels with low values are invaded first
#--------------------------------------------------------------
orderbest <- order(ecovalTOT,decreasing=FALSE)
tabbest <- tabRA[orderbest,]
ecobest <- ecovalTOT[t(orderbest)]
ecobestreal <- ecovalTOTreal[t(orderbest)]
resbest <- sum(ecobestreal[1:ninv])
nIDbest <- nIDselect[t(orderbest)]; nIDbest <- nIDbest[1:ninv]		# ID of the lines of invaded cells
vectbest[nIDbest] <- 1							# update the file for invaded cells
	plot(xygrid[,2],xygrid[,1],pch=19,cex=0.3,col="grey",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Best case",asp=1.3,cex.axis=1,cex.lab=1)
	points(tab[EI>0,2],tab[EI>0,1],pch=19,cex=0.3,col="orange")
	points(tabRA[,2],tabRA[,1],pch=19,cex=0.3,col="orange")
	points(tabbest[1:ninv,2],tabbest[1:ninv,1],pch=19,cex=0.5,col="red")
	plot(shp,pch=".",add=T)

	hist(ecobest[1:ninv],main="Best case",xlab="Economic value (euros/km²)",breaks=breakvalues)


# Medium case scenario: random selection of the pixels
#-----------------------------------------------------

# we make 100 replicate simulations:
memreal <- NULL
mem <- NULL
for (i in 1:sim)
{
ni <- sample(1:length(tabRA[,1]),size=ninv,replace=F)
resmed <- sum(ecovalTOTreal[ni])
nIDrand <- nIDselect[ni]
vectrand[nIDrand] <- 1	
memreal <- rbind(memreal,cbind(i,resmed))
mem <- rbind(mem,rbind(t(t(ecovalTOT[ni]))))
}
	plot(xygrid[,2],xygrid[,1],pch=19,cex=0.3,col="grey",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Random case",asp=1.3,cex.axis=1,cex.lab=1)
	points(tab[EI>0,2],tab[EI>0,1],pch=19,cex=0.3,col="orange")
	points(tabRA[,2],tabRA[,1],pch=19,cex=0.3,col="orange")
	points(tabRA[ni,2],tabRA[ni,1],pch=19,cex=0.5,col="red")
	plot(shp,pch=".",add=T)

	nom <- paste("Random case - ",sim," replicate simulations",sep="")
	hist(mem,main=nom,xlab="Economic value (euros/km²)",breaks=breakvalues)
	resrand <- summary(memreal[,2])

# Worst case scenario: pixels with high values are invaded first
#-----------------------------------------------------------------
orderworst <- order(ecovalTOT,decreasing=TRUE)
tabworst <- tabRA[orderworst,]
ecoworst <- ecovalTOT[t(orderworst)]
ecoworstreal <- ecovalTOTreal[t(orderworst)]
resworst <- sum(ecoworstreal[1:ninv])
nIDworst <- nIDselect[t(orderworst)]; nIDworst <- nIDworst[1:ninv]		# ID of the lines of invaded cells
vectworst[nIDworst] <- 1							# update the file for invaded cells
	plot(xygrid[,2],xygrid[,1],pch=19,cex=0.3,col="grey",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Worst case",asp=1.3,cex.axis=1,cex.lab=1)
	points(tab[EI>0,2],tab[EI>0,1],pch=19,cex=0.3,col="orange")
	points(tabRA[,2],tabRA[,1],pch=19,cex=0.3,col="orange")
	points(tabworst[1:ninv,2],tabworst[1:ninv,1],pch=19,cex=0.5,col="red")
	plot(shp,pch=".",add=T)

	hist(ecoworst[1:ninv],main="Worst case",xlab="Economic value (euros/km²)",breaks=breakvalues)

par(mfrow=c(1,1))
}
if(ninv==0){ 
	plot(xygrid[,2],xygrid[,1],pch=19,cex=0.7,col="grey",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="LG-Econ",asp=1.3,cex.axis=1,cex.lab=1)
	points(tab[EI>0,2],tab[EI>0,1],pch=19,cex=0.7,col="orange")
	plot(shp,pch=".",add=T)
	resworst=0
	resbest=0
	resrand=0
	}

	if(export)
	{
	vectbest[exprNATOT|exprEconNA] <- NA
	xyz = data.frame(x=xygrid[,2], y=xygrid[,1], z=vectbest)
	coordinates(xyz) = ~x+y
	rasterbest=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterbest,paste(name,"_best.tiff",sep=""),format="GTiff")

	vectworst[exprNATOT|exprEconNA] <- NA
	xyz = data.frame(x=xygrid[,2], y=xygrid[,1], z=vectworst)
	coordinates(xyz) = ~x+y
	rasterworst=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterworst,paste(name,"_worst.tiff",sep=""),format="GTiff")

	vectrand[exprNATOT|exprEconNA] <- NA
	xyz = data.frame(x=xygrid[,2], y=xygrid[,1], z=vectrand)
	coordinates(xyz) = ~x+y
	rasterrand=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterrand,paste(name,"_rand.tiff",sep=""),format="GTiff")

	tmp <- paste("NOTE that three output files called ", name," were exported",sep="")
	print(tmp)
	}

print("LEGEND:                                          ")
print("red dots: invaded cells, orange dots: non invaded but suitable cells, grey dots: EI=0, white: no data on climate suitability")
return(list(econvalues=econmem,nRA=nRA,ninv=ninv,ntot=dim(tab)[1],pniche=val,worst=resworst,best=resbest,rand=resrand))
}


##############################################
# RADIAL RANGE EXPANSION: model B
##############################################


radial <- function(RR=20,t=1,coord=NULL,export=F,name="Radial-DD-output")
{
	expr0 <- (EI==0)

# ENTRY POINT:
#-------------
# if random selection of entry point:
if(is.null(coord))
{
entry <- sample(1:sum(!expr0),size=1)
tmp <- tab[!expr0,]
X <- tmp[entry,2]; Xtot <- X
Y <- tmp[entry,1]; Ytot <- Y
exprentry <- (tab[,1]==Y)&(tab[,2]==X)
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
y <- coord[2*i]	        # latitude,  not necessarily from the grid cell=> will find the closest grid cell
tmp <- tab[!expr0,]
expr <- (min((tmp[,2]-x)^2+(tmp[,1]-y)^2)==(tmp[,2]-x)^2+(tmp[,1]-y)^2)
X <- tmp[expr,2]
Y <- tmp[expr,1]
if(sum(expr)>1){X <- X[1]; Y <- Y[1]}
exprentry <- (tab[,1]==Y)&(tab[,2]==X)
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
	tmp <- round(sin((pi/180)*xygrid[i,1])*sin((pi/180)*xygrid[nnon,1])+cos((pi/180)*xygrid[i,1])*cos((pi/180)*xygrid[nnon,1])*cos((pi/180)*(xygrid[nnon,2]-xygrid[i,2])),7)
	ds <- acos(tmp)
	R <- 20000*ds/pi
	expr <- (R <= RR*t)
	ninvaded <- nnon[expr]
	vect2[ninvaded] <- 1
	}
	resinvEI <- xygrid[vect2==1,]
	resnonEI <- xygrid[vect2==0,]

	if(export)
	{
	vect2[exprNATOT] <- NA
	xyz = data.frame(x=xygrid[,2], y=xygrid[,1], z=vect2)
	coordinates(xyz) = ~x+y
	rasterrad=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterrad,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}


longref <- 20
Xkm <- (40000/360)*(tab[,2]-longref)*cos(tab[,1]*pi/180) # X-coordinates in metric system
Ykm <- (40000/360)*tab[,1]				 # Y-coordinates in metric system
Xkminv <- (40000/360)*(resinvEI[,2]-longref)*cos(resinvEI[,1]*pi/180) 
Ykminv <- (40000/360)*resinvEI[,1]
Xtotkm <- (40000/360)*(Xtot-longref)*cos(Ytot*pi/180) 	 # X-coordinates in metric system
Ytotkm <- (40000/360)*Ytot				 # Y-coordinates in metric system

plot(tab[,2],tab[,1],col="white",pch=".",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Radial expansion",asp=1.3)
points(tab[expr0,2],tab[expr0,1],col="grey",pch=19,cex=0.7)
points(tab[!expr0,2],tab[!expr0,1],col="orange",pch=19,cex=0.7)
points(resinvEI[,2],resinvEI[,1],col="red",pch=19,cex=0.7)
points(Xtot,Ytot,pch=19)
plot(shp,pch=".",add=T)

pniche=100*dim(resinvEI)[1]/sum(!expr0)
print("LEGEND:                                                   ")
print("red dots: invaded cells, orange dots: non invaded but suitable cells, grey dots: EI=0, white: no data")
return(list(radial=resinvEI,ninv=dim(resinvEI)[1],nRA=sum(!expr0),ntot=dim(tab)[1],pniche=pniche))
}



##############################################
# COMBINED MODEL: model A random + model B
##############################################

radialrand <- function(N0=1,r=1,t=1,RR=20,coord=NULL,export=F,name="RadialRand-DD-output")
{
	expr0 <- (EI==0)

# ENTRY POINT:
#-------------
# if random selection of entry point:
if(is.null(coord))
{
entry <- sample(1:sum(!expr0),size=1)
tmp <- tab[!expr0,]
X <- tmp[entry,2]
Y <- tmp[entry,1]
exprentry <- (tab[,1]==Y)&(tab[,2]==X)
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
tmp <- tab[!expr0,]
expr <- (min((tmp[,2]-x)^2+(tmp[,1]-y)^2)==(tmp[,2]-x)^2+(tmp[,1]-y)^2)
X <- tmp[expr,2]
Y <- tmp[expr,1]
if(sum(expr)>1){X <- X[1]; Y <- Y[1]}
exprentry <- (tab[,1]==Y)&(tab[,2]==X)
tmpientry <- (1:nxy)[exprentry] 
ientry <- rbind(ientry,tmpientry)
Xtot <- rbind(Xtot,X)
Ytot <- rbind(Ytot,Y)
}
}

# calculation for time t
#------------------------
	vect <- (EI<0); vect[EI==0] <- (-1); vect[exprNATOT] <- (-9999) 	# 0 -non invaded- everywhere at the beginning,  -1 for non suitable (EI=0) and then -9999 for missing data
	vect[ientry] <- 1	# 1 where the cell is invaded. 
	nID <- 1:nxy    	# ID of the line 
	ninv <- nID[vect==1] 
	nnon <- nID[vect==0]	
	vect2 <- vect

resinv <- tab[ientry,]
resnon <- tab[-ientry,]

	mem <- NULL
	for (i in ninv)
	{
	tmp <- round(sin((pi/180)*xygrid[i,1])*sin((pi/180)*xygrid[nnon,1])+cos((pi/180)*xygrid[i,1])*cos((pi/180)*xygrid[nnon,1])*cos((pi/180)*(xygrid[nnon,2]-xygrid[i,2])),7)
	ds <- acos(tmp)
	R <- 20000*ds/pi
	expr <- (R <= RR*t)
	ninvaded <- nnon[expr]
	vect2[ninvaded] <- 1
	}

	resinvEI <- xygrid[vect2==1,]
	resnonEI <- xygrid[vect2==0,]


# Random selection of the pixels
#-----------------------------------------------------
	nRA <- sum(EI>0) # number of pixels included in the risk area (where EI>0)
	val <- 100*N0*exp(r*t)/(100+N0*(exp(r*t)-1)) # is the percentage of assets invaded
	ninv <- round(val*nRA/100) # number of pixels invaded 

	if (ninv>sum(vect2==1)){vect3 <- vect2; tabINV <- xygrid[vect2==1,]}
	if (ninv<=sum(vect2==1))
	{
	ni <- sample(1:sum(vect2==1),size=ninv,replace=F)
	nvect2 <- nID[vect2==1]
	vect3 <- vect
	vect3[nvect2[ni]] <- 1
	tabINV <- xygrid[nvect2[ni],]
	}


	if(export)
	{
	vect3[exprNATOT] <- NA
	xyz = data.frame(x=xygrid[,2], y=xygrid[,1], z=vect3)
	coordinates(xyz) = ~x+y
	rasterradrand=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterradrand,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
	}


longref <- 20
Xkm <- (40000/360)*(tab[,2]-longref)*cos(tab[,1]*pi/180) # X-coordinates in metric system
Ykm <- (40000/360)*tab[,1]				 # Y-coordinates in metric system
Xkminv <- (40000/360)*(tabINV[,2]-longref)*cos(tabINV[,1]*pi/180) 
Ykminv <- (40000/360)*tabINV[,1]
Xtotkm <- (40000/360)*(Xtot-longref)*cos(Ytot*pi/180) 	 # X-coordinates in metric system
Ytotkm <- (40000/360)*Ytot				 # Y-coordinates in metric system

plot(tab[,2],tab[,1],col="white",pch=".",xlab="Longitude",ylab="Latitude",cex.lab=1.5,cex.axis=1.5,main="Random radial expansion",asp=1.3)
points(tab[expr0,2],tab[expr0,1],col="grey",pch=19,cex=0.7)
points(tab[!expr0,2],tab[!expr0,1],col="orange",pch=19,cex=0.7)
points(tabINV[,2],tabINV[,1],col="red",pch=19,cex=0.7)
points(Xtot,Ytot,pch=19)
plot(shp,pch=".",add=T)

pniche=100*dim(tabINV)[1]/sum(!expr0)
print("LEGEND:                                                   ")
print("red dots: invaded cells, orange dots: non invaded but suitable cells, grey dots: EI=0, white: no data")
return(list(radialrand=tabINV,ninv=dim(tabINV)[1],nRA=sum(!expr0),ntot=dim(tab)[1],pniche=pniche))
}



####################################################
# DISPERSAL KERNEL: DETERMINISTIC VERSION - model D
####################################################

plotkernel <- function(vectpop0,presence=NULL,t=NULL,threshold=NULL,legend=T)
{
if(is.null(threshold))
{
# plot the figure:
expr0 <- vectpop0<(10^(-6)) & (!vectsea)
expr1 <- (vectpop0>=(10^(-6)) ) & (vectpop0<(10^(-5)))
expr2 <- (vectpop0>=(10^(-5)) ) & (vectpop0<(10^(-4)))
expr3 <- (vectpop0>=(10^(-4)) ) & (vectpop0<(10^(-3)))
expr4 <- (vectpop0>=(10^(-3)) ) & (vectpop0<(10^(-2)))
expr5 <- (vectpop0>=(10^(-2)) ) & (vectpop0<(10^(-1)))
expr6 <- (vectpop0>=(10^(-1)) ) & (vectpop0<1)
expr7 <- (vectpop0>=1) & (vectpop0<10)
expr8 <- (vectpop0>=10)
coulpal <- colorRampPalette(c("white","yellow","orange","red"),space="rgb")
coul <- coulpal(9)

if(!is.null(t)){nom <- paste("Dispersal kernel for time t = ",t,sep="")}
if(is.null(t)){nom <- paste("Dispersal kernel")}


plot(vectxy[,1],vectxy[,2],pch=19,asp=1.3,xlab="Longitude", ylab="Latitude",col="white",main=nom,cex.lab=1.5,cex.axis=1.5)
points(vectxy[vectsea,1],vectxy[vectsea,2],col="grey",pch=19)
points(vectxy[expr0,1],vectxy[expr0,2],col=coul[1],pch=19,cex=0.5)
points(vectxy[expr1,1],vectxy[expr1,2],col=coul[2],pch=19)
points(vectxy[expr2,1],vectxy[expr2,2],col=coul[3],pch=19)
points(vectxy[expr3,1],vectxy[expr3,2],col=coul[4],pch=19)
points(vectxy[expr4,1],vectxy[expr4,2],col=coul[5],pch=19)
points(vectxy[expr5,1],vectxy[expr5,2],col=coul[6],pch=19)
points(vectxy[expr6,1],vectxy[expr6,2],col=coul[7],pch=19)
points(vectxy[expr7,1],vectxy[expr7,2],col=coul[8],pch=19)
points(vectxy[expr8,1],vectxy[expr8,2],col=coul[9],pch=19)
plot(shp,pch=".",add=T)
points(vectxy[presence>0,1],vectxy[presence>0,2],col="black",pch=19)

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

plot(vectxy[,1],vectxy[,2],pch=19,asp=1.3,xlab="Longitude", ylab="Latitude",col="white",main=nom,cex.lab=1.5,cex.axis=1.5)
points(vectxy[vectsea,1],vectxy[vectsea,2],col="grey",pch=19)
points(vectxy[vectpop0>threshold,1],vectxy[vectpop0>threshold,2],col="red",pch=19,cex=0.7)
if(!is.null(presence)){points(vectxy[presence>0,1],vectxy[presence>0,2],col="black",pch=19)}
plot(shp,pch=".",add=T)

if(legend)
{
print("LEGEND:                                                                ")
print("Nt: population abundance (%)                                           ")
print("red dots: cells where Nt > = threshold, white: cells where Nt < threshold. Grey means no data")
}
}

}


exportkernel <- function(vectpop0,name="kernelOutput")
{
	vectpop0[exprNATOT] <- NA 
	xyz = data.frame(x=vectxy[,1], y=vectxy[,2], z=vectpop0)
	coordinates(xyz) = ~x+y
	rasterkernel=pointsToRaster(extrast,xyz,values=xyz$z,fun=mean)
	writeRaster(rasterkernel,paste(name,".tiff",sep=""),format="GTiff")
	tmp <- paste("NOTE that an output file called ", name," was exported",sep="")
	print(tmp)
}



dispk <- function(N0=NULL, t=1,lmax=7, P=1, p=50, u=100,presencefile=T,nentry=NULL)
{

#--------------------------------------------
if(!presencefile)
{# case 3
entry <- sample(1:sum(!vectexpr0),size=nentry)
tmp <- vectxy[!vectexpr0,]
xentry <- tmp[entry,1]
yentry <- tmp[entry,2]
presence <- cbind(xentry,yentry)
presencegrid <- rep(F,length(vectxy[,1]))
for (i in 1:nentry)
{
	expr <- (min( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ) == ( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ))
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
presencegrid <- rep(F,length(vectxy[,1]))
for (i in 1:dim(presence)[1])
{
	expr <- (min( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ) == ( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ))
	presencegrid[expr] <- T
}
presence <- presencegrid
# vector of TRUE/FALSE
}

if(ncolumn==3)
{ # case 2
presencegrid <- rep(0,length(vectxy[,1]))
count  <- rep(0,length(vectxy[,1]))
for (i in 1:dim(presence)[1])
{
	expr <- (min( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ) == ( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ))
	presencegrid[expr] <- presencegrid[expr]+presence[i,3]
	count[expr] <- count[expr]+1
}
presence <- presencegrid
presence[!presence==0] <- presence[!presence==0]/count[!presence==0] 
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
popdisp <- rep(0,length(vectxy[,1]))


for (i in (1:length(vectxy[,1]))[(!vectexpr0)|(presence>0)]) # we assume that a source point is necessarily where the EI>0 (to reduce the potential source points and calculations)
{ 
	x0 <- vectxy[i,1]
	y0 <- vectxy[i,2]
	tmp <- round(sin((pi/180)*y0)*sin((pi/180)*vectxy[,2])+cos((pi/180)*y0)*cos((pi/180)*vectxy[,2])*cos((pi/180)*(vectxy[,1]-x0)),7)
	ds <- acos(tmp)
	R <- 20000*ds/pi
	s <- sin(R/(20000/pi))/(R/(20000/pi))
	exprNA <- is.na(s)
	s[exprNA] <- 1 # for a distance=0, we have s=NA but it should be 1 (sin(x)~x)
	f <- (1/s)*(1/((u^2)*p*pi)) * gamma((p+1)/2)/( gamma((p-1)/2)*(1+(R^2)/(p*(u^2)))^((p+1)/2) )
	probaf <- f*area

	popdisp <- popdisp + probaf*P*vectpop0[i] # we consider only individuals engaged in dispersal
	popdisp[vectEI==0] <- 0
}


popdisp <- popdisp + (1-P)*vectpop0 # we add individuals staying at the same place

popdisp[vectexpr0] <- 0 # we assume that individuals die where the EI=0

growth <- vectGI*log(lmax)/100
vectpop0 <- (100*popdisp*exp(growth)/(100+popdisp*(exp(growth)-1)))

plotkernel(vectpop0=vectpop0,presence=presence,t=k,legend=F)
}

print("LEGEND:                                                                ")
print("Nt: population abundance (%)                                           ")
print(" graduated colors from white (Nt < 10^(-6) %) to yellow, orange and red (Nt > = 10%). Grey means no data")

# percentage of the niche invaded:
pniche <- (100*sum(vectpop0>1)/sum(!vectexpr0))

return(list(dispk=vectpop0,pniche=pniche,ntot=length(vectpop0),nRA=sum(!vectexpr0),ninv= sum(vectpop0>1),presence=presence))
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


#################################################
# DISPERSAL KERNEL: STOCHASTIC VERSION - model D
#################################################


sto <- function(N0=1, t=1, P=1, p=50, u=100,nentry=1,presencefile=T,lmax=7)
{

#------------------------------
if(!presencefile)
{# case 3
entry <- sample(1:sum(!vectexpr0),size=nentry)
tmp <- vectxy[!vectexpr0,]
xentry <- tmp[entry,1]
yentry <- tmp[entry,2]
presence <- cbind(xentry,yentry)
presencegrid <- rep(F,length(vectxy[,1]))
for (i in 1:nentry)
{
	expr <- (min( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ) == ( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ))
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
presencegrid <- rep(F,length(vectxy[,1]))
for (i in 1:dim(presence)[1])
{
	expr <- (min( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ) == ( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ))
	presencegrid[expr] <- T
}
presence <- presencegrid
# vector of TRUE/FALSE
}

if(ncolumn==3)
{ # case 2
presencegrid <- rep(0,length(vectxy[,1]))
count  <- rep(0,length(vectxy[,1]))
for (i in 1:dim(presence)[1])
{
	expr <- (min( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ) == ( ((presence[i,1]-vectxy[,1])^2)+((presence[i,2]-vectxy[,2])^2) ))
	presencegrid[expr] <- presencegrid[expr]+presence[i,3]
	count[expr] <- count[expr]+1
}
presence <- presencegrid
presence[!presence==0] <- presence[!presence==0]/count[!presence==0] # we calculate the mean population density in case there is more than one presence point close to the grid point 
# vector of population density
}

}
#------------------------------

if(is.logical(presence)){presence <- as.numeric(presence)*N0}
vectpop0 <- presence

for (k in 1:t)
{
print(c("t = ",k))
popdisp <- rep(0,length(vectxy[,1]))

for (i in (1:length(vectxy[,1]))[(!vectexpr0)|(presence>1)]) # we assume that a source point is necessarily where the EI>0 (to reduce the potential source points and calculations)
{ 
x0 <- vectxy[i,1]
y0 <- vectxy[i,2]

	tmp <- round(sin((pi/180)*y0)*sin((pi/180)*vectxy[,2])+cos((pi/180)*y0)*cos((pi/180)*vectxy[,2])*cos((pi/180)*(vectxy[,1]-x0)),7)
	ds <- acos(tmp)
	R <- 20000*ds/pi
	s <- sin(R/(20000/pi))/(R/(20000/pi))
	exprNA <- is.na(s)
	s[exprNA] <- 1 
f <- (1/s)*(1/((u^2)*p*pi)) * gamma((p+1)/2)/( gamma((p-1)/2)*(1+(R^2)/(p*(u^2)))^((p+1)/2) )
probaf <- f*area

#-------------------

cumulsup <- cumsum(probaf)
cumulinf <- cumulsup[-length(cumulsup)]; cumulinf <- rbind(0,t(t(cumulinf)))

dispi <- round(P*vectpop0[i]) # number of "units" that will disperse
if (dispi>0)
{
cellp <- runif(n=dispi,min=0,max=max(cumulsup))
for (ii in 1:dispi)
{
exprcell <- ( (cellp[ii] >= cumulinf) & (cellp[ii] < cumulsup) )
popdisp[exprcell] <- popdisp[exprcell] + 1 
}
}


#-------------------
}
popdisp <- popdisp + (1-P)*vectpop0 # we add individuals staying at the same place

popdisp[vectexpr0] <- 0 # we assume that individuals die where the EI=0


growth <- vectGI*log(lmax)/100
vectpop0 <- (100*popdisp*exp(growth)/(100+popdisp*(exp(growth)-1)))

}

return(list(res=vectpop0,presence=presence))
}



dispksto <- function(N0=1, t=1, P=1, p=50, u=100,nentry=1,presencefile=T,lmax=7,sim=200)
{

ressto <- rep(0,length(vectxy[,1]))
for (n in 1:sim)
{
print(c("sim = ",n))
tmp <- sto(N0=N0, t=t, P=P, p=p, u=u,nentry=nentry,lmax=lmax)
ressto <- ressto+tmp$res
}
resfinal <- ressto/sim

presence <- tmp$presence

plotkernel(vectpop0=resfinal,presence=presence,t=t,legend=F)


# percentage of the niche invaded:
pniche <- (100*sum(resfinal>1)/sum(!vectexpr0))

return(list(dispksto=resfinal,pniche=pniche,ntot=length(resfinal),nRA=sum(!vectexpr0),ninv= sum(resfinal>1) ))
}




##########################################################################
# supplementary function(s)

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
tmp <- tab[EI>0,]
expr <- (min((tmp[,2]-x)^2+(tmp[,1]-y)^2)==(tmp[,2]-x)^2+(tmp[,1]-y)^2)
X <- tmp[expr,2]
Y <- tmp[expr,1]
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



