### Simulate spatially-explicit captures of stoats based on the distance to the centre of the home range
### Based on real data from South Ronaldsay - need updating the spatial extent of SR

## Packages
library(MASS)
library(MBESS)
library(adehabitatHR)
library(sp)
library(spatstat)
library(raster)
library(geosphere)
library(maptools)
library(rgdal)
library(SDMTools)
library(sf)
library(beepr)
library(rgeos)


set.seed(2017)

##### STOAT-SPECIFIC SIMULATIONS
### Load the shape of South Ronaldsay
sr<-readOGR("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SR_area.shp")

sr

plot(sr)

### Convert to raster for further manipulation later
#sr.raster<-raster(xmn=extent(sr)[1], xmx=extent(sr)[2], ymn=extent(sr)[3], ymx=extent(sr)[4], crs=crs(sr), resolution=5)

#sr.raster

#raster.sr<-rasterize(sr, sr.raster, field="OBJECTID", FUN=sum)

#values(raster.sr)[values(raster.sr)!=0]<-1

#values(raster.sr)[is.na(values(raster.sr))==TRUE]<-0

#### Load the trap data (trap positions)

trap.pos<-readOGR("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/TrapLocations.shp")

trap.pos

plot(trap.pos, add=TRUE)

### Extent
#ext.sr<-extent(sr)

#mask.sr<-matrix(as.logical(values(raster.sr)), ncol=ncol(raster.sr), nrow=nrow(raster.sr), byrow=TRUE)

#win.sr<-owin(xrange=c(ext.sr[1], ext.sr[2]), yrange=c(ext.sr[3], ext.sr[4]), mask=mask.sr)  ### Spatial window - for the point process

##### Total area of the island (ha)

tot.area<-sum(sr$Shape_Area)/10000

## Total number of stoats
init.density<-10/100			#### Density in individuals per ha

tot.stoat<-rpois(1, init.density*tot.area)

tot.stoat

### All potential home range centres
#x.stoat<-sample(win.sr$xcol, length(win.sr$xcol))

#y.stoat<-sample(win.sr$yrow, length(win.sr$yrow))

#tot.pos<-expand.grid(x.stoat=x.stoat, y.stoat=y.stoat)

#nrow(tot.pos)

#pot.stoat<-ppp(tot.pos$x.stoat, tot.pos$y.stoat, window=win.sr)

#pot.stoat$n

# All potential home ranges within polygon
HRcenters <- spsample(sr, n = tot.stoat, "random", iter = 20)
HRcentersDF <- as.data.frame(HRcenters)

#### Choose at random within the window defined
#samp.hr<-sample(c(1:pot.stoat$n), tot.stoat, replace=FALSE)                    ### Select home range centres at random

#centre.stoat<-ppp(pot.stoat$x[samp.hr], pot.stoat$y[samp.hr], window=win.sr)            ### Home range centres

##################################### TRAPS
n.trap<-length(trap.pos$ID)		### Number of traps

n.day<-20				### Number of days

win.trap<-owin(xrange=c(extent(trap.pos)[1], extent(trap.pos)[2]), yrange=c(extent(trap.pos)[3], extent(trap.pos)[4]))

sim.trap<-ppp(coordinates(trap.pos)[,1], coordinates(trap.pos)[, 2], win=win.trap, check=FALSE)  	#### Real traps

mean(nndist(sim.trap))  ### Check mean distance between traps (metres)

### Plot both the traps and the stoats

#points(centre.stoat$x, centre.stoat$y,col="black", pch=19, cex=2)
points(HRcentersDF$x, HRcentersDF$y,col="red", pch=19, cex=1)       ### Plot home range centres
points(sim.trap$x, sim.trap$y, col="blue", pch=19, cex=2)           ### Plot the traps

###### Calculate distances between stoat home range centres and each trap
#dist.trap<-crossdist(centre.stoat, sim.trap, win=win.sr)

trap.pos <- as.data.frame(trap.pos)
xtrap <- trap.pos[, 32]                      
ytrap <- trap.pos[, 33]

animX <- HRcentersDF[, 1]
animY <- HRcentersDF[, 2]

# a matrix with rows = animal and columns = trap, the value is the distance between them
distmat <- function(x1,y1,x2,y2){                     
  xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)     
  yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
  d <- t(sqrt(xd + yd)) 
  return(d)
}

dist.trap <- distmat(xtrap, ytrap, animX, animY)

#alternative way to calculate distance matrix using library pracma #gives identical matrix as above
trap.pos <- as.data.frame(trap.pos)
xtrap <- trap.pos[, 32]                      
ytrap <- trap.pos[, 33]
trap.loc <-cbind(xtrap,ytrap)

animX <- HRcentersDF[, 1]
animY <- HRcentersDF[, 2]
anim.loc <- cbind(animX,animY)

library(pracma)
dist.trap<-pdist2(trap.loc,anim.loc) #pdist2 is used instead of distmat as the later one is masked from other libraries
              
#### Simulate spatially-explicit stoat captures
g0<-0.04            ## Probability of capture at the centre of the home range

sigma<-400          ## Spatial decay parameter

p.bycatch<-0.01          ### Probability of by-catch

capt.trap<-matrix(0, ncol=tot.stoat, nrow=n.trap)

for (i in 1:tot.stoat){

        p.capt<-g0*exp((-dist.trap[i, ]*dist.trap[i, ])/(2*sigma^2))                   ## Probability of capture per stoat per night per trap

        #capt.trap[, i]<-p.capt

        tot.trap<-1-prod(1-(p.capt*(1-p.bycatch)))^n.day                     ### Probability of capture in any given trap after the total number of nights
        
        #tot.trap<-1-prod(1-(p.capt*(1-p.bycatch)^n.day))
        
        capt.trap[, i]<-rbinom(n.trap, size=1, prob=tot.trap)

         }


apply(capt.trap, 2, max)

tot.stoat-sum(apply(capt.trap, 2, max)==1)                      #### Stoats not caught

sum(apply(capt.trap, 2, max))/tot.stoat                         #### Proportion stoats within the active range of a trap
