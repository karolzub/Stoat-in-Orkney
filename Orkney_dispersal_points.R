### Simulate spatially-explicit capture based on the distance to the centre of the home range
### Stoat distribution is random in entire area of Mainlad
### Trap are distributed either randomly in entire area (file: main)
### or randomly in areas with secured access (file: access)
### and only on the properties borders in no-access areas(file: trap_random)
### There is also design with traps on the properties borders only (file: trap_bord)
### Parameters: init.density: stoat density per ha
### tot.area: area (ha)
### n.trap: number of traps
### n.day: number of trapping days
### g0: probability of capture at the centre of the home range
### sigma: spatial decay parameter
### p.bycatch: probability of by-catch per trap day
### tot.area: total area of the site (ha)
### spat.layer: spatial layer of the site; shapefile format 
### traps_rand: x and y coordinates of the trap positions

#n.day<-28
#init.density<-0.03
#trap.dist<-250
#g0<-0.05
#sigma<-400
#p.bycatch<-0.01

### This script can be used for optimal camera-traps deployment using existing design
### This script can be applied for any area with known trap deployment

stoat.dist<-function(CRW, step, horizon, time, init.density)
   
        {

### Packages required.
        library(MASS)
        library(MBESS)
        library(adehabitatHR)
        library(sp)
        library(spatstat)
        library(raster)
        library(geosphere)
        library(maptools)
        library(rgdal)
        library(sf)
        library(beepr)
        library(rgeos)
        library(pracma)
        library(dplyr)
        library(SiMRiv)
        library(mc2d)
  
  
        ## Total number of stoats
        #main <- readOGR("C:/Orkney/Model/Main.shp")
        tot.area<-sum(Orkney$Area)/10000 
        tot.stoat<-rpois(1, init.density*tot.area)
        
        #CRW<-0.98
        #step<-100
        #horizon<-300
        #time<-1000
        HRcenters <- spsample(Orkney, n = tot.stoat, "random", iter = 20)
        Coordinates <- as.data.frame(HRcenters,rownames = NULL)
        colnames(Coordinates) <- NULL
        Coordinates <- as.matrix(Coordinates, rownames = FALSE)
        #resistanceOrkney <- resistanceFromShape("C:/Orkney/Model/Orkney habitats.shp", res = 10, field = "Habitat", mapvalues = c("coast" = 0.1, "grassland" = 0.8, "meadow" = 0.8, "heatland" = 0.8,"marine" = 0.8, "road" = 0.8, "town" = 0.9,"water" = 1.00, "wood" = 0.8 ), background = 1.00, margin = 1000)
        
        
        walk<- vector("list") 
        
        for (i in 1:tot.stoat){
          
          sim <- simulate((species(state.RW() + state.CRW(CRW),
                                   trans = transitionMatrix(0.005, 0.01))+step)*horizon, 
                          time, resist = resistanceOrkney, coords = Coordinates[i,1:2])
          
          walk[[i]] <- sim[ ,1:2]
        }
        
        coordinates<-do.call(rbind.data.frame, walk)
        
        return(coordinates)
}

