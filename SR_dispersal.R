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

stoat.dist<-function(CRW, step, horizon, time, 
                     init.density, n.day, g0, sigma, p.bycatch)
   
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
        tot.area<-sum(SR$SR_area)/10000 
        tot.stoat<-rpois(1, init.density*tot.area)
        
        #CRW<-0.98
        #step<-100
        #horizon<-300
        #time<-1000
        time1<- round(time*0.9)
        HRcenters <- spsample(SR, n = tot.stoat, "random", iter = 20)
        Coordinates <- as.data.frame(HRcenters,rownames = NULL)
        colnames(Coordinates) <- NULL
        Coordinates <- as.matrix(Coordinates, rownames = FALSE)
        #resistanceSR <- resistanceFromShape("C:/Orkney/Model/SR_habitats.shp", res = 10, field = "Habitat", mapvalues = c("coast" = 0.2, "grassland" = 0.5, "meadow" = 0.8, "heatland" = 0.25,"marine" = 0.1, "road" = 0.1, "town" = 0.9,"water" = 1.00, "wood" = 0.75 ), background = 1.00, margin = 1000)
        
        
        #levy.walker <- species(state.RW() + state.CRW(CRW), trans = transitionMatrix(0.005, 0.01)) + 25
        #levy.walker <- (levy.walker + step) * horizon
        #sim <- simulate(levy.walker, time, resist = resistanceSR, coords = Coordinates[1:n.stoat,1:2])
        
        #walk <- matrix(0, ncol=2, nrow=tot.stoat)
        
        #for (i in 1:tot.stoat){
          
          #sim <- simulate((species(state.RW() + state.CRW(CRW),trans = transitionMatrix(0.005, 0.01))+step)*horizon, time, resist = resistanceSR, coords = Coordinates[i,1:2])
          #walk[i,1:2]<-sim[time,1:2]
        #}
        
        walk<- vector("list") 
        
        for (i in 1:tot.stoat){
          
          sim <- simulate((species(state.RW() + state.CRW(CRW),
                                  trans = transitionMatrix(0.005, 0.01))+step)*horizon, 
                         time, resist = resistanceSR, coords = Coordinates[i,1:2])
          
          walk[[i]] <- sim[time1:time,1:2]
        }
        
        #coordinates<-do.call(rbind.data.frame, walk)
        
        ### All potential home range centres
        walk<-as.data.frame(walk)
        #HRcenters_w <- walk
        HRcenters_w <- sample_n(walk, tot.stoat, replace = TRUE)
        #HRcenters <- coordinates
        HRcentersDF <- as.data.frame(HRcenters_w,rownames = NULL)
        
        ### All potential trap locations
        #centre.trap.rand <- spsample(access, n = n.trap, "random", iter = 20)
        centre.trap.rand <- SR_traps
        centre.trap <- as.data.frame(centre.trap.rand)
        xtrap <- centre.trap$coords.x1                     
        ytrap <- centre.trap$coords.x2
        trap.loc <-cbind(xtrap,ytrap)
        n.trap<-length(centre.trap$coords.x1)
        
        ### Stoat home range cencters locations
        animX <- HRcentersDF[, 1]
        animY <- HRcentersDF[, 2]
        anim.loc <- cbind(animX,animY)
        
        #### Calculate distances between stoat home range centres and each trap
        
        dist.trap<-pdist2(anim.loc, trap.loc) #used pdist2 instead of distmat, which is disabled by other functions
        
        #### Simulate spatially-explicit stoat captures
        
        capt.trap<-matrix(0, ncol=tot.stoat, nrow=n.trap)
        
        for (i in 1:tot.stoat){
          
          p.capt<-g0*exp((-dist.trap[i, ]^2)/2/(sigma^2))                   ## Probability of capture per stoat per night per trap
          
          #capt.trap[, i]<-p.capt
          
          tot.trap<-1-prod(1-(p.capt*(1-p.bycatch)^n.day))                     ### Probability of capture in any given trap after the total number of nights
          
          #tot.trap<-1-((1-p.capt)^n.day) 
          
          capt.trap[, i]<-rbinom(n.trap, size=1, prob=tot.trap)
          
        }
     
        apply(capt.trap, 2, max)
        
        tot.stoat-sum(apply(capt.trap, 2, max)==1)                      #### Stoats not caught
        ### Percentage of stoat population within the range of each camera trap
        
        numb.catch <- sum(apply(capt.trap, 2, max))                        #### Proportion stoats within the active range of a trap
        
        ### Output of the function - returns the mean home range size (ha) and the percentage of stoats caught
        
        #out.sim<-round(c(mean(dist.trap), perc.catch), digits=2)
        #out.sim <- perc.catch
        #return(out.sim)
        
        ### Trap locations in 500 m buffer along the coast line
        centre.trap.rand2 <- SR_traps_out
        centre.trap2 <- as.data.frame(centre.trap.rand2)
        xtrap2 <- centre.trap2$coords.x1                     
        ytrap2 <- centre.trap2$coords.x2
        trap.loc2 <-cbind(xtrap2,ytrap2)
        n.trap2<-length(centre.trap2$coords.x1)
        
        dist.trap2<-pdist2(anim.loc, trap.loc2) #used pdist2 instead of distmat, which is disabled by other functions
        
        capt.trap2<-matrix(0, ncol=tot.stoat, nrow=n.trap2)
        
        for (i in 1:tot.stoat){
          
          p.capt2<-g0*exp((-dist.trap2[i, ]^2)/2/(sigma^2))                   ## Probability of capture per stoat per night per trap
          
          tot.trap2<-1-prod(1-(p.capt2*(1-p.bycatch)^n.day))                     ### Probability of capture in any given trap after the total number of nights
          
          capt.trap2[, i]<-rbinom(n.trap2, size=1, prob=tot.trap2)
          
        }
        
        apply(capt.trap2, 2, max)
        
        numb.catch2 <- sum(apply(capt.trap2, 2, max)) 
        
        #prop.500 <- (numb.catch2/numb.catch)*100
        
        #perc.catch2 <- sum(apply(capt.trap2, 2, max))/tot.stoat *100                        #### Proportion stoats within the active range of a trap
        
        
        out.sim<-round(c(perc.catch2,perc.catch), digits=2)
        out.sim <- perc.catch
        return(out.sim)
        
        
        ### Check the captures of stoats
        #id.capt<-apply(capt.trap2, 2, max)
        
        ### Extract the id of those stoats captured
        
        #id.capt<-which(id.capt==1)
        ### Home range centres of stoats caught
        
        #post.capt<-anim.loc[id.capt, ]
        ### Output of the function
       
        #return(post.capt)
}

