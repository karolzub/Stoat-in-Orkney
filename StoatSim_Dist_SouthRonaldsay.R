### R script - Model to simulate home ranges and traps within the home range.
### Data on trap deployment for South Ronaldsay (fixed number of traps = 406)
### Parameters: init.density: stoat density per ha
### tot.area: area (ha)
### n.trap: number of traps
### n.day: number of trapping days
### cor.mov: correlation in the movements of stoats between one day and the next
### min.det: minimum probability of capture
### max.det: maximum probability of capture

stoat.dist<-function(init.density, n.trap, trap.dist, n.day, g0, sigma, p.bycatch)
        
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
        
        ## Total number of stoats
        #init.density<-0.02			#### Density in individuals per ha
        tot.stoat<-rpois(1, init.density*tot.area)
        
        ### All potential home range centres
        HRcenters <- spsample(sr, n = tot.stoat, "random", iter = 20)
        HRcentersDF <- as.data.frame(HRcenters)
        
        #### TRAPS
        n.trap<-length(trap.pos$ID)		### Number of traps
        #n.trap
        
        n.day <- 20
        
        #### Real traps
        trap.pos <- as.data.frame(trap.pos)
        xtrap <- trap.pos$coords.x1                      
        ytrap <- trap.pos$coords.x2
        
        sim.trap<-ppp(xtrap, ytrap, check=FALSE)  	#### Real traps
        
        trap.dist <- mean(nndist(sim.trap))  ### Check mean distance between traps (metres)
        
        #### Calculate distances between stoat home range centres and each trap
        
        xtrap <- trap.pos$coords.x1                      
        ytrap <- trap.pos$coords.x2
        trap.loc <-cbind(xtrap,ytrap)
        
        animX <- HRcentersDF[, 1]
        animY <- HRcentersDF[, 2]
        anim.loc <- cbind(animX,animY)
        
        
        dist.trap<-pdist2(anim.loc, trap.loc) #used pdist2 instead of distmat, which is disabled by other functions
        
        #a matrix with rows = animal and columns = trap, the value is the distance between them
        #distmat <- function(x1,y1,x2,y2){                     
               # xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)     
                #yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
               # d <- t(sqrt(xd + yd)) 
                #return(d)
        #}
        
        #dist.trap <- distmat(xtrap, ytrap, animX, animY)    
        
        #### Simulate spatially-explicit stoat captures
        g0<-0.04            ## Probability of capture at the centre of the home range
        
        sigma<-600          ## Spatial decay parameter
        
        p.bycatch<-0.05         ### Probability of by-catch per trap day
        
        capt.trap<-matrix(0, ncol=tot.stoat, nrow=n.trap)
        
        for (i in 1:tot.stoat){
                
                p.capt<-g0*exp((-dist.trap[i, ]^2)/2/(sigma^2))                   ## Probability of capture per stoat per night per trap
                
                #capt.trap[, i]<-p.capt
                
                tot.trap<-1-prod(1-(p.capt*(1-p.bycatch)^n.day))                     ### Probability of capture in any given trap after the total number of nights
                
                capt.trap[, i]<-rbinom(n.trap, size=1, prob=tot.trap)
                
        }
     
        #apply(capt.trap, 2, max)
        
        #tot.stoat-sum(apply(capt.trap, 2, max)==1)                      #### Stoats not caught
        ### Percentage of stoat population within the range of each camera trap
        
        perc.catch <- sum(apply(capt.trap, 2, max))/tot.stoat *100                        #### Proportion stoats within the active range of a trap
        
        
        
        
        ### Output of the function - returns the mean home range size (ha) and the percentage of stoats caught

        #out.sim<-round(c(mean(dist.trap), perc.catch), digits=2)
        out.sim <-perc.catch
        return(out.sim)

}
