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
        HRcenters <- spsample(main, n = tot.stoat, "random", iter = 20)
        HRcentersDF <- as.data.frame(HRcenters)
        #plot(HRcenters, add = TRUE, col="red", pch=19, cex=1)
        
        #### TRAPS
        centre.trap.rand <- spsample(main, n = 5000, "random", iter = 20)
        #centre.trap <- traps_bord[sample(nrow(traps_bord), 5000), ]
        centre.trap <- as.data.frame(centre.trap.rand)
        coordinates(centre.trap)  <-  c("x", "y")
        #plot(centre.trap, add = TRUE, col="blue", pch=19, cex=1)
        ### Number of traps
        n.trap<-length(centre.trap$x)
        #n.trap<-length(centre.trap$coords.x1)
        n.trap
        
        #n.day <- 20
        
        #### Real traps
        trap.pos <- as.data.frame(centre.trap)
        xtrap <- trap.pos$x                      
        ytrap <- trap.pos$y
        
        sim.trap<-ppp(xtrap, ytrap, check=FALSE)  	#### Real traps
        
        trap.dist <- mean(nndist(sim.trap))  ### Check mean distance between traps (metres)
        
        #### Calculate distances between stoat home range centres and each trap
        
        xtrap <- trap.pos$x                     
        ytrap <- trap.pos$y
        trap.loc <-cbind(xtrap,ytrap)
        
        animX <- HRcentersDF[, 1]
        animY <- HRcentersDF[, 2]
        anim.loc <- cbind(animX,animY)
        
        
        dist.trap<-pdist2(anim.loc, trap.loc) #used pdist2 instead of distmat, which is disabled by other functions
        
        #a matrix with rows = animal and columns = trap, the value is the distance between them
        #distmat <- function(x1,y1,x2,y2){                     
               #xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)     
                #yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
               #d <- t(sqrt(xd + yd)) 
                #return(d)
        #}
        
        #dist.trap <- distmat(xtrap, ytrap, animX, animY)    
        
        #### Simulate spatially-explicit stoat captures
        #g0<-0.04            ## Probability of capture at the centre of the home range
        
        #sigma<-600          ## Spatial decay parameter
        
        #p.bycatch<-0.008        ### Probability of by-catch per trap day
        
        capt.trap<-matrix(0, ncol=tot.stoat, nrow=n.trap)
        
        for (i in 1:tot.stoat){
                
                p.capt<-g0*exp((-dist.trap[i, ]^2)/2/(sigma^2))                   ## Probability of capture per stoat per night per trap
                
                #capt.trap[, i]<-p.capt
                
                tot.trap<-1-prod(1-(p.capt*(1-p.bycatch)^n.day))                     ### Probability of capture in any given trap after the total number of nights
                
                #tot.trap<-1-((1-p.capt)^n.day) 
                
                capt.trap[, i]<-rbinom(n.trap, size=1, prob=tot.trap)
                
        }
     
        #apply(capt.trap, 2, max)
        
        #tot.stoat-sum(apply(capt.trap, 2, max)==1)                      #### Stoats not caught
        ### Percentage of stoat population within the range of each camera trap
        
        perc.catch <- sum(apply(capt.trap, 2, max))/tot.stoat *100                        #### Proportion stoats within the active range of a trap
        
        perc.catch
        
        
        
        ### Output of the function - returns the mean home range size (ha) and the percentage of stoats caught

        #out.sim<-round(c(mean(dist.trap), perc.catch), digits=2)
        out.sim <- perc.catch
        return(out.sim)

}

perc.catch



