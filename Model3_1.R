### Simulate spatially-explicit capture based on the distance of trap to the center of the stoat home range
### Data on trap deployment for Orkney
### Additionally there could be simulated random distribution of traps either equal number to the real one or two-times higher
### Parameters: init.density: initial stoat density per ha
### tot.area: area (ha)
### n.trap: number of traps
### n.day: number of trapping days
### g0: probability of capture at the centre of the home range
### sigma: spatial decay parameter
### p.bycatch: probability of by-catch per trap day

### This script can be used for optimal trap deployment using existing design
### This script can be applied for any area with known trap deployment

stoat.dist<-function(init.density, n.day, g0, sigma, p.bycatch)
        
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
        tot.stoat<-rpois(1, init.density*tot.area)
        
        ### All potential home range centres
        HRcenters <- spsample(orkney, n = tot.stoat, "random", iter = 20)
        HRcentersDF <- as.data.frame(HRcenters)
        
        #### Real traps locations only within areas with access
        trap.pos <- as.data.frame(trap.pos)
        xtrap <- trap.pos$coords.x1                      
        ytrap <- trap.pos$coords.x2
        trap.loc <-cbind(xtrap,ytrap)
        n.trap<-length(trap.pos$ID)
        
        ### All potential trap locations
        #centre.trap.rand <- spsample(orkney, n = n.trap, "random", iter = 20)
        #centre.trap <- as.data.frame(centre.trap.rand)
        #xtrap <- centre.trap$x                      
        #ytrap <- centre.trap$y
        #trap.loc <-cbind(xtrap,ytrap)
        
        ### Stoat home range cencters locations
        animX <- HRcentersDF[, 1]
        animY <- HRcentersDF[, 2]
        anim.loc <- cbind(animX,animY)
        
        #### Calculate distances between stoat home range centres and each trap
        
        dist.trap<-pdist2(anim.loc, trap.loc) #used pdist2 instead of distmat, which is disabled by other functions
        
        #### Alternative method to calculate distances between stoat home range centres and each trap
        #a matrix with rows = animal and columns = trap, the value is the distance between them
        #distmat <- function(x1,y1,x2,y2){                     
                #xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)     
                #yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
                #d <- t(sqrt(xd + yd)) 
                #return(d)
        #}
        
        #dist.trap <- distmat(xtrap, ytrap, animX, animY)    
        
        #### Simulate spatially-explicit stoat captures
        
        capt.trap<-matrix(0, ncol=tot.stoat, nrow=n.trap)
        
        for (i in 1:tot.stoat){
          
          p.capt<-g0*exp((-dist.trap[i, ]^2)/2/(sigma^2))                   ## Probability of capture per stoat per night per trap
          
          #capt.trap[, i]<-p.capt
          
          tot.trap<-1-prod((1-(p.capt*(1-p.bycatch)))^n.day)                     ### Probability of capture in any given trap after the total number of nights
          
          #tot.trap<-1-((1-p.capt)^n.day) 
          
          capt.trap[, i]<-rbinom(n.trap, size=1, prob=tot.trap)
          
        }
     
        ### Percentage of stoat population within the range of all traps
        
        #perc.catch <- sum(apply(capt.trap, 2, max))/tot.stoat *100                        #### Proportion stoats within the active range of a trap
        
        ### Output of the function - returns the percentage of stoats caught

        #out.sim <- perc.catch
        #return(out.sim)


### Output - initial number of staots and number of not captured stoats

### Check the captures of stoats
#id.capt<-apply(capt.trap, 2, max)

### Extract the id of those stoats not captured
#id.notcapt<-which(id.capt==0)

### Home range centres of stoats not caught (coordinates)
#<-anim.loc[id.notcapt, ]
#post.notcapt<-as.data.frame(post.notcapt)
#tot.stoat2<-length(post.notcapt$animX)


### Output of the function

#return(tot.stoat2)

#out.sim<-round(c(perc.catch, tot.stoat2), digits=2)


#return(out.sim)



#Alternatively the output are coordinates of not captured stoats

### Check the captures of stoats
id.capt<-apply(capt.trap, 2, max)

### Extract the id of those stoats not captured
id.notcapt<-which(id.capt==0)

### Home range centres of stoats not caught
post.notcapt<-anim.loc[id.notcapt, ]

### Output of the function
return(post.notcapt)
}