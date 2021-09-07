### Simulate spatially-explicit capture based on the distance to the centre of the home range
### Stoat distribution is random in entire area of South Ronaldsay
### Trap are distributed according the the deployment of March 2021
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
  tot.area<-sum(SR$SR_area)/10000 
  tot.stoat<-rpois(1, init.density*tot.area)
  
  #Home range centers of stoats
  HRcenters <- spsample(SR, n = tot.stoat, "random", iter = 20)
  Coordinates <- as.data.frame(HRcenters,rownames = NULL)
  colnames(Coordinates) <- NULL
  Coordinates <- as.matrix(Coordinates, rownames = FALSE)
  
  
  # Randomly distributed stoats are allowed to move prior to interacting with traps
  # Correlated Random Walk is defined by CRW- correlation between current and previous step, step - length of step,
  # horizon - scanning horizon (distance from which animal can "see" resistance of habitat), time - time of walk, 
  # resistance - taken values from 0 to1, indicating how a given habitat is suitable for stoat movement (1 means complete resistance, e.g. sea)
  # there is defined probability of transition between random walk (RW) and correlated random walk (CRW)
  
  
  walk<- vector("list") 
  
  
  for (i in 1:tot.stoat){
    
    sim <- simulate((species(state.RW() + state.CRW(CRW),
                             trans = transitionMatrix(0.005, 0.01))+step)*horizon, 
                    time, resist = resistanceSR, coords = Coordinates[i,1:2])
    
    walk[[i]] <- sim[time-1:time,1:2]  #function retains last position of each individual from correated random walk
  }
  
  ### All potential home range centers after movement defined as final points of random walk
  walk<-as.data.frame(walk)
  HRcenters_w <- walk
  HRcentersDF <- as.data.frame(HRcenters_w,rownames = NULL)
  
  ### Trap locations
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
    
    tot.trap<-1-prod((1-(p.capt*(1-p.bycatch)))^n.day)                     ### Probability of capture in any given trap after the total number of nights
    
    #tot.trap<-1-((1-p.capt)^n.day) 
    
    capt.trap[, i]<-rbinom(n.trap, size=1, prob=tot.trap)
    
  }
  
  ### Output 1 - Percentage of stoat population within the range of all traps
  
  perc.catch <- sum(apply(capt.trap, 2, max))/tot.stoat *100                        #### Proportion stoats within the active range of a trap
  
  ### Output 2 - number of not captured stoats
  
  ### Check the captures of stoats
  id.capt<-apply(capt.trap, 2, max)
  
  ### Extract the id of those stoats not captured
  id.notcapt<-which(id.capt==0)
  
  ### Number of stoats not caught
  post.notcapt<-anim.loc[id.notcapt, ]
  post.notcapt<-as.data.frame(post.notcapt)
  tot.stoat2<-length(post.notcapt$animX)
  
  
  ### Output of the function
  
  #return(tot.stoat2)
  
  out.sim<-round(c(perc.catch, tot.stoat2), digits=2)
  
  ### Function return percent of captured staots and number of not captured stoats
  return(out.sim)
  
  #Alternatively the output are coordinates of not captured stoats
  
  ### Check the captures of stoats
  #id.capt<-apply(capt.trap, 2, max)
  
  ### Extract the id of those stoats not captured
  #id.notcapt<-which(id.capt==0)
  
  ### Home range centres of stoats not caught
  #post.notcapt<-anim.loc[id.notcapt, ]
  
  ### Output of the function
  #return(post.notcapt)
 
}


