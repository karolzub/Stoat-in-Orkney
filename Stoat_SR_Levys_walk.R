### init.pos: a matrix containing the initial x and y position of the stoats at time 0 (column 1 = x, column 2 = y)
### nest.pos: a matrix containing the x and y positions of the kiwi nests (column 1 = x, column 2 = y)
### total.time: total time (days/months, etc)
### mean.dist: mean distance travelled by a stoat in a time step
### p.return: probability of returning to the previous position
### pref.index: stoat selectivity for kiwi nests (from 0 to 1): vector with values for each stost in the landscape
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
library(SiMRiv)

sronaldsay<-readOGR("C:/Orkney/Model/SRonaldsay.shp")
plot(sronaldsay)
tot.area<-sum(sronaldsay$SR_area)/10000

cover <- resistanceFromShape(sronaldsay, res = 50, background = 1, buffer = NA, margin = 200)

#init.density<-0.02
tot.stoat<-rpois(1, init.density*tot.area)
HRcenters <- spsample(sronaldsay, n = tot.stoat, "random", iter = 20)
init.pos <- as.data.frame(HRcenters)
write.table(init.pos, file="C:/Orkney/Model/dat2.txt", append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

#n.trap <- 500
centre.trap.rand <- spsample(sronaldsay, n = n.trap, "random", iter = 20)
nest.pos <- as.data.frame(centre.trap.rand)

total.time <- 60

#p.return <- 1

#mean.dist <- 300

stoat.ibm<-function(init.pos, nest.pos, total.time, init.density, n.kiwi, mean.dist, p.return){
        require(adehabitatHR)
        require(sp)
        require(raster)
	
		output<-vector("list")

		### Total number of stoats simulated
		n.stoat<-nrow(init.pos)
	
        #### Total number of kiwi nests
        n.kiwi<-nrow(nest.pos)

        ####### KIWI NESTS EACH STEP
        surv.nest<-array(0, c(n.kiwi, 2, total.time))
        surv.nest[, ,1]<-nest.pos		### The first matrix is the position of the nests

        ###### STOAT MOVEMENT
    	### An array containing the (x,y) coordinates of each stoat at each time step. Columns 1 and 2 are the x and y-coordinates, respectively. 
        ### Column 3 is the status of the individual (1=moving, 2=established, 3=dead)
        movement<-array(0, c(total.time, 2, n.stoat))
		    movement[1,1:2,]<-round(as.vector(t(init.pos)), digits=0)			### The first row of the movement array are the initial coordinates

		### Loop to make the simulate the search dynamics
		### n.inds: number of stoats
		
		for (j in 2:total.time){			### total.time: total time steps


                        for (i in 1:n.stoat){				
					
                                   ### Make the stoat move  using a Levy walk
											
									step.size<-rexp(1, 1/mean.dist)     		### Movement per time step in metres
									theta<-runif(1, 0, 2*pi)						### Movement angle
									movement[j,1, i]<-round(movement[j-1, 1, i]+step.size*sin(theta), digits=0)			### New x-position
									movement[j,2, i]<-round(movement[j-1, 2, i]+step.size*cos(theta), digits=0)			### New y-position
											
									}

                        }

}

write.table(movement, file="C:/Orkney/Model/dat1.txt", append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

                    
####### Estimate 95% MCP stoat home ranges
                    data.homerange<-data.frame(x=as.vector(movement[,1,]), y=as.vector(movement[,2,]), id=rep(c(1:n.stoat), each=total.time))

                    coordinates(data.homerange)<-data.homerange[, c('x', 'y')]

                    mcp.95<-mcp(data.homerange[, 3], percent = 95)

                    hr.area<-as.data.frame(mcp.95)

                  
                    ###### Loop for estimating overlap between each home range and each kiwi nest

                    over.nest<-matrix(0, ncol=n.kiwi, nrow=n.stoat)

                    for (i in 1:n.stoat){

                            stoat.coord<-mcp.95@polygons[[i]]@Polygons[[1]]@coords				### Coordinates of the home range of stoat 1

					        over.nest1<-point.in.polygon(nest.pos[,1], nest.pos[,2], as.vector(stoat.coord[,1]), as.vector(stoat.coord[,2]), mode.checked=FALSE)			#### Check which kiwi nests are within the HR

					        over.nest[i,]<-ifelse(over.nest1>0, 1, 0)						#### Nests within the home range are 1 other 0

                         }

                    #### Stoat specialisation
                    over<-apply(over.nest, 1, max)

                    id.overlap<-which(over==1)

                    tot.overlap<-length(id.overlap)

                    ##### Selectivity index & depredation

                    pref.index<-rbeta(tot.overlap, 5, 5)

                    prob.pred<-over.nest

                    prob.pred[id.overlap, ]<-over.nest[id.overlap, ]*pref.index


                    ######## TOTAL PROBABILITY OF PREDATION FOR EACH NEST

                    p.surv<-apply(1-prob.pred, 2, prod)           #### Probability of nest survival


                    nest.surv<-rbinom(length(p.surv), size=1, prob=p.surv)      #### The nest survived or not (from a Binomial)
         
                    #### Export results


                    output<-c(n.kiwi, sum(nest.surv), n.stoat, median(p.return), median(hr.area$area), median(pref.index), sum(pref.index>0.74),        
                                       median(mean.dist), sum(apply(over.nest, 2, max))/n.kiwi)         
    
                #### Return output: initial number of nests; number of nests not predated; proportion surviving; stoat density; probability of returning; detection distance; median home range size; proportion of stoats whose home range contain at least one kiwi nest
   

                return(output)
                    
                    
                    output<-c(n.kiwi, n.stoat, median(hr.area$area), median(mean.dist), sum(apply(over.nest, 2, max))/n.kiwi)         
                    




