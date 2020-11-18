### Stoat data - Spatially-explicit model first occasion captures for South Ronaldsay
### Model 1- Trials with one trapping occasion only and no covariates for the location of the centre of the home range

set.seed(2012)

### Load the libraries
library(coda)
library(MASS) 
library(MCMCpack)
library(MuMIn)
library(jagsUI) 
library(ggmcmc)
library(corrplot)
library(nimble)
library(parallel)
library(doParallel)
library(foreach)
library(beepr)
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
library(rgeos)
library(abc)
library(lhs)
library(RColorBrewer)


########## Load the trapping data

stoat.trapping<-read.table("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SpatExp2/Trapping-data.csv", header=TRUE, sep=",")

summary(stoat.trapping)

stoat.trapping$Trap.ID<-as.numeric(stoat.trapping$TRAP_ID)

### Order by trapping session

stoat.trapping<-stoat.trapping[order(stoat.trapping$VISIT_NO), ]


#### Create a new variable indicating whether there was a stoat capture
stoat.trapping$stoat.capt<-ifelse(stoat.trapping$Catch_type=="S", 1, 0) 

stoat.trapping$stoat.capt

stoat.trapping$Catch_type

n.session<-max(stoat.trapping$VISIT_NO)

n.session

### Shape and coordinates of the island
sr<-readOGR("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SR_British_Grid.shp")

plot(sr)
   
########## Load the data on trapping effort per session

rem.data<-read.table("d:/CONTAIN/Stoats RSPB/Removal-SR.csv", sep=",", header=TRUE)

rem.data$Days<-round(rem.data$Mean_days, digits=0)

#### Effective trapping effort for stoats after accounting for captured rats and closed traps
eff.stoat<-rem.data$Days		 

eff.stoat


#### Establish the maximum number of stoats
area.sr<-49.8		#### Area (km2) of South Ronaldsay

#n.tot<-round(10*area.sr, digits=0)			### Maximum stoat abundance (8 is a maximum density/km2 in NZ)


#### ABC MODEL

se.stoat<-function(x, stoat.trapping, eff.stoat, sr){
 
	### Priors
	g0<-plogis(x[1])			### Probability of capture at the centre of the home range
	
	sigma<-x[2]		### Decay

	n.tot<-exp(x[3])

    ### Stoat presence in the population

	z<-rpois(1, n.tot)

    stoat.pos<-matrix(0, ncol=2, nrow=sum(z))

	#### Process model - stoat presence in the population and home range centre
    sim.hrcentre<-spsample(sr, n=z, "random", iter=20)

    stoat.pos<-as.matrix(coordinates(sim.hrcentre))
    
    ### Removed stoats - either 0 or 1

    stoat.rem<-matrix(0, nrow=sum(z), ncol=n.session)

    ### List to store traps with captures

    stoat.cap1<-vector("list")

    #### First session
    stoat.analys<-stoat.trapping[stoat.trapping$VISIT_NO==1, ]

    trap.coord<-matrix(NA, ncol=2, nrow=nrow(stoat.analys))

    trap.coord[, 1]<-stoat.analys$X

    trap.coord[, 2]<-stoat.analys$Y

    #### Number of traps

    n.trap<-nrow(stoat.analys)		### Number of traps (considering two traps per actual trap)

    ### Number of stoats available

    av.stoat<-sum(z)-sum(stoat.rem[ , 1])

    ##### Observation model - distances between home ranges and each trap and probability of capture

            d<-p.capt<-comp.eff<-matrix(NA, ncol=n.trap, nrow=av.stoat)
    
	        for (i in 1:av.stoat){

			    for (j in 1:n.trap){
			
				### Distances between stoat home range centres and each trap
				d[i, j]<-sqrt((stoat.pos[i, 1]-trap.coord[j, 1])^2+(stoat.pos[i, 2]-trap.coord[j, 2])^2)

				### Probability of capture
				p.capt[i, j]<-g0*exp(-d[i, j]*d[i, j]/(2*sigma*sigma))

				### Effective probability of capture (p.capt) in that trap after all trapping nights nights
    			comp.eff[i, j]<-1-p.capt[i, j]	#### Complement of the effective capture prob

				}
			}

    
        #### Probability of a capture in each trap

            p.trap<-p.night<-stoat.cap<-rep(0, n.trap)

	        for (t in 1:n.trap){
		
				p.night[t]<-1-prod(comp.eff[ , t])		### Effective probability of capturing at least one stoat in any trap one night

                p.trap[t]<-1-(1-p.night[t])^eff.stoat[1]

		        stoat.cap[t]<-rbinom(1, 1, p.trap[t])

	            }

            stoat.cap1[[1]]<-sum(stoat.cap)

            #### Probability of removal of each stoat

            for (n in 1:av.stoat){

                p.stoat1<-1-prod(comp.eff[n, ])

                p.stoat<-1-(1-p.stoat1)^eff.stoat[1]

                stoat.rem[n, 1]<-rbinom(1, 1, p.stoat)
            }

    ### Subsequent sessions
    for (r in 2:n.session){

            ### Choose the session data

            stoat.analys<-stoat.trapping[stoat.trapping$VISIT_NO==r, ]

            trap.coord<-matrix(NA, ncol=2, nrow=nrow(stoat.analys))

            trap.coord[, 1]<-stoat.analys$X

            trap.coord[, 2]<-stoat.analys$Y

            #### Number of traps

            n.trap<-nrow(stoat.analys)		### Number of traps (considering two traps per actual trap)

            ### Number of stoats available

            av.stoat<-sum(z)-sum(stoat.rem[ , 1:(r-1)])

            if(av.stoat!=0){

	        ##### Observation model - distances between home ranges and each trap and probability of capture

            d<-p.capt<-comp.eff<-matrix(NA, ncol=n.trap, nrow=av.stoat)
    
	        for (i in 1:av.stoat){

			    for (j in 1:n.trap){
			
				### Distances between stoat home range centres and each trap
				d[i, j]<-sqrt((stoat.pos[i, 1]-trap.coord[j, 1])^2+(stoat.pos[i, 2]-trap.coord[j, 2])^2)

				### Probability of capture
				p.capt[i, j]<-g0*exp(-d[i, j]*d[i, j]/(2*sigma*sigma))

                ### Effective probability of capture (p.capt) in that trap after all trapping nights nights
    			comp.eff[i, j]<-1-p.capt[i, j] 	#### Complement of the effective capture prob


				}
			}

	        #### Probability of a capture in each trap

            p.trap<-stoat.cap<-rep(0, n.trap)

	        for (t in 1:n.trap){
		
				p.night[t]<-1-prod(comp.eff[ , t])		### Effective probability of capturing at least one stoat in any trap one night

                p.trap[t]<-1-(1-p.night[t])^eff.stoat[r]

		        stoat.cap[t]<-rbinom(1, 1, p.trap[t])

	            }

            stoat.cap1[[r]]<-stoat.cap

            #### Probability of removal of each stoat

            for (n in 1:av.stoat){

                p.stoat1<-1-prod(comp.eff[n, ])

                p.stoat<-1-(1-p.stoat1)^eff.stoat[r]

                stoat.rem[n, 1]<-rbinom(1, 1, p.stoat)
            }
        
        }else{

            stoat.cap1[[r]]<-0

            stoat.rem[, r]<-0
        }

    }

    ### Return the capture per trap and session

    return(stoat.pos)
}


### Load the data
prior.vector<-read.table("d:/CONTAIN/Stoats RSPB/Modelling results/ParamEstimates-ABC.csv", header=TRUE, sep=",")

prior.vector

n.its<-nrow(prior.vector)

n.its

### First set a progress bar
pb<-txtProgressBar(min = 0, max = n.its, style = 3)

### Empty matrices to store the results of the simulations
pred.vector<-vector("list")

#### Timing the function
start.time<-Sys.time()

### Run the simulations
for (i in 1:n.its){
     
    pred.vector[[i]]<-tryCatch(se.stoat(as.numeric(prior.vector[i, ]), stoat.trapping, eff.stoat, sr=sr), error=function(e) NA)   ### Storing the summary statistics

    setTxtProgressBar(pb, i)

}

close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

pred.vector

#### All coordinates of home range centres together
tot.coord<-as.data.frame(do.call(rbind, pred.vector))

coordinates(tot.coord)<-c("x", "y")

kernel.ref<-kernelUD(tot.coord, h="href")

image(kernel.ref)

points(stoat.trapping$X[stoat.trapping$stoat.capt==1], stoat.trapping$Y[stoat.trapping$stoat.capt==1], pch=19, col="blue",
        cex.size=5)

### Kernel density estimation
tot.coord<-do.call(rbind, pred.vector)

kn.2d<-kde2d(tot.coord[, 1], tot.coord[, 2])

### Plotting 
filled.contour(kn.2d)

contour(kn.2d)

plot(sr, add=TRUE)

points(stoat.trapping$X, stoat.trapping$Y, pch=19, col="black")

points(stoat.trapping$X[stoat.trapping$stoat.capt==1], stoat.trapping$Y[stoat.trapping$stoat.capt==1], pch=19, col="blue")


### Transform to raster
t1<-raster(kn.2d$z*100, xmn=min(kn.2d$x), xmx=max(kn.2d$x), ymn=min(kn.2d$y), ymx=max(kn.2d$y))

pt.cuts<-seq(from=min(c(unlist(kn.2d$z*100))), to=max(c(unlist(kn.2d$z*100))), length.out=10)

col.pt<-colorRampPalette(c("yellow", "orange", "red"))

plot(t1, breaks=pt.cuts, col=col.pt(length(pt.cuts)))

points(stoat.trapping$X[stoat.trapping$stoat.capt==1], stoat.trapping$Y[stoat.trapping$stoat.capt==1], pch=19)


##### Half-normal detection & abundance estimate
par(mfrow=c(1, 2))

hist(exp(prior.vector[, 3]), xlab="Posterior abundance estimate", main="Estimated stoat abundance - South Ronaldsay")

abline(v=42, lwd=4)

### Halfnormal
dist<-seq(from=0, to=2000, length.out=2000)

av.eff<-plogis(colMeans(prior.vector)[1])*exp(-dist*dist/(2*colMeans(prior.vector)[2]*colMeans(prior.vector)[2]))

plot(av.eff~dist, type="l", lwd=4, ylab="Probability of capture per trap and night", xlab="Distance to the centre of the home range",
        main="Mean probability of capture")


