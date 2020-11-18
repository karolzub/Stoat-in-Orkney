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
library(EasyABC)
library(lhs)
library(abc)

    ########## Load the trapping data

stoat.trapping<-read.table("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SpatExp2/Trapping-data.csv", header=TRUE, sep=",")

stoat.trapping$Trap.ID<-as.numeric(stoat.trapping$TRAP_ID)

### Order by trapping session

stoat.trapping<-stoat.trapping[order(stoat.trapping$VISIT_NO), ]


#### Create a new variable indicating whether there was a stoat capture
stoat.trapping$stoat.capt<-ifelse(stoat.trapping$Catch_type=="S", 1, 0) 

n.session<-max(stoat.trapping$VISIT_NO)

### Shape and coordinates of the island
sr<-readOGR("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SR_British_Grid.shp")

########## Load the data on trapping effort per session

rem.data<-read.table("d:/CONTAIN/Stoats RSPB/Removal-SR.csv", sep=",", header=TRUE)

rem.data$Days<-round(rem.data$Mean_days, digits=0)


#### Effective trapping effort for stoats after accounting for captured rats and closed traps
eff.stoat<-rem.data$Days		 

#### Establish the maximum number of stoats
area.sr<-49.8		#### Area (km2) of South Ronaldsay

#n.tot<-round(10*area.sr, digits=0)			### Maximum stoat abundance (8 is a maximum density/km2 in NZ)


#### ABC MODEL

se.stoat<-function(x){
 
	### Priors
	g0<-plogis(x[1])			### Probability of capture at the centre of the home range
	
	sigma<-x[2]		### Decay

	n.tot<-exp(x[3])

    ### Stoat presence in the population

	z<-rpois(1, n.tot)

	#### Process model - stoat presence in the population and home range centre
    sim.hrcentre<-spsample(sr, n=z, "random", iter=20)

    stoat.pos<-as.matrix(coordinates(sim.hrcentre))
    
    ### Removed stoats - either 0 or 1

    stoat.rem<-matrix(0, nrow=sum(z), ncol=n.session)

    ### List to store traps with captures

    stoat.cap1<-stoat.cap2<-vector("list")

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

            stoat.cap1[[1]]<-stoat.cap

            stoat.cap2[[1]]<-sum(stoat.cap)

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

            stoat.cap2[[r]]<-sum(stoat.cap)

            #### Probability of removal of each stoat

            for (n in 1:av.stoat){

                p.stoat1<-1-prod(comp.eff[n, ])

                p.stoat<-1-(1-p.stoat1)^eff.stoat[r]

                stoat.rem[n, 1]<-rbinom(1, 1, p.stoat)
            }
        
        }else{

            stoat.cap1[[r]]<-rep(0, n.trap)

            stoat.cap2[[r]]<-0

            stoat.rem[, r]<-0
        }

    }

    ### Return the capture per trap and session
    out.sim<-c(c(unlist(stoat.cap1)), c(unlist(stoat.cap2)))

    norm.sim<-(out.sim1-mean(out.sim1))/(sd(out.sim1)/sqrt(length(out.sim1)))

    return(norm.sim)
}



########## Summary statistics to target
sum.stat<-c(stoat.trapping$stoat.capt, as.numeric(table(stoat.trapping$stoat.capt, stoat.trapping$VISIT_NO)[2,]))

sum.stat

norm.stat<-(sum.stat-mean(sum.stat))/(sd(sum.stat)/sqrt(length(sum.stat)))

norm.stat

### Vector with the priors in order: parameters in the regressions
prior.list<-list(c("unif", -10, 0), c("normal", 255, 70), c("unif", log(45), log(500)))

### Define the number of ABC iterations
n.its<-100

### ABC algorithm
abc.rej<-ABC_rejection(model=se.stoat, prior=prior.list, nb_simul=n.its, summary_stat_target=sum.stat, tol=1, 
            progress_bar=TRUE)

colMeans(abc.rej$stats)

plot(sum.stat, apply(abc.rej$stats, 2, median), xlab="True value of the summary statistics", ylab="ABC simulated value of the summary statistics")

abline(a=0, b=1, lwd=4)


### Sequential (particle filter)
abc.pf<-ABC_sequential(method="Lenormand", model=se.stoat, prior=prior.list, nb_simul=n.its, summary_stat_target=sum.stat,
            progress_bar=TRUE, alpha=0.1)

colMeans(abc.pf$stats)

plot(sum.stat, apply(abc.pf$stats, 2, median), xlab="True value of the summary statistics", 
                ylab="ABC simulated value of the summary statistics", pch=19, cex=2)

abline(a=0, b=1, lwd=4)


############################ Bayesian p-value
bpval.model<-gfit(target=sum.stat, sumstat=abc.pf$stats, statistic=median, nb.replicate=200)

summary(bpval.model)
