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

### Coordinates of the island
x.min<-min(stoat.trapping$X)-10

x.max<-max(stoat.trapping$X)+10

y.min<-min(stoat.trapping$Y)-10

y.max<-max(stoat.trapping$Y)+10

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

se.stoat<-function(x, stoat.trapping, eff.stoat){
 
	### Priors
	g0<-plogis(x[1])			### Probability of capture at the centre of the home range
	
	sigma<-exp(x[2])		### Decay

	n.tot<-exp(x[3])

    ### Stoat presence in the population

	z<-rpois(1, n.tot)

    stoat.pos<-matrix(0, ncol=2, nrow=sum(z))

	#### Process model - stoat presence in the population and home range centre
	for (i in 1:sum(z)){

		#### Stoat home range centres in space
		stoat.pos[i, 1]<-runif(1, x.min, x.max)

		stoat.pos[i, 2]<-runif(1, y.min, y.max)

		}
    
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

            stoat.cap1[[1]]<-stoat.cap

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

            stoat.cap1[[r]]<-rep(0, n.trap)

            stoat.rem[, r]<-0
        }

    }

    ### Return the capture per trap and session

    return(c(unlist(stoat.cap1)))
}


### Define the number of simulations to run
n.its<-50000

lhs1<-randomLHS(n.its, 4)

g0<-qunif(lhs1[,1], -10, 0)

sigma<-round(qunif(lhs1[,2], -5, 7), digits=2)

n.tot<-round(qunif(lhs1[,3], log(45), log(500)), digits=1)           ### Minimum distance between traps

### First set a progress bar
pb<-txtProgressBar(min = 0, max = n.its, style = 3)

### Empty matrices to store the results of the simulations
pred.vector<-matrix(0, nrow=length(stoat.trapping$stoat.capt), ncol=n.its)

prior.vector<-matrix(0, nrow=3, ncol=n.its)

#### Timing the function
start.time<-Sys.time()

### Run the simulations
for (i in 1:n.its){
     
    prior.vector[, i]<-c(g0[i], sigma[i], n.tot[i])       #### Priors

    pred.vector[ , i]<-tryCatch(se.stoat(prior.vector[, i], stoat.trapping, eff.stoat), error=function(e) NA)   ### Storing the summary statistics

    setTxtProgressBar(pb, i)

}

close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken


#### Identifying combinations of prior values that produced NAs in the summary statistics

sum(is.na(colMeans(pred.vector)))/n.its                     ### Proportion not plausible

id.na<-which(is.na(colMeans(pred.vector)))

pred2<-pred.vector[ , -id.na]                                  ### Exclude simulations that produce NAs

ncol(pred2)

prior2<-prior.vector[ , -id.na]

### Exclude all 0s
id.zero<-which(colSums(pred.vector)==0)

id.zero

pred3<-pred.vector[ , -id.zero]                                  ### Exclude simulations that produce NAs

ncol(pred3)

ncol(pred3)/n.its

prior3<-prior.vector[ , -id.zero]

#### Cross-validation using ridge regression

### Tolerances
tol<-c(0.01, 0.05, 0.1, 0.5, 1)

loclinear.mean<-cv4abc(param=t(prior3), sumstat=t(pred3), nval=100, tols=tol, statistic = "median", method="loclinear")

ridge.mean<-cv4abc(param=t(prior3), sumstat=t(pred3), nval=100, tols=tol, statistic = "median", method="ridge")

rejection.mean<-cv4abc(param=t(prior3), sumstat=t(pred3), nval=100, tols=tol, statistic = "median", method="rejection")

nnet.mean<-cv4abc(param=t(prior3), sumstat=t(pred3), nval=100, tols=tol, statistic = "median", method="neuralnet")

#### Summaries

summary(loclinear.mean)

summary(rejection.mean)

summary(ridge.mean)

summary(nnet.mean)

#### Choosing one method and plotting
############################ BAYESIAN P-VALUES
sum.stat<-stoat.trapping$stoat.capt

bpval.model<-gfit(target=sum.stat, sumstat=t(pred3), statistic=median, nb.replicate=200)

summary(bpval.model)

#### Inference

inf.abc<-abc(target=sum.stat, param=t(prior3), sumstat=t(pred3), tol=0.5, method="rejection")

##### SUMMARIES
summary(inf.abc)

post.pred<-inf.abc$unadj.values

summary(post.pred)

write.table(post.pred, "d:/CONTAIN/Stoats RSPB/Modelling results/Inference-ModelSETrap.csv", sep=",")
