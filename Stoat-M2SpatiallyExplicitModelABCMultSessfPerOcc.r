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
x.min<-min(stoat.trapping$X)-1000

x.max<-max(stoat.trapping$X)+1000

y.min<-min(stoat.trapping$Y)-1000

y.max<-max(stoat.trapping$Y)+1000


#### Establish the maximum number of stoats
area.sr<-49.8		#### Area (km2) of South Ronaldsay

#n.tot<-round(10*area.sr, digits=0)			### Maximum stoat abundance (8 is a maximum density/km2 in NZ)


#### ABC MODEL

se.stoat<-function(x, stoat.trapping){
 
	### Priors
	g0<-x[1]			### Probability of capture at the centre of the home range
	
	sigma<-x[2]		### Decay

	psi<-x[3]			### Presence in the population

	n.tot<-rpois(1, x[4])

    ### Stoat presence in the population

	z<-rbinom(n.tot, 1, psi)

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

				### Effective probability of capture (presence and p.capt)
    			comp.eff[i, j]<-1-p.capt[i, j]		#### Complement of the effective capture prob

				}
			}

    
        #### Probability of a capture in each trap

            p.trap<-stoat.cap<-rep(0, n.trap)

	        for (t in 1:n.trap){
		
				p.trap[t]<-1-prod(comp.eff[ , t])		### Effective probability of capturing at least one stoat in that trap

		        stoat.cap[t]<-rbinom(1, 1, p.trap[t])

	            }

            stoat.cap1[[1]]<-sum(stoat.cap)

            #### Probability of removal of each stoat

            for (n in 1:av.stoat){

                p.stoat<-1-prod(comp.eff[n, ])

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

	        ##### Observation model - distances between home ranges and each trap and probability of capture

            d<-p.capt<-comp.eff<-matrix(NA, ncol=n.trap, nrow=av.stoat)
    
	        for (i in 1:av.stoat){

			    for (j in 1:n.trap){
			
				### Distances between stoat home range centres and each trap
				d[i, j]<-sqrt((stoat.pos[i, 1]-trap.coord[j, 1])^2+(stoat.pos[i, 2]-trap.coord[j, 2])^2)

				### Probability of capture
				p.capt[i, j]<-g0*exp(-d[i, j]*d[i, j]/(2*sigma*sigma))

				### Effective probability of capture (presence and p.capt)
    			comp.eff[i, j]<-1-p.capt[i, j]		#### Complement of the effective capture prob

				}
			}

	        #### Probability of a capture in each trap

            p.trap<-stoat.cap<-rep(0, n.trap)

	        for (t in 1:n.trap){
		
				p.trap[t]<-1-prod(comp.eff[ , t])		### Effective probability of capturing at least one stoat in that trap

		        stoat.cap[t]<-rbinom(1, 1, p.trap[t])

	            }

            stoat.cap1[[r]]<-sum(stoat.cap)

            #### Probability of removal of each stoat

            for (n in 1:av.stoat){

                p.stoat<-1-prod(comp.eff[n, ])

                stoat.rem[n, r]<-rbinom(1, 1, p.stoat)

            }

    }

    ### Return the capture per trap and session

    return(c(unlist(stoat.cap1)))
}



### Define the number of simulations to run & values following a LHS design
n.sim<-1000

lhs1<-randomLHS(n.sim, 4)

g0<-qbeta(lhs1[,1], 1, 1)

sigma<-round(qunif(lhs1[,2], 0.1, 1000), digits=2)

psi<-qbeta(lhs1[,3], 1, 1)   

n.tot<-round(qunif(lhs1[,4], 40, 500), digits=1)           ### Minimum distance between traps

### First set a progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Empty matrices to store the results of the simulations
pred.vector<-matrix(0, nrow=n.session, ncol=n.sim)

prior.vector<-matrix(0, nrow=4, ncol=n.sim)

#### Timing the function
start.time<-Sys.time()

### Run the simulations
for (i in 1:n.sim){
     
    prior.vector[, i]<-c(g0[i], sigma[i], psi[i], n.tot[i])       #### Priors

    pred.vector[ , i]<-tryCatch(se.stoat(prior.vector[, i], stoat.trapping), error=function(e) NA)   ### Storing the summary statistics

    setTxtProgressBar(pb, i)

}

close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken


#### Identifying combinations of prior values that produced NAs in the summary statistics

sum(is.na(colMeans(pred.vector)))/n.sim                     ### Proportion not plausible

id.na<-which(is.na(colMeans(pred.vector)))

pred2<-pred.vector[ , -id.na]                                  ### Exclude simulations that produce NAs

ncol(pred2)

prior2<-prior.vector[ , -id.na]


#### Cross-validation using ridge regression

### Tolerances
tol<-c(0.05, 0.1, 0.5, 1)

loclinear.mean<-cv4abc(param=t(prior2), sumstat=t(pred2), nval=100, tols=c(0.1, 0.5, 1), statistic = "median", method="loclinear")

#ridge.mean<-cv4abc(param=t(prior2), sumstat=pred2, nval=100, tols=tol, statistic = "mean", method="ridge")

rejection.mean<-cv4abc(param=t(prior2), sumstat=t(pred2), nval=100, tols=tol, statistic = "median", method="rejection")

nnet.mean<-cv4abc(param=t(prior2), sumstat=t(pred2), nval=100, tols=tol, statistic = "median", method="neuralnet")

#### Summaries

summary(loclinear.mean)

summary(rejection.mean)

summary(nnet.mean)

#### Choosing one method and plotting
############################ BAYESIAN P-VALUES
sum.stat<-as.numeric(table(stoat.trapping$stoat.capt, stoat.trapping$VISIT_NO)[2,])

bpval.model<-gfit(target=sum.stat, sumstat=t(pred2), statistic=median, nb.replicate=200)

summary(bpval.model)

#### Inference

inf.abc<-abc(target=sum.stat, param=t(prior2), sumstat=t(pred2), tol=0.1, method="loclinear")

##### SUMMARIES
summary(inf.abc)

post.pred<-inf.abc$adj.values

summary(post.pred)

