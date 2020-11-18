### Data simulation to verify the MCMC model
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

########## Load the trapping data

stoat.trapping<-read.table("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SpatExp2/Trapping-data.csv", header=TRUE, sep=",")

summary(stoat.trapping)

### For this example, use first trapping session only

stoat.analys<-stoat.trapping[stoat.trapping$VISIT_NO==1, ]

summary(stoat.analys)

#### Sample size, number of traps, and other for the Bayesian analysis

n.trap<-nrow(stoat.analys)		### Number of traps (considering two traps per actual trap)

n.trap

### Create a matrix to store the trap coordinates
trap.coord<-matrix(0, ncol=2, nrow=n.trap)

trap.coord[, 1]<-stoat.analys$X

trap.coord[, 2]<-stoat.analys$Y

x.min<-min(trap.coord[, 1])-1000

x.max<-max(trap.coord[, 1])+1000

y.min<-min(trap.coord[, 2])-1000

y.max<-max(trap.coord[, 2])+1000

##### Simulate stoats and capture data

n.stoat<-50                     ### Number of stoats

pos.init<-matrix(0, ncol=2, nrow=n.stoat)       #### Location of their home range centres (random uniform)

pos.init[, 1]<-runif(n.stoat, x.min, x.max)

pos.init[, 2]<-runif(n.stoat, y.min, y.max)

### Transform to point patters
trap.loc<-ppp(trap.coord[, 1], trap.coord[, 2], check=FALSE)

centre.stoat<-ppp(pos.init[, 1], pos.init[, 2], check=FALSE)

###### Calculate distances between stoat home range centres and each trap
dist.trap<-crossdist(centre.stoat, trap.loc)


#### Simulate spatially-explicit stoat captures
g0<-0.04            ## Probability of capture at the centre of the home range

sigma<-600          ## Spatial decay parameter

p.bycatch<-0.27          ### Probability of by-catch

capt.trap<-matrix(0, ncol=n.stoat, nrow=n.trap)

for (i in 1:n.stoat){

        p.capt<-g0*exp((-dist.trap[i, ]*dist.trap[i, ])/(2*sigma^2))                   ## Probability of capture per stoat per night per trap

        #capt.trap[, i]<-p.capt

        capt.trap[, i]<-rbinom(n.trap, size=1, prob=p.capt*p.bycatch)

         }

stoat.capt<-ifelse(rowSums(capt.trap)!=0, 1, 0)          #### Strike a capture if at least one has been caught in the trap


############################### MODELLING
### Simulate the process of fitting the model as if we didn't just simulate the data
n.tot<-100      #### Maximum total of stoats

#### NIMBLE BAYESIAN MODEL

se.stoat<-nimbleCode({
 
	### Priors
	g0~dbeta(1, 1)			### Probability of capture at the centre of the home range
	
	sigma~dunif(0.001, 1000)		### Decay

	psi~dbeta(1, 1)			### Presence in the population

	
	#### Process model - stoat presence in the population and home range centre
	for (i in 1:n.tot){

		### Stoat presence in the population

		z[i]~dbern(psi)

		#### Stoat home range centres in space
		stoat.pos[i, 1]~dunif(min.x, max.x)

		stoat.pos[i, 2]~dunif(min.y, max.y)

		}

	##### Observation model - distances between home ranges and each trap and probability of capture

	for (i in 1:n.tot){

			for (j in 1:n.trap){
			
				### Distances between stoat home range centres and each trap
				d[i, j]<-sqrt(pow(stoat.pos[i, 1]-trap.coord[j, 1], 2)+pow(stoat.pos[i, 2]-trap.coord[j, 2],2))

				### Probability of capture
				p.capt[i, j]<-g0*exp(-d[i, j]*d[i, j]/(2*sigma*sigma))

				### Effective probability of capture (presence and p.capt)
				p.eff[i, j]<-p.capt[i, j]*z[i]

				comp.eff[i, j]<-1-p.eff[i, j]		#### Complement of the effective capture prob

				}
			}

	#### Probability of a capture in each trap

	for (t in 1:n.trap){
		
		p.trap[t]<-1-prod(comp.eff[1:n.tot, t])		### Effective probability of capturing at least one stoat in that trap

		stoat.cap[t]~dbern(p.trap[t])
	}

	N<-sum(z[1:n.tot])
})
	
      

## Initial values - where to start the MCMC chains


inits<-list(sigma=runif(1, 0.001, 1000), psi=rbeta(1, 1, 1), g0=rbeta(1, 1, 1), 
	z=c(rep(1, sum(stoat.capt)), rep(0, n.tot-sum(stoat.capt))))

### Define constant values
const.rem<-list(n.tot=n.tot, n.trap=n.trap, min.x=x.min, max.x=x.max, min.y=y.min, max.y=y.max)

#### Put the data together
data.rem<-list(stoat.cap=stoat.capt, trap.coord=trap.coord)


###### THE MODEL        
Rmodel<-nimbleModel(code=se.stoat, constants=const.rem, data=data.rem, inits=inits)

rem.conf<-configureMCMC(Rmodel,  monitors=list('sigma', 'psi', 'g0', 'N'), thin=1)

Rmcmc<-buildMCMC(rem.conf)

Cmodel<-compileNimble(Rmodel)

Cmcmc<-compileNimble(Rmcmc, project = Rmodel) 
 
### Run the model with three chain and check how it behaves
rem.m2<-runMCMC(Cmcmc, niter=70000, nchains=3, nburnin=10000, inits=inits, samplesAsCodaMCMC=FALSE)

out.post<-do.call(rbind.data.frame, rem.m2)

summary(out.post)

beep(sound=8)

#### Plot estimated vs actual
par(mfrow=c(1, 2))

hist(c(unlist(out.post[, 1])), ylab="Frequency", xlab="Sigma", main="Distance decay", xlim=c(0, 1000))

abline(v=sigma, lwd=5)

hist(c(unlist(out.post[, 2])), ylab="Frequency", xlab="g0", main="Pob. capt centre HR", xlim=c(0, 1))

abline(v=g0, lwd=5)

hist(c(unlist(out.post[, 2])), ylab="Frequency", xlab="N", main="Stoat population size", xlim=c(0, n.tot))

abline(v=n.stoat, lwd=5)
