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

stoat.trapping$Trap.ID<-as.numeric(stoat.trapping$TRAP_ID)

#### Create a new variable indicating whether there was a stoat capture
stoat.trapping$stoat.capt<-ifelse(stoat.trapping$Catch_type=="S", 1, 0) 

stoat.trapping$stoat.capt

stoat.trapping$Catch_type

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


#### Establish the maximum number of stoats
area.sr<-49.8		#### Area (km2) of South Ronaldsay

n.tot<-round(10*area.sr, digits=0)			### Maximum stoat abundance (8 is a maximum density/km2 in NZ)

n.tot<-300

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

pos.init<-matrix(0, ncol=2, nrow=n.tot)

pos.init[, 1]<-runif(n.tot, x.min, x.max)

pos.init[, 2]<-runif(n.tot, y.min, y.max)

inits<-list(sigma=runif(1, 0.001, 1000), psi=rbeta(1, 1, 1), g0=rbeta(1, 1, 1), 
	z=c(rep(1, sum(stoat.trapping$stoat.capt)), rep(0, n.tot-sum(stoat.trapping$stoat.capt))))

### Define constant values
const.rem<-list(n.tot=n.tot, n.trap=n.trap, min.x=x.min, max.x=x.max, min.y=y.min, max.y=y.max)

#### Put the data together
data.rem<-list(stoat.cap=stoat.analys$stoat.capt, trap.coord=trap.coord)


###### THE MODEL        
Rmodel<-nimbleModel(code=se.stoat, constants=const.rem, data=data.rem, inits=inits)

rem.conf<-configureMCMC(Rmodel,  monitors=list('N'), thin=1)

Rmcmc<-buildMCMC(rem.conf)

Cmodel<-compileNimble(Rmodel)

Cmcmc<-compileNimble(Rmcmc, project = Rmodel) 
 
### Run the model with three chain and check how it behaves
rem.m2<-runMCMC(Cmcmc, niter=70000, nchains=3, nburnin=10000, inits=inits, samplesAsCodaMCMC=FALSE)

beep(sound=8)

### Summary of the posterior & exclude NAs
out.post<-do.call(rbind.data.frame, rem.m2)

hist(c(unlist(out.post)), ylab="Frequency", xlab="Stoat population size", main="Estimated stoat population size", xlim=c(0, 100))

abline(v=42, lwd=5)

### Posterior estimates of population size
### Mean
mean(c(unlist(out.post)))

### Standard deviation
sd(c(unlist(out.post)))

### 95% Credible Intervals
quantile(c(unlist(out.post)), c(0.0275, 0.975))


