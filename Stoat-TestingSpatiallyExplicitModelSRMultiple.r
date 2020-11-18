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
library(lhs)

########## Load the trapping data

stoat.trapping<-read.table("d:/CONTAIN/Stoats RSPB/Simulations Movement/Spatial data/SpatExp2/Trapping-data.csv", header=TRUE, sep=",")

summary(stoat.trapping)

### For this example, use first trapping session only

stoat.analys<-stoat.trapping[stoat.trapping$VISIT_NO==1, ]

summary(stoat.analys)


### Data-generation function to simulate data
### stoat data: capture data at the trap-level and including coordinate
### n.stoat: stoat abundance
### g0: probability of capture at the centre of the home range
### sigma: distance decay in the probability of capture
### p.bycath: probability of bycatch

data.gen<-function(stoat.data, n.stoat, g0, sigma, p.bycatch){

	#### Sample size, number of traps, and other for the Bayesian analysis

	n.trap<-nrow(stoat.data)		### Number of traps (considering two traps per actual trap)

	### Create a matrix to store the trap coordinates
	trap.coord<-matrix(0, ncol=2, nrow=n.trap)

	trap.coord[, 1]<-stoat.data$X

	trap.coord[, 2]<-stoat.data$Y

	x.min<-min(trap.coord[, 1])-1000

	x.max<-max(trap.coord[, 1])+1000

	y.min<-min(trap.coord[, 2])-1000

	y.max<-max(trap.coord[, 2])+1000

	##### Simulate stoats and capture data

	pos.init<-matrix(0, ncol=2, nrow=n.stoat)       #### Location of their home range centres (random uniform)

	pos.init[, 1]<-runif(n.stoat, x.min, x.max)

	pos.init[, 2]<-runif(n.stoat, y.min, y.max)

	### Transform to point patterns
	trap.loc<-ppp(trap.coord[, 1], trap.coord[, 2], check=FALSE)

	centre.stoat<-ppp(pos.init[, 1], pos.init[, 2], check=FALSE)

	###### Calculate distances between stoat home range centres and each trap
	dist.trap<-crossdist(centre.stoat, trap.loc)

	#### Simulate spatially-explicit stoat captures

	capt.trap<-matrix(0, ncol=n.stoat, nrow=n.trap)

		for (i in 1:n.stoat){

        	p.capt<-g0*exp((-dist.trap[i, ]*dist.trap[i, ])/(2*sigma^2))      ## Probability of capture per stoat per night per trap

        	capt.trap[, i]<-rbinom(n.trap, size=1, prob=p.capt*p.bycatch)

         }

		stoat.capt<-ifelse(rowSums(capt.trap)!=0, 1, 0)          #### Strike a capture if at least one has been caught in the trap

	##### Output data
	out.sim<-vector("list")

	out.sim$stoat.capt<-stoat.capt

	out.sim$g0<-g0

	out.sim$sigma<-sigma

	out.sim$n.stoat<-n.stoat

	out.sim$n.trap<-n.trap

	out.sim$trap.coord<-trap.coord

	out.sim$x.min<-x.min

	out.sim$x.max<-x.max

	out.sim$y.min<-y.min

	out.sim$y.max<-y.max

	return(out.sim)
}

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


	#### Errors compared to true values
	err.n<-true.n-N

	err.g0<-true.g0-g0

	err.sigma<-true.sigma-sigma

})
	
      
#### Function to fit the model to the simulated data

mod.fit<-function(sim.data, n.tot){

	## Initial values - where to start the MCMC chains
	inits<-list(sigma=runif(1, 0.001, 1000), psi=rbeta(1, 1, 1), g0=rbeta(1, 1, 1), 
		z=c(rep(1, sum(sim.data$stoat.capt)), rep(0, n.tot-sum(sim.data$stoat.capt))))

	### Define constant values
	const.rem<-list(n.tot=n.tot, n.trap=sim.data$n.trap, min.x=sim.data$x.min, max.x=sim.data$x.max, min.y=sim.data$y.min, max.y=sim.data$y.max)

	#### Put the data together
	data.rem<-list(stoat.cap=sim.data$stoat.capt, trap.coord=sim.data$trap.coord, true.n=sim.data$n.stoat, true.g0=sim.data$g0, true.sigma=sim.data$sigma)

	###### THE MODEL        
	Rmodel<-nimbleModel(code=se.stoat, constants=const.rem, data=data.rem, inits=inits)

	rem.conf<-configureMCMC(Rmodel,  monitors=list('err.sigma', 'err.g0', 'err.n'), thin=1)

	Rmcmc<-buildMCMC(rem.conf)

	Cmodel<-compileNimble(Rmodel)

	Cmcmc<-compileNimble(Rmcmc, project = Rmodel) 
 
	### Run the model with three chain and check how it behaves
	rem.m2<-runMCMC(Cmcmc, niter=70000, nchains=3, nburnin=10000, inits=inits, samplesAsCodaMCMC=TRUE)

	return(summary(rem.m2)$statistics[, 1])

}


### Latin hypercube design

n.sim<-100

lhs1<-randomLHS(n.sim, 4)

stoat.ab<-qpois(lhs1[,1], 50)          ### Stoat abundance

g0.sim<-round(qunif(lhs1[ ,2], min=0, max=1), digits=2)		#### Simulated g0

sigma.sim<-round(qunif(lhs1[,3], min=0.01, max=1000), digits=2)		### Simulated sigma

bycatch.sim<-round(qunif(lhs1[,4], min=0.001, max=1), digits=2)		### Simulated rat bycatch

### Check number of cores
n.core<-detectCores()

n.core

#### Parallel
registerDoParallel(n.core)

#### Timing the function
start.time<-Sys.time()

foreach(i=1:n.sim, .combine = rbind, .packages=c("nimble", "spatstat")) %dopar% mod.fit(data.gen(stoat.analys, stoat.ab[i], g0.sim[i], sigma.sim[i], bycatch.sim[i]), n.tot=100)

#### Stop cluster
stopImplicitCluster()

end.time<-Sys.time()

end.time-start.time     ### Time to run

beep(sound=8)
