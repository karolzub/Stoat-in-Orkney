### Stoat and rat data - Removal population & trappability estimates for South Ronaldsay
### Evaluating the convergence and mixing of the chains
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

########## Load the data

rem.data<-read.table("C:/Orkney/Model/Removal.csv", sep=";", header=TRUE)

rem.data

#### Effective trapping effort for stoats after accounting for captured rats and closed traps
eff.stoat<-(rem.data$Effort-rem.data$rat-rem.data$closed.trap)		 

eff.stoat

#### Effective trapping effort for rat after accounting for captured stoats and closed traps
eff.rat<-(rem.data$Effort-rem.data$stoat-rem.data$closed.trap)		 

eff.rat

#### Number of trapping occasions
n.occ<-max(rem.data$Session)

n.occ

#### NIMBLE BAYESIAN MODEL

rem.model<-nimbleCode({
 
	 ## Prior for the stoat population size before any removal - whole of SR
        mean.stoat~dunif(1, 500)

        n.stoat~dpois(mean.stoat)

        ## Prior for the rat population size before any removal - whole of SR
        mean.rat~dunif(100, 10000)

        n.rat~dpois(mean.rat)


        ## Prior for the probability of capture of stoats per trap and night
        alpha.stoat~dnorm(0, sd=3.16)

	beta.stoat~dnorm(0, sd=3.16)

	## Prior for the probability of capture of rats per trap and night
        alpha.rat~dnorm(0, sd=3.16)

	beta.rat~dnorm(0, sd=3.16)

	### Probability of not capturing a stoat after a certain trapping effort (trap nights)
	### And including rat by-catch effects	
	#### First occasion
	logit(p.stoat[1])<-alpha.stoat+beta.stoat*eff.stoat[1]
	
	stoat.cap[1]<-p.stoat[1]

	### Probability of not capturing a rat after a certain trapping effort (trap nights)
	### And including stoat by-catch effects	

	logit(p.rat[1])<-alpha.rat+beta.rat*eff.rat[1]
	
	rat.cap[1]<-p.rat[1]

        # Likelihood - Stoats and rats removed per occasion
	### First occasion
	av.stoat[1]<-n.stoat		### Stoats not yet caught and available for capture in the new occasion

	y.stoat[1]~dbin(stoat.cap[1], av.stoat[1])	#### Stoats

	av.rat[1]<-n.rat		### Rats not yet caught and available for capture in the new occasion

	y.rat[1]~dbin(rat.cap[1], av.rat[1])		#### Rats
	
	#### Subsequent occasions
		
        for (j in 2:n.occ){

		## Stoats removed per occasion
		### Stoats not yet caught and available for capture in the new occasion
		av.stoat[j]<-n.stoat-sum(y.stoat[1:(j-1)])
		
		### Stoats removed
		### Probability of not capturing a stoat after a certain trapping effort (trap nights)
		### And including rat by-catch effects	
	
		logit(p.stoat[j])<-alpha.stoat+beta.stoat*eff.stoat[j]
	
		stoat.cap[j]<-p.stoat[j]

		y.stoat[j]~dbin(stoat.cap[j], av.stoat[j])		
	
	
		## Rats removed per occasion
        	### Rats not yet caught and available for capture in the new occasion
		av.rat[j]<-n.rat-sum(y.rat[1:(j-1)])

		### Probability of not capturing a rat after a certain trapping effort (trap nights)
		### And including stoat by-catch effects	

		logit(p.rat[j])<-alpha.rat+beta.rat*eff.rat[j]
	
		rat.cap[j]<-p.rat[j]
	
		y.rat[j]~dbin(rat.cap[j], av.rat[j])		
	

               }
      
})


## Initial values - where to start the MCMC chains

inits<-list(alpha.stoat=rnorm(1, 0, 3.16), beta.stoat=rnorm(1, 0, 3.16), alpha.rat=rnorm(1, 0, 3.16),
		beta.rat=rnorm(1, 0, 3.16), n.stoat=rpois(1, sum(rem.data$stoat)*2), n.rat=rpois(1, sum(rem.data$rat)*2))

### Define constant values
const.rem<-list(n.occ=n.occ)

#### Put the data together
data.rem<-list(eff.stoat=log10(eff.stoat), eff.rat=log10(eff.rat), y.stoat=rem.data$stoat, y.rat=rem.data$rat)


###### THE MODEL        
Rmodel<-nimbleModel(code=rem.model, constants=const.rem, data=data.rem, inits=inits)

rem.conf<-configureMCMC(Rmodel,  monitors=list('alpha.stoat', 'beta.stoat', 'alpha.rat', 'beta.rat', 'n.stoat', 'n.rat'), thin=1)

Rmcmc<-buildMCMC(rem.conf)

Cmodel<-compileNimble(Rmodel)

Cmcmc<-compileNimble(Rmcmc, project = Rmodel) 
 
### Run the model with one chain
rem.m1<-runMCMC(Cmcmc, niter=1000000, nchains=1, nburnin=500000, inits=inits, samplesAsCodaMCMC=TRUE)

summary(rem.m1)

beep(sound=8)

raftery.diag(rem.m1)

plot(rem.m1)

### Run the model with three chains
rem.m2<-runMCMC(Cmcmc, niter=1000000, nchains=3, nburnin=500000, inits=inits, samplesAsCodaMCMC=TRUE)

summary(rem.m2)

beep(sound=8)

gelman.diag(rem.m2)

plot(rem.m2)
