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

rem.data<-read.table("d:/CONTAIN/Stoats RSPB/Removal-Tot.csv", sep=",", header=TRUE)

rem.data$Days<-round(rem.data$Mean_days, digits=0)

rem.data$Site2<-as.numeric(rem.data$Site)

rem.data

### Re-order the data frame by site

new.rem<-rem.data[order(rem.data$Site2), ]

new.rem

#### Effective trapping effort for stoats after accounting for captured rats and closed traps
eff.stoat<-(new.rem$Total.traps-new.rem$rat-new.rem$closed.trap)*new.rem$Days		 

eff.stoat

#### Number of sites
n.site<-max(new.rem$Site2)

n.site

#### Number of trapping occasions per site
n.occ<-tapply(new.rem$Session, new.rem$Site2, max)

n.occ

### Matrices (multiple site)
### Trapping effort per session and site
trap.eff<-matrix(NA, ncol=max(n.occ), nrow=n.site)

for (i in 1:n.site){

	trap.eff[i, 1:length(eff.stoat[new.rem$Site2==i])]<-eff.stoat[new.rem$Site2==i]

}

trap.eff

### Stoats captured
stoat.capt<-matrix(NA, ncol=max(n.occ), nrow=n.site)

for (i in 1:n.site){

	stoat.capt[i, 1:length(new.rem$stoat[new.rem$Site2==i])]<-new.rem$stoat[new.rem$Site2==i]

}

stoat.capt

#### NIMBLE BAYESIAN MODEL

rem.model<-nimbleCode({
 
    
        ## Prior for the probability of capture of stoats per trap and night
        alpha.stoat~dnorm(0, sd=3.16)

		beta.stoat~dnorm(0, sd=3.16)


    for (s in 1:n.site){

	 	## Prior for the stoat population size before any removal - whole of SR
        mean.stoat[s]~dunif(1, 1000)

        n.stoat[s]~dpois(mean.stoat[s])

		### Probability of not capturing a stoat after a certain trapping effort (trap nights)
		### And including rat by-catch effects	
		#### First occasion
		logit(p.stoat[s, 1])<-alpha.stoat+beta.stoat*eff.stoat[s, 1]
	
		stoat.cap[s, 1]<-p.stoat[s, 1]

        ### Likelihood - Stoats and rats removed per occasion
		### First occasion
		av.stoat[s, 1]<-n.stoat[s]		### Stoats not yet caught and available for capture in the new occasion

		y.stoat[s, 1]~dbin(stoat.cap[s, 1], av.stoat[s, 1])	#### Stoats
	
		#### Subsequent occasions
		
        for (j in 2:n.occ[s]){

			## Stoats removed per occasion
			### Stoats not yet caught and available for capture in the new occasion
			av.stoat[s, j]<-n.stoat[s]-sum(y.stoat[s, 1:(j-1)])
		
			### Stoats removed
			### Probability of not capturing a stoat after a certain trapping effort (trap nights)
			### And including rat by-catch effects	
	
			logit(p.stoat[s, j])<-alpha.stoat+beta.stoat*eff.stoat[s, j]
	
			stoat.cap[s, j]<-p.stoat[s, j]

			y.stoat[s, j]~dbin(stoat.cap[s, j], av.stoat[s, j])		
	
			
               }
    
    }
      
})


## Initial values - where to start the MCMC chains

inits<-list(alpha.stoat=rnorm(1, 0, 3.16), beta.stoat=rnorm(1, 0, 3.16), n.stoat=rpois(n.site, max(apply(stoat.capt, 1, function(x) sum(x, na.rm=TRUE)))*2))

### Define constant values
const.rem<-list(n.occ=n.occ, n.site=n.site)

#### Put the data together
data.rem<-list(eff.stoat=log10(trap.eff), y.stoat=stoat.capt)


###### THE MODEL        
Rmodel<-nimbleModel(code=rem.model, constants=const.rem, data=data.rem, inits=inits)

rem.conf<-configureMCMC(Rmodel,  monitors=list('alpha.stoat', 'beta.stoat', 'n.stoat'), thin=1)

Rmcmc<-buildMCMC(rem.conf)

Cmodel<-compileNimble(Rmodel)

Cmcmc<-compileNimble(Rmcmc, project = Rmodel) 
 
### Run the model with one chain
rem.m1<-runMCMC(Cmcmc, niter=100000, nchains=1, nburnin=1000, inits=inits, samplesAsCodaMCMC=TRUE)

summary(rem.m1)

beep(sound=8)

raftery.diag(rem.m1)

plot(rem.m1)

### Run the model with three chains
rem.m2<-runMCMC(Cmcmc, niter=500000, nchains=3, nburnin=2000, inits=inits, samplesAsCodaMCMC=TRUE)

summary(rem.m2)

beep(sound=8)

gelman.diag(rem.m2)

plot(rem.m2)
