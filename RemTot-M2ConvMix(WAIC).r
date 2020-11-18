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

rem.data<-read.table("d:/CONTAIN/Stoats RSPB/Removal-Updated3.csv", sep=",", header=TRUE)

rem.data$Days<-round(rem.data$Mean_days, digits=0)

rem.data$Site2<-as.numeric(rem.data$Site)

rem.data

### Re-order the data frame by site

new.rem<-rem.data[order(rem.data$Site2), ]

new.rem

#### Effective trapping effort for stoats after accounting for captured rats and closed traps
eff.stoat<-(new.rem$Total.traps-new.rem$closed.trap)		 

eff.stoat

#### Number of sites
n.site<-max(new.rem$Site2)

n.site

#### Number of trapping occasions per site
n.occ<-tapply(new.rem$Session, new.rem$Site2, max)

n.occ

### Matrices (multiple site)
### Trapping effort per night
trap1<-matrix(NA, ncol=max(n.occ), nrow=n.site)

for (i in 1:n.site){

	trap1[i, 1:length(eff.stoat[new.rem$Site2==i])]<-eff.stoat[new.rem$Site2==i]

}

trap1

### Number of nights
trap2<-matrix(NA, ncol=max(n.occ), nrow=n.site)

for (i in 1:n.site){

	trap2[i, 1:length(new.rem$Days[new.rem$Site2==i])]<-new.rem$Days[new.rem$Site2==i]

}

trap2

### Stoats captured
stoat.capt<-matrix(NA, ncol=max(n.occ), nrow=n.site)

for (i in 1:n.site){

	stoat.capt[i, 1:length(new.rem$stoat[new.rem$Site2==i])]<-new.rem$stoat[new.rem$Site2==i]

}

stoat.capt

#### NIMBLE BAYESIAN MODEL

rem.model<-nimbleCode({
 
    	### Prior for the abundance

		alpha.pop~dnorm(0, sd=3.16)

		beta.pop~dnorm(0, sd=3.16)

        ## Prior for the probability of capture of stoats per trap and night
        p.mean~dbeta(1, 1)

		c1~dunif(-1, 1)

		c2~dunif(-1, 1)

    for (s in 1:n.site){

	 	## Prior for the stoat population size before any removal
        log(mean.stoat[s])<-alpha.pop+beta.pop*area[s]

        n.stoat[s]~dpois(mean.stoat[s])

		stoat.const[s]~dconstraint(n.stoat[s]<max.stoat[s])

		### Probability of not capturing a stoat after a certain trapping effort (trap nights)
		### And including rat by-catch effects	
		#### First occasion
		p.stoat1[s, 1]<-1-pow(1-p.mean, c1*trap.eff1[s, 1])

		p.stoat[s, 1]<-1-pow(1-p.stoat1[s, 1], c2*trap.eff2[s, 1])
	
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
	
			p.stoat1[s, j]<-1-pow(1-p.mean, c1*trap.eff1[s, j])

			p.stoat[s, j]<-1-pow(1-p.stoat1[s, j], c2*trap.eff2[s, j])
	
			stoat.cap[s, j]<-p.stoat[s, j]

			y.stoat[s, j]~dbin(stoat.cap[s, j], av.stoat[s, j])		
	
			
               }
    
    }
      
})


## Initial values - where to start the MCMC chains
init.stoat<-apply(stoat.capt, 1, function(x) sum(x, na.rm=TRUE))

init.stoat

inits<-list(p.mean=rbeta(1, 1, 1), c1=runif(1,0, 1), c2=runif(1, 0, 1), alpha.pop=rnorm(1, 0, 3.16), beta.pop=rnorm(1, 0, 3.16),
		n.stoat=init.stoat+1)

### Define constant values
const.rem<-list(n.occ=n.occ, n.site=n.site)

#### Put the data together
data.rem<-list(trap.eff1=log10(trap1), trap.eff2=log10(trap2), y.stoat=stoat.capt, area=log10(unique(new.rem$Area)),
		max.stoat=round(20*unique(new.rem$Area), digits=0), stoat.const=rep(1, n.site))

###### THE MODEL        
Rmodel<-nimbleModel(code=rem.model, constants=const.rem, data=data.rem, inits=inits)

rem.conf<-configureMCMC(Rmodel,  monitors=list('p.mean', 'n.stoat', 'c1', 'c2', 'alpha.pop', 'beta.pop'), thin=1, enableWAIC = TRUE)

Rmcmc<-buildMCMC(rem.conf)

Cmodel<-compileNimble(Rmodel)

Cmcmc<-compileNimble(Rmcmc, project = Rmodel) 
 

### Run the model with three chains
rem.m2<-runMCMC(Cmcmc, niter=500000, nchains=3, nburnin=10000, inits=inits, samplesAsCodaMCMC=TRUE, WAIC=TRUE)

summary(rem.m2)

rem.m2$WAIC
