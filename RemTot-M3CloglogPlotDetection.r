### Stoat and rat data - Removal population & trappability estimates for South Ronaldsay
### Evaluating the convergence and mixing of the chains
set.seed(2012)


########## Load the data

rem.data<-read.table("d:/CONTAIN/Stoats RSPB/Removal-Updated2.csv", sep=",", header=TRUE)

rem.data$Days<-round(rem.data$Mean_days, digits=0)

rem.data$Site2<-as.numeric(rem.data$Site)

rem.data

### Re-order the data frame by site

new.rem<-rem.data[order(rem.data$Site2), ]

new.rem

#### Effective trapping effort for stoats after accounting for captured rats and closed traps
eff.stoat<-sort((new.rem$Total.traps-new.rem$closed.trap)*new.rem$Days)		 

eff.stoat

##### Simulate
n.sim<-20000

### Probability of capture
pred.prob<-matrix(0, ncol=n.sim, nrow=length(eff.stoat))

for (i in 1:n.sim){

    alpha<-rnorm(1, -5.18, 0.51)

    beta<-rnorm(1, 0.60, 0.16)

    y<-alpha+beta*(log10(eff.stoat))

    pred.prob[ , i]<-1-exp(-exp(y))

}


### Summary of the posterior
av.det<-rowMeans(pred.prob)

quant.det<-apply(pred.prob, 1, function(x) quantile(x, c(0.0275, 0.975)))

quant.det

plot(av.det~log10(eff.stoat), type="l", lwd=4, ylim=c(0, 0.4), main="Stoat probability of capture", ylab="Probability of capture",
		xlab="Trapping effort (log10 transformed)")

lines(log10(eff.stoat), quant.det[1, ], lwd=2, lty=2)

lines(log10(eff.stoat), quant.det[2, ], lwd=2, lty=2)


##### Relationship area and population size
area=sort(log10(unique(new.rem$Area)))

pred.size<-matrix(0, ncol=n.sim, nrow=length(area))

for (i in 1:n.sim){

    alpha<-rnorm(1, 2.55, 0.19)

    beta<-rnorm(1, 1.20, 0.15)

    pred.size1<-exp(alpha+beta*area)

    pred.size[, i]<-rpois(length(y), pred.size[,i])
}


### Summary of the posterior & exclude NAs
av.det<-rowMeans(pred.size)

quant.det<-apply(pred.size, 1, function(x) quantile(x, c(0.0275, 0.975)))

quant.det

plot(av.det~area, type="l", lwd=4,  main="Stoat abundance", ylab="Estimated abundance",
		xlab="Estimated effective trapping area (log10 transformed)", ylim=c(0, max(quant.det[2,])))

lines(area, quant.det[1, ], lwd=2, lty=2)

lines(area, quant.det[2, ], lwd=2, lty=2)
