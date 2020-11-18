### R script to run the model multiple times

set.seed(200)

library(compiler)

#### This is just calling the simulation function.
source("C:/Orkney/Model/StoatSimHR_SouthRonaldsay.R")

stoat.hr<-cmpfun(stoat.dist)

#### Load the trap data (trap positions)
trap.pos<-readOGR("C:/Orkney/Model/TrapLocations.shp")
n.trap<-length(trap.pos$ID)
trap.pos<-as.data.frame(trap.pos)

############# SIMULATIONS
### Define the study area (same across simulations)
#sr<-readOGR("C:/Orkney/Model/SR_area.shp")
#tot.area<-sum(sr$Shape_Area)/10000                     ### Area covered in ha (South Ronaldsay)

############# NUMERICAL SIMULATIONS
n.sim<-1000     ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 7)

n.day<-round(qunif(lhs1[,1], 7, 40), digits=0)
init.density<-round(qunif(lhs1[,2], 0.01, 0.1), digits=3)
n.trap<-round(qunif(lhs1[,3], 100, 500), digits=0)           ### Total number of traps
trap.dist<-round(qunif(lhs1[,4], 50, 250), digits=0)           ### Minimum distance between traps
cor.mov<-round(qunif(lhs1[,5], 0.001, 0.8), digits=3)           ### Position correlation between days across stoats
max.det<-round(qunif(lhs1[,6], 0.11, 0.8), digits=3)           ### Maximum probability of capture per night (across stoats)
p.bycatch<-round(qunif(lhs1[,7], 0.005, 0.015), digits=5) 

### Empty matrix to store the results
stoat.sim<-matrix(0, nrow=3, ncol=n.sim)    #### The first column is the %caught & the third are the remaining stoats

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)


### Run the simulations
for (i in 1:n.sim){
        
        stoat.sim[ , i]<-stoat.hr(n.day=n.da[i], init.density=init.density[i], n.trap = n.trap[i], trap.dist=trap.dist[i], 
                                   cor.mov=cor.mov[i], max.det=max.det[i], p.bycatch=p.bycatch[i])

        setTxtProgressBar(pb,i)
      
}

close(pb)

library(pse)

factors <- c("n.day", "init.density", "n.trap", "trap.dist","cor.mov", "max.det", "p.bycatch")
q <- c("qunif","qunif", "qunif", "qunif","qunif", "qunif", "qunif")
q.arg <- list(list(min=7, max=40),list(min=0.01, max=0.1), list(min=250, max=1000),list(min=50, max=250),
              list(min=0.001, max=0.8), list(min=0.11, max=0.8),
              list(min=0.005, max=0.015) )


modelRun <- function (init.density, n.trap, trap.dist, cor.mov, max.det, p.bycatch) {
  return(mapply(stoat.dist, lhs1[,1], lhs1[,2], lhs1[,3], lhs1[,4], lhs1[,5], lhs1[,6], lhs1[,7]))
}

myLHS <- LHS(modelRun, factors, 1000, q, q.arg, nboot = 100, repetitions = 10)

myLHS

get.data(myLHS)

get.results(myLHS)

plotecdf(myLHS)

plotprcc(myLHS)


### First part taken from our existing scripts

### Latin hypercube desing to explore a wide variety & good coverage of values

library(lhs)
n.sim <- 10000
lhs1<-randomLHS(n.sim, 7)
n.day<-round(qunif(lhs1[,1], 7, 40), digits=0)
init.density<-round(qunif(lhs1[,2], 0.01, 0.1), digits=3)
n.trap<-round(qunif(lhs1[,3], 100, 500), digits=0)           ### Total number of traps
trap.dist<-round(qunif(lhs1[,4], 50, 250), digits=0)           ### Minimum distance between traps
cor.mov<-round(qunif(lhs1[,5], 0.001, 0.8), digits=3)           ### Position correlation between days across stoats
max.det<-round(qunif(lhs1[,6], 0.11, 0.8), digits=3)           ### Maximum probability of capture per night (across stoats)
p.bycatch<-round(qunif(lhs1[,7], 0.005, 0.015), digits=5) 

### Empty matrix to store the results

stoat.sim<-matrix(0, nrow=3, ncol=n.sim)  #### The first column is the %caught & the third are the remaining stoats

#### Timing the function

start.time<-Sys.time()

### Progress bar

pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulations

for (i in 1:n.sim){
  
  stoat.sim[ , i]<-stoat.hr(n.day=n.day[i], init.density=init.density[i], n.trap = n.trap[i], trap.dist=trap.dist[i], 
                            cor.mov=cor.mov[i], max.det=max.det[i], p.bycatch=p.bycatch[i])
  
  setTxtProgressBar(pb,i)
  
}     

close(pb)

#### Now the boosted regression tree

library(dismo)
library(gbm)


### Put all the data together

sens.data <- data.frame(prop.stoat=stoat.sim[2, ], prop.stoat=stoat.sim[1, ],n.day = n.day, init.density = init.density, n.trap = n.trap, trap.dist= trap.dist, cor.mov=cor.mov, max.det = max.det, p.bycatch = p.bycatch)

sens.mod<-gbm.step(sens.data, gbm.x=c(2:9), gbm.y=1, tree.complexity=3, learning.rate=0.0001, bag.fraction=0.7, step.size=100, max.trees=10000,  n.folds=10, family="gaussian")

summary(sens.mod)

plot(sens.mod)




end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

### Exploratory plotting
par(mfrow=c(2, 3))

plot(stoat.sim[2, ]~init.density, xlab="Initial density", ylab="Percentage caught")

plot(stoat.sim[2, ]~trap.dist, xlab="Average distance between traps (metres)", ylab="Percentage caught")

plot(stoat.sim[2, ]~cor.mov, xlab="Stoat movement correlation", ylab="Percentage caught")

#plot(stoat.sim[2, ]~max.det, xlab="Maximum probability of capture per night", ylab="Percentage caught")

plot(stoat.sim[2, ]~n.trap, xlab="Number of traps", ylab="Percentage caught")

plot(stoat.sim[2, ]~p.bycatch, xlab="Probability of by-catch per night", ylab="Percentage caught")

plot(stoat.sim[2, ]~stoat.sim[1,], xlab="Average stoat home range (ha)", ylab="Percentage caught")


### Plotting with fitted plynomial regression lines and 95% CI
library(ggplot2)

# Plot initial density vs. percentage caught

Perc_caught <- stoat.sim[2,]
Init_dens <- init.density

Perc_caught <- as.data.frame(Perc_caught)
Init_dens <- as.data.frame(Init_dens)
DataDF <- cbind(data.frame(Init_dens), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Init_dens + I(Init_dens^2), data = DataDF))
prd <- data.frame(Init_dens = seq(from = range(DataDF$Init_dens)[1], to = range(DataDF$Init_dens)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Init_dens, y = fit)) +
  labs(x="Initial density", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Init_dens, y = Perc_caught))

# Plot average distance vs. percentage caught

Perc_caught <- stoat.sim[2,]
Aver_dist <- trap.dist

Perc_caught <- as.data.frame(Perc_caught)
Aver_dist <- as.data.frame(Aver_dist)
DataDF <- cbind(data.frame(Aver_dist), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Aver_dist + I(Aver_dist^2), data = DataDF))
prd <- data.frame(Aver_dist = seq(from = range(DataDF$Aver_dist)[1], to = range(DataDF$Aver_dist)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Aver_dist, y = fit)) +
  labs(x="Average distance between traps (metres)", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Aver_dist, y = Perc_caught))

# Plot movement correlation vs. percentage caught

Perc_caught <- stoat.sim[2,]
Mov_corr <- cor.mov

Perc_caught <- as.data.frame(Perc_caught)
Mov_corr <- as.data.frame(Mov_corr)
DataDF <- cbind(data.frame(Mov_corr), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Mov_corr + I(Mov_corr^2), data = DataDF))
prd <- data.frame(Mov_corr = seq(from = range(DataDF$Mov_corr)[1], to = range(DataDF$Mov_corr)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Mov_corr, y = fit)) +
  labs(x="Stoat movement correlation", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Mov_corr, y = Perc_caught))

# Plot max probability of capture vs. percentage caught

Perc_caught <- stoat.sim[2,]
Prob_cap <- max.det

Perc_caught <- as.data.frame(Perc_caught)
Prob_cap <- as.data.frame(Prob_cap)
DataDF <- cbind(data.frame(Prob_cap), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Prob_cap + I(Prob_cap^2), data = DataDF))
prd <- data.frame(Prob_cap = seq(from = range(DataDF$Prob_cap)[1], to = range(DataDF$Prob_cap)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Prob_cap, y = fit)) +
  labs(x="Maximum probability of capture per night", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Prob_cap, y = Perc_caught))

# Plot average home range vs. percentage caught

Perc_caught <- stoat.sim[2,]
hr <- stoat.sim[1,]

Perc_caught <- as.data.frame(Perc_caught)
hr <- as.data.frame(hr)
DataDF <- cbind(data.frame(hr), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ hr + I(hr^2), data = DataDF))
prd <- data.frame(hr = seq(from = range(DataDF$hr)[1], to = range(DataDF$hr)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = hr, y = fit)) +
  labs(x="Average stoat home range (ha)", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = hr, y = Perc_caught))

# Plot proportion of by-catch vs. percentage caught

Perc_caught <- stoat.sim[2,]
Bycatch <- p.bycatch

Perc_caught <- as.data.frame(Perc_caught)
Bycatch <- as.data.frame(Bycatch)
DataDF <- cbind(data.frame(Bycatch), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Bycatch + I(Bycatch^2), data = DataDF))
prd <- data.frame(Bycatch = seq(from = range(DataDF$Bycatch)[1], to = range(DataDF$Bycatch)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Bycatch, y = fit)) +
  labs(x="Probability of a by-catch per day", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Bycatch, y = Perc_caught))

# Plot number of traps vs. percentage caught

Perc_caught <- stoat.sim[2,]
Traps <- n.trap

Perc_caught <- as.data.frame(Perc_caught)
Traps <- as.data.frame(Traps)
DataDF <- cbind(data.frame(Traps), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Traps + I(Traps^2), data = DataDF))
prd <- data.frame(Traps = seq(from = range(DataDF$Traps)[1], to = range(DataDF$Traps)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Traps, y = fit)) +
  labs(x="Number of traps", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Traps, y = Perc_caught))


