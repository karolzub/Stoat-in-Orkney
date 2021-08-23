### R script to run the model multiple times
### with sensitivity analysis

set.seed(200)

library(compiler)

#### This is just calling the simulation function.
source("C:/Orkney/Frontiers/Model3_1.R")

stoat.hr<-cmpfun(stoat.dist)

#### Load the trap data (trap positions)
trap.pos<-readOGR("C:/Orkney/Model/Trap_locations_Orkney.shp")
n.trap<-length(trap.pos$ID)
trap.pos<-as.data.frame(trap.pos)
#plot(trap.pos, add = TRUE)

############# SIMULATIONS
### Define the study area (same across simulations)
orkney<-readOGR("C:/Orkney/Model/Orkney.shp")
tot.area<-sum(orkney$Area)/10000                     ### Area covered in ha (Orkney)
#plot(orkney)

############# NUMERICAL SIMULATIONS
n.sim<-1000     ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 5)

n.day<-round(qunif(lhs1[,1], 1, 28), digits=0)    ### number of trapping days in a session

init.density<-round(qunif(lhs1[,2], 0.005, 0.1), digits=4)    ### initial density of stoats per ha

g0<-round(qunif(lhs1[,3], 0.01, 0.1), digits=4)     ### Probability of capture at the centre of the home range

sigma<-round(qunif(lhs1[,4], 100, 700), digits=0)    ### Spatial decay parameter

p.bycatch<-round(qunif(lhs1[,5], 0.005, 0.015), digits=5)    ### Probability of by-catch per trap day

### Empty matrix to store the results
#stoat.sim<-matrix(0, nrow=2, ncol=n.sim)
stoat.sim<-vector("list")

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulations
for (i in 1:n.sim){
        
        stoat.sim[[i]]<-stoat.hr(n.day=n.day[i],init.density=init.density[i], 
                                   g0=g0[i], sigma = sigma[i], p.bycatch = p.bycatch[i])

        setTxtProgressBar(pb,i)
      
}

close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

### Convert to data frame
export.coordinates<-do.call(rbind.data.frame, stoat.sim)

### Export
write.table(export.coordinates,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")


## Export
#write.table(stoat.sim,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")


### Exploratory plotting
par(mfrow=c(2, 3))

plot(stoat.sim[1, ]~n.day, xlab="Number of days", ylab="Percentage caught")

plot(stoat.sim[1, ]~init.density, xlab="Initial density", ylab="Percentage caught")

plot(stoat.sim[1, ]~n.trap, xlab="Number of traps", ylab="Percentage caught")

plot(stoat.sim[1, ]~g0, xlab="Probability of capture", ylab="Percentage caught")

plot(stoat.sim[1, ]~sigma, xlab="Spatial decay parameter", ylab="Percentage caught")

plot(stoat.sim[1, ]~p.bycatch, xlab="Probability of a by-catch", ylab="Percentage caught")



### Plotting with fitted plynomial regression lines and 95% CI
library(ggplot2)

# Plot initial density vs. percentage caught

Perc_caught <- stoat.sim[1,]
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

Perc_caught <- stoat.sim[1,]
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
  #geom_line(color="red", size = 4) +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Aver_dist, y = Perc_caught))


# Plot probability of capture vs. percentage caught

Perc_caught <- stoat.sim[1,]
Prob_cap <- g0

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
  labs(x="Probability of capture per night", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Prob_cap, y = Perc_caught))

# Plot spatial decay parameter vs. percentage caught

Perc_caught <- stoat.sim[1,]
Spat_dec <- sigma

Perc_caught <- as.data.frame(Perc_caught)
Mov_corr <- as.data.frame(Spat_dec)
DataDF <- cbind(data.frame(Spat_dec), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Spat_dec + I(Spat_dec^2), data = DataDF))
prd <- data.frame(Spat_dec = seq(from = range(DataDF$Spat_dec)[1], to = range(DataDF$Spat_dec)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Spat_dec, y = fit)) +
  labs(x="Spatial decay parameter", y="Percentage caught")+
  theme_bw() +
  #geom_line(color="red", size = 4) +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Spat_dec, y = Perc_caught))

# Plot probability of by-catch vs. percentage caught

Perc_caught <- stoat.sim[1,]
Bycatch <- p.bycatch

Perc_caught <- as.data.frame(Perc_caught)
Mov_corr <- as.data.frame(p.bycatch)
DataDF <- cbind(data.frame(Bycatch), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Bycatch + I(Bycatch^2), data = DataDF))
prd <- data.frame(Bycatch = seq(from = range(DataDF$Bycatch)[1], to = range(DataDF$Bycatch)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Bycatch, y = fit)) +
  labs(x="Probability of a by-catch", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Bycatch, y = Perc_caught))

# Plot number of days vs. percentage caught

Perc_caught <- stoat.sim[1,]
Days <- n.day

Perc_caught <- as.data.frame(Perc_caught)
Days <- as.data.frame(Days)
DataDF <- cbind(data.frame(Days), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ Days + I(Days^2), data = DataDF))
prd <- data.frame(Days = seq(from = range(DataDF$Days)[1], to = range(DataDF$Days)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = Days, y = fit)) +
  labs(x="Number of days", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Days, y = Perc_caught))


### Sensitivity analysis

### First part taken from our existing scripts

### Latin hypercube desing to explore a wide variety & good coverage of values
n.sim<-1000

library(lhs)

lhs1<-randomLHS(n.sim, 7)

n.day<-round(qunif(lhs1[,1], 1, 28), digits=0)

init.density<-round(qunif(lhs1[,2], 0.005, 0.1), digits=4)

trap.dist<-round(qunif(lhs1[,4], 50, 250), digits=0)           ### Minimum distance between traps

g0<-round(qunif(lhs1[,5], 0.01, 0.1), digits=4)           ### Probability of capture at the centre of the home range

sigma<-round(qunif(lhs1[,6], 100, 700), digits=0)           ### Spatial decay parameter

p.bycatch<-round(qunif(lhs1[,6], 0.005, 0.015), digits=5)           ### Probability of by-catch per trap day


### Empty matrix to store the results
stoat.sim<-matrix(0, nrow=1, ncol=n.sim)

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)


### Run the simulations
for (i in 1:n.sim){
  
  stoat.sim[ , i]<-stoat.hr(n.day=n.day[i],init.density=init.density[i], trap.dist=trap.dist[i], 
                            g0=g0[i], sigma = sigma[i], p.bycatch = p.bycatch[i])
  
  setTxtProgressBar(pb,i)
  
}

close(pb)


#### Now the boosted regression tree

library(dismo)
library(gbm)


### Put all the data together

sens.data <- data.frame(prop.stoat=stoat.sim[1, ], n.day = n.day, init.density = init.density, g0=g0, sigma = sigma, p.bycatch=p.bycatch)

sens.mod<-gbm.step(sens.data, gbm.x=c(2:6), gbm.y=1, tree.complexity=3, learning.rate=0.001, bag.fraction=0.7, step.size=100, max.trees=10000,  n.folds=10, family="gaussian")

summary(sens.mod)

plot(sens.mod,4)  ### plot for different parameters indicated by number after comma

