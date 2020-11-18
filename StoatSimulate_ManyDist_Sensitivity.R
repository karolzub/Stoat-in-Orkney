### R script to run the model multiple times

set.seed(200)

library(compiler)

#### This is just calling the simulation function.
source("C:/Orkney/Model/StoatSim_Dist_Sensiviti.R")

stoat.hr<-cmpfun(stoat.dist)

#### Load the trap data (trap positions)
trap.pos<-readOGR("C:/Orkney/Model/TrapLocations.shp")
n.trap<-length(trap.pos$ID)
trap.pos<-as.data.frame(trap.pos)

############# SIMULATIONS
### Define the study area (same across simulations)
sr<-readOGR("C:/Orkney/Model/SR_area.shp")
tot.area<-sum(sr$Shape_Area)/10000                     ### Area covered in ha (South Ronaldsay)

############# NUMERICAL SIMULATIONS
n.sim<-1000     ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 6)

init.density<-round(qunif(lhs1[,1], 0.01, 0.1), digits=4)

n.trap<-round(qunif(lhs1[,2], 250, 1000), digits=0)           ### Total number of traps

trap.dist<-round(qunif(lhs1[,3], 50, 250), digits=0)           ### Minimum distance between traps

g0<-round(qunif(lhs1[,4], 0.01, 0.1), digits=4)           ### Probability of capture at the centre of the home range

sigma<-round(qunif(lhs1[,5], 200, 700), digits=0)           ### Spatial decay parameter

p.bycatch<-round(qunif(lhs1[,6], 0.005, 0.015), digits=5)           ### Probability of by-catch per trap day


### Empty matrix to store the results
stoat.sim<-matrix(0, nrow=1, ncol=n.sim)

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)


### Run the simulations
for (i in 1:n.sim){
  
  stoat.sim[ , i]<-stoat.hr(n.day=20, init.density=init.density[i], n.trap = n.trap[i], trap.dist=trap.dist[i], 
                             g0=g0[i], sigma = sigma[i], p.bycatch = p.bycatch[i])
  
  setTxtProgressBar(pb,i)
  
}

close(pb)     

library(pse)

factors <- c("init.density", "n.trap", "trap.dist","g0", "sigma", "p.bycatch")
q <- c("qunif", "qunif", "qunif","qunif", "qunif", "qunif")
q.arg <- list(list(min=0.01, max=0.1), list(min=250, max=1000),list(min=50, max=250),
              list(min=0.01, max=0.1), list(min=200, max=700),
               list(min=0.005, max=0.015) )


modelRun <- function (init.density, n.trap, trap.dist, g0, sigma, p.bycatch) {
  return(mapply(stoat.dist, lhs1[,1], lhs1[,2], lhs1[,3], lhs1[,4], lhs1[,5], lhs1[,6]))
}

myLHS <- LHS(modelRun, factors, 1000, q, q.arg, nboot = 1000, repetitions = 10)

myLHS

get.data(myLHS)

get.results(myLHS)

plotecdf(myLHS)

corPlot(myLHS)

plotprcc(myLHS)

p <- pic(myLHS)

print(p$pic)


uncoupledLHS <- LHS(model=NULL, factors, 200, q, q.arg,nboot = 100, repetitions = 10)
write.csv(get.data(uncoupledLHS), file="C:/Orkney/Model/mydata.csv")

library(data.table)
stoat.sim <- transform(stoat.sim)

coupledLHS <- tell(uncoupledLHS, stoat.sim)

coupledLHS

plotprcc(coupledLHS)

### First part taken from our existing scripts

### Latin hypercube desing to explore a wide variety & good coverage of values

library(lhs)

lhs1<-randomLHS(n.sim, 6)
init.density<-round(qunif(lhs1[,1], 0.01, 0.1), digits=4)
n.trap<-round(qunif(lhs1[,2], 250, 1000), digits=0)           ### Total number of traps
trap.dist<-round(qunif(lhs1[,3], 50, 250), digits=0)           ### Minimum distance between traps
g0<-round(qunif(lhs1[,4], 0.01, 0.1), digits=4)           ### Probability of capture at the centre of the home range
sigma<-round(qunif(lhs1[,5], 200, 700), digits=0)           ### Spatial decay parameter
p.bycatch<-round(qunif(lhs1[,6], 0.005, 0.015), digits=5)           ### Probability of by-catch per trap day

### Empty matrix to store the results

stoat.sim<-matrix(0, nrow=1, ncol=n.sim)

#### Timing the function

start.time<-Sys.time()

### Progress bar

pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulations

for (i in 1:n.sim){
  
  stoat.sim[ , i]<-stoat.hr(n.day=20, init.density = init.density[i], n.trap = n.trap[i], trap.dist = trap.dist[i],
                            g0=g0[i], sigma = sigma[i], p.bycatch = p.bycatch[i])
  
  setTxtProgressBar(pb,i)
  
  
}      

close(pb)

#### Now the boosted regression tree

library(dismo)
library(gbm)


### Put all the data together

sens.data <- data.frame(prop.stoat=stoat.sim[1, ], init.density = init.density, n.trap = n.trap, trap.dist= trap.dist, g0=g0, sigma = sigma, p.bycatch = p.bycatch)

sens.mod<-gbm.step(sens.data, gbm.x=c(2:7), gbm.y=1, tree.complexity=3, learning.rate=0.001, bag.fraction=0.7, step.size=10, max.trees=10000,  n.folds=10, family="gaussian")

summary(sens.mod)

plot(sens.mod)















   
close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

### Exploratory plotting
par(mfrow=c(2, 3))

plot(stoat.sim[1, ]~init.density, xlab="Initial density", ylab="Percentage caught")

plot(stoat.sim[1, ]~trap.dist, xlab="Average distance between traps (metres)", ylab="Percentage caught")

plot(stoat.sim[1, ]~g0, xlab="Probability of capture", ylab="Percentage caught")

plot(stoat.sim[1, ]~sigma, xlab="Spatial decay parameter", ylab="Percentage caught")

plot(stoat.sim[1, ]~p.bycatch, xlab="Probability of by-catch", ylab="Percentage caught")

plot(stoat.sim[1, ]~n.trap, xlab="Number of traps", ylab="Percentage caught")


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


# Plot probability of capture vs. percentage caught

Perc_caught <- stoat.sim[2,]
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

Perc_caught <- stoat.sim[2,]
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
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Spat_dec, y = Perc_caught))

# Plot probability of by-catch vs. percentage caught

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
  labs(x="Probability of by-catch per night", y="Percentage caught")+
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

###Alternative way of plotting (without 95% CI)
newdat = data.frame(Bycatch = seq(min(DataDF$Bycatch), max(DataDF$Bycatch), length.out = 100))
newdat$pred = predict(fit, newdata = newdat)

plot(Perc_caught ~ Bycatch, xlab="Probability of by-catch per night", ylab="Percentage caught", data = DataDF)
with(newdat, lines(x = Bycatch, y = pred, col = "blue"))

