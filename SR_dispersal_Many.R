### R script to run the model multiple times
### with sensitivity analysis

set.seed(200)

library(compiler)

#### This is just calling the simulation function.
source("C:/Orkney/Model/SR_dispersal.R")

############# SIMULATIONS
### Define the study area (same across simulations)
SR <- readOGR("C:/Orkney/Model/SRonaldsay.shp")
tot.area<-sum(SR$SR_area)/10000
SR_traps <- readOGR("C:/Orkney/Model/SR_traps_500in.shp")
SR_traps_out <- readOGR("C:/Orkney/Model/SR_traps_500out.shp")
resistanceSR <- resistanceFromShape("C:/Orkney/Model/SR_habitats.shp", res = 10, field = "Habitat", mapvalues = c("coast" = 0.1, "grassland" = 0.3, "meadow" = 0.8, "heatland" = 0.1,"marine" = 0.1, "road" = 0.1, "town" = 0.9,"water" = 1.00, "wood" = 0.75 ), background = 1.00, margin = 1000)
#resistanceSR <- resistanceFromShape("C:/Orkney/Model/SR_habitats.shp", res = 10, field = "Habitat", mapvalues = c("coast" = 0.0, "grassland" = 0.0, "meadow" = 0.0, "heatland" = 0.0,"marine" = 0.0, "road" = 0.0, "town" = 0.0,"water" = 1.00, "wood" = 0.0 ), background = 1.00, margin = 1000)
#tot.area<-sum(main$Area3)/10000                     ### Area covered in ha (Mainland)
#tot.area <- 52325



############# NUMERICAL SIMULATIONS
n.sim<-10  ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 9)

CRW<-round(qunif(lhs1[,1], 0.75, 0.99), digits=4)           ### Probability of by-catch per trap day

step<-round(qunif(lhs1[,2], 25, 50), digits=2)           ### Probability of by-catch per trap day

horizon<-round(qunif(lhs1[,3], 50, 100), digits=2)           ### Probability of by-catch per trap day

time<-round(qunif(lhs1[,4], 200, 500), digits=2)           ### Probability of by-catch per trap day

n.day<-round(qunif(lhs1[,5], 14, 28), digits=0)

init.density<-round(qunif(lhs1[,6], 0.005, 0.1), digits=4)

g0<-round(qunif(lhs1[,7], 0.01, 0.1), digits=4)           ### Probability of capture at the centre of the home range

sigma<-round(qunif(lhs1[,8], 150, 300), digits=0)           ### Spatial decay parameter

p.bycatch<-round(qunif(lhs1[,9], 0.005, 0.015), digits=5)           ### Probability of by-catch per trap day


### Empty list to store the results
#stoat.sim<-vector("list")
stoat.sim<-matrix(0, nrow=2, ncol=n.sim)

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulations

for (i in 1:n.sim){
  
    stoat.sim[ ,i]<-stoat.dist(CRW = CRW[i], step = step[i], horizon=horizon[i], time = time[i],
                               
                          n.day=n.day[i],init.density=init.density[i],  
                           
                           g0=g0[i], sigma = sigma[i], p.bycatch = p.bycatch[i])
  
    setTxtProgressBar(pb,i)
    }



close(pb)
end.time<-Sys.time()
time.taken<-end.time-start.time
time.taken

mean(stoat.sim [1, ])
mean(stoat.sim [2, ])


### Convert to data frame
#export.coordinates<-do.call(rbind.data.frame, stoat.sim)

### Export
#write.table(export.coordinates,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")

#write.table(HRcentersDF,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")





### Exploratory plotting
par(mfrow=c(2, 3))

plot(stoat.sim[1, ]~time, xlab="Initial density", ylab="Percentage caught")

plot(stoat.sim[1, ]~CRW, xlab="Average distance between traps (metres)", ylab="Percentage caught")

plot(stoat.sim[1, ]~g0, xlab="Probability of capture", ylab="Percentage caught")

plot(stoat.sim[1, ]~sigma, xlab="Spatial decay parameter", ylab="Percentage caught")

plot(stoat.sim[1, ]~p.bycatch, xlab="Probability of by-catch", ylab="Percentage caught")

plot(stoat.sim[1, ]~step, xlab="Number of traps", ylab="Percentage caught")

plot(stoat.sim[1, ]~horizon, xlab="Number of days", ylab="Percentage caught")


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

# Plot CRW vs. percentage caught

Perc_caught <- stoat.sim[1,]
CRW <- CRW

Perc_caught <- as.data.frame(Perc_caught)
CRW <- as.data.frame(CRW)
DataDF <- cbind(data.frame(CRW), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ CRW + I(CRW^2), data = DataDF))
prd <- data.frame(CRW = seq(from = range(DataDF$CRW)[1], to = range(DataDF$CRW)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = CRW, y = fit)) +
  labs(x="CRW", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = CRW, y = Perc_caught))


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
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = Spat_dec, y = Perc_caught))

# Plot probability of by-catch vs. percentage caught

Perc_caught <- stoat.sim[1,]
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

# Plot Number of steps vs. percentage caught

Perc_caught <- stoat.sim[1,]
time <- time

Perc_caught <- as.data.frame(Perc_caught)
Traps <- as.data.frame(time)
DataDF <- cbind(data.frame(time), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ time + I(time^2), data = DataDF))
prd <- data.frame(time = seq(from = range(DataDF$time)[1], to = range(DataDF$time)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = time, y = fit)) +
  labs(x="Number of steps", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = time, y = Perc_caught))

# Plot step length vs. percentage caught

Perc_caught <- stoat.sim[1,]
step <- step

Perc_caught <- as.data.frame(Perc_caught)
Days <- as.data.frame(step)
DataDF <- cbind(data.frame(step), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ step + I(step^2), data = DataDF))
prd <- data.frame(step = seq(from = range(DataDF$step)[1], to = range(DataDF$step)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = step, y = fit)) +
  labs(x="Step length", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = step, y = Perc_caught))

# Plot perceptual range vs. percentage caught

Perc_caught <- stoat.sim[1,]
horizon <- horizon

Perc_caught <- as.data.frame(Perc_caught)
Days <- as.data.frame(horizon)
DataDF <- cbind(data.frame(horizon), data.frame(Perc_caught))

fit <- (lm(Perc_caught ~ horizon + I(horizon^2), data = DataDF))
prd <- data.frame(horizon = seq(from = range(DataDF$horizon)[1], to = range(DataDF$horizon)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = horizon, y = fit)) +
  labs(x="Perceptual range", y="Percentage caught")+
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = DataDF, aes(x = horizon, y = Perc_caught))



### Sensitivity analysis

### First part taken from our existing scripts

#### Now the boosted regression tree

library(dismo)
library(gbm)


### Put all the data together

sens.data <- data.frame(post.capt=stoat.sim[1, ], CRW = CRW, step = step, horizon=horizon, time = time,n.day = n.day, init.density = init.density, g0=g0, sigma = sigma, p.bycatch = p.bycatch)

sens.mod<-gbm.step(sens.data, gbm.x=c(2:9), gbm.y=1, tree.complexity=3, learning.rate=0.0001, bag.fraction=0.7, step.size=100, max.trees=10000,  n.folds=10, family="gaussian")

summary(sens.mod)

plot(sens.mod)

plot(sens.mod, i.var = 1)
