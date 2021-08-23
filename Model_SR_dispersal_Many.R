### R script to run the model multiple times
### with sensitivity analysis

set.seed(200)

library(compiler)

#### This is just calling the simulation function.
source("C:/Orkney/Model/Model _SR_dispersal.R")

############# SIMULATIONS
### Define the study area (same across simulations)
SR <- readOGR("C:/Orkney/Model/SRonaldsay.shp")
tot.area<-sum(SR$SR_area)/10000
#Location of traps
SR_traps <- readOGR("C:/Orkney/Model/SR_traps.shp")


#resistanceSR <- resistanceFromShape("C:/Orkney/Model/SR_habitats.shp", res = 10, field = "Habitat", mapvalues = c("coast" = 0.1, "grassland" = 0.3, "meadow" = 0.8, "heatland" = 0.1,"marine" = 0.1, "road" = 0.1, "town" = 0.9,"water" = 1.00, "wood" = 0.75 ), background = 1.00, margin = 1000)
# Resistance of habitats for wlaking staots. Below resistance for all ghabitats is set to zero except for water and sea (set to 1 so not possible to cross)
resistanceSR <- resistanceFromShape("C:/Orkney/Model/SR_habitats.shp", res = 10, field = "Habitat", mapvalues = c("coast" = 0.0, "grassland" = 0.0, "meadow" = 0.0, "heatland" = 0.0,"marine" = 0.0, "road" = 0.0, "town" = 0.0,"water" = 1.00, "wood" = 0.0 ), background = 1.00, margin = 1000)



############# NUMERICAL SIMULATIONS
n.sim<-1000  ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 9)

CRW<-round(qunif(lhs1[,1], 0.95, 0.99), digits=4)               ### correlation between current and previous step, in this case highly correlated thus directional walk

step<-round(qunif(lhs1[,2], 100, 300), digits=2)                ### length of one step

horizon<-round(qunif(lhs1[,3], 100, 500), digits=2)             ### scanning horizon (distance when individual sees resistance of habotat on teh front of it)

time<-round(qunif(lhs1[,4], 200, 500), digits=2)                ### walking time (number times teh steps will be walked) )

n.day<-round(qunif(lhs1[,5], 14, 28), digits=0)                 ### number of trapping days

init.density<-round(qunif(lhs1[,6], 0.014, 0.019), digits=4)    ### initial density of stoats

g0<-round(qunif(lhs1[,7], 0.01, 0.1), digits=4)                 ### probability of capture at the center of the home range

sigma<-round(qunif(lhs1[,8], 50, 100), digits=0)                ### spatial decay parameter

p.bycatch<-round(qunif(lhs1[,9], 0.005, 0.015), digits=5)       ### Probability of by-catch per trap day


### Empty list to store the results
#stoat.sim<-vector("list")                                      ### list to store coordinates
stoat.sim<-matrix(0, nrow=2, ncol=n.sim)                        ### matrix to store values of the output 

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulations                                         ### should be defined as stoat.sim[[i]] if coordinates will be stored

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

mean(stoat.sim [1, ])     ### averaged percent of captured stoats
mean(stoat.sim [2, ])     ### averaged number of not captured staots


### Convert to data frame (for coordinates)
#export.coordinates<-do.call(rbind.data.frame, stoat.sim)

### Export (for coordinates)
#write.table(export.coordinates,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")


### Export (for percentage of captured stoats and number of not captured stoats)
write.table(stoat.sim,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")

### Sensitivity analysis
#### Now the boosted regression tree

library(dismo)
library(gbm)


### Put all the data together

sens.data <- data.frame(prop.stoat=stoat.sim[1, ], CRW = CRW, step = step, horizon=horizon, time = time,
                        
                        n.day=n.day,init.density=init.density,g0=g0, sigma = sigma, p.bycatch = p.bycatch)

sens.mod<-gbm.step(sens.data, gbm.x=c(2:10), gbm.y=1, tree.complexity=3, learning.rate=0.001, bag.fraction=0.7, step.size=100, max.trees=10000,  n.folds=10, family="gaussian")

summary(sens.mod)

plot(sens.mod,4)  ### plot for different parameters indicated by number after comma


