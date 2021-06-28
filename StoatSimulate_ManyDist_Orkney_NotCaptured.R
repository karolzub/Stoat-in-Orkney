### R script to run the model multiple times
### with sensitivity analysis

set.seed(200)

library(compiler)

#### This is just calling the simulation function.
source("C:/Orkney/Model/StoatSim_Dist_Orkney_NotCaptured.R")

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
n.sim<-10  ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 6)

n.day<-round(qunif(lhs1[,1], 1, 28), digits=0)

init.density<-round(qunif(lhs1[,2], 0.0026, 0.0035), digits=4)

trap.dist<-round(qunif(lhs1[,3], 50, 250), digits=0)           ### Minimum distance between traps

g0<-round(qunif(lhs1[,4], 0.01, 0.1), digits=4)           ### Probability of capture at the centre of the home range

sigma<-round(qunif(lhs1[,5], 100, 900), digits=0)           ### Spatial decay parameter

p.bycatch<-round(qunif(lhs1[,6], 0.005, 0.015), digits=5)           ### Probability of by-catch per trap day


### Empty list to store the results
stoat.sim<-vector("list")

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulations

for (i in 1:n.sim){
  
  stoat.sim[[i]]<-stoat.dist(n.day=n.day[i],init.density=init.density[i], trap.dist=trap.dist[i], 
                             
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

#write.table(HRcentersDF,"C:/Orkney/Model/Simulated-Stoats.csv", sep=",")


