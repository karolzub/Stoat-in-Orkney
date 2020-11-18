### R script to run the model multiple times

set.seed(200)

library(compiler)

#### This is just calling the simulation function - change it! This is the 'StoatSimCapt.r' script
source("C:/Orkney/Model/Pablo_scripts/StoatSimCapt.r")

stoat.hr<-cmpfun(stoat.dist)

############# SIMULATIONS
### Define the study area (same across simulations)
tot.area<-4971.145                      ### Area covered in ha (South Ronaldsay)

border.size<-sqrt(tot.area*10000)

x.coord<-c(0:border.size)

y.coord<-c(0:border.size)

############# NUMERICAL SIMULATIONS
n.sim<-200     ### Number of simulations

### Latin hypercube desing to explore a wide variety & good coverage of values
library(lhs)

lhs1<-randomLHS(n.sim, 7)

init.density<-round(qunif(lhs1[,1], 0.01, 0.1), digits=3)  

n.trap<-round(qunif(lhs1[,2], 100, 1000), digits=0)           ### Total number of traps

cor.mov<-round(qunif(lhs1[,3], 0.001, 0.8), digits=2)           ### Position correlation between days across stoats

min.det<-round(qunif(lhs1[,4], 0, 0.001), digits=2)           ### Minimum probability of capture per night (across stoats)

max.det<-round(qunif(lhs1[,5], 0.11, 0.80), digits=2)           ### Maximum probability of capture per night (across stoats)

trap.dist<-round(qunif(lhs1[,6], 50, 250), digits=0)           ### Minimum distance between traps

p.bycatch<-round(qunif(lhs1[,7], 0.01, 0.95), digits=2)         #### Probability of by-catch per night


### Empty matrix to store the results
stoat.sim<-matrix(0, nrow=20, ncol=n.sim)

#### Timing the function
start.time<-Sys.time()

### Progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

### Run the simulatoins
for (i in 1:n.sim){
  
  stoat.sim[ , i]<-stoat.hr(init.density=init.density[i], tot.area=tot.area, border.size=border.size, x.coord=x.coord, y.coord=y.coord, 
                            p.bycatch=p.bycatch, trap.dist=trap.dist[i], n.trap=n.trap[i], n.day=30, cor.mov=cor.mov[i], min.det=min.det[i], max.det=max.det[i])
  
  setTxtProgressBar(pb,i)
  
}

close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

### Exploratory plotting
par(mfrow=c(2, 3))

plot(stoat.sim[2, ]~n.trap, xlab="Number of traps", ylab="Percentage caught")

plot(stoat.sim[2, ]~cor.mov, xlab="Stoat movement correlation", ylab="Percentage caught")

plot(stoat.sim[2, ]~max.det, xlab="Maximum probability of capture per night", ylab="Percentage caught")

plot(stoat.sim[2, ]~stoat.sim[1,], xlab="Average stoat home range (ha)", ylab="Percentage caught")

#plot(stoat.sim[2, ]~stoat.sim[3,], xlab="Average distance between traps (metres)", ylab="Percentage caught")

plot(stoat.sim[2, ]~trap.dist, xlab="Average distance between traps (metres)", ylab="Percentage caught")

plot(stoat.sim[2, ]~init.density, xlab="Initial density (N/ha)", ylab="Percentage caught")


