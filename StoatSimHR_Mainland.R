### R script - Model to simulate home ranges and traps within the home range.
### Data on trap deployment for South Ronaldsay (fixed number of traps = 406)
### Parameters: init.density: stoat density per ha
### tot.area: area (ha)
### n.trap: number of traps
### n.day: number of trapping days
### cor.mov: correlation in the movements of stoats between one day and the next
### min.det: minimum probability of capture
### max.det: maximum probability of capture

stoat.dist<-function(n.day, init.density, n.trap, trap.dist,
                     cor.mov, max.det, p.bycatch){

### Packages required.
        require(MASS)
        require(MBESS)
        require(adehabitatHR)
        require(sp)
        require(spatstat)
        require(raster)
        require(geosphere)
        require(maptools)
        require(rgdal)
        require(sf)
        require(beepr)
        require(rgeos)   
        
        #n.day <- 20
        #cor.mov <- 0.3
        #min.det <- 0.001
        #max.det <- 0.5
        #p.bycatch <- 0.008
        #n.trap <- 5000
        
        ## Total number of stoats
        #init.density<-2/100
        tot.area <- 52325
        tot.stoat<-rpois(1, init.density*tot.area)
        #tot.stoat
        
        ### All potential home range centres
        HRcenters <- spsample(main, n = tot.stoat, "random", iter = 20)
        HRcentersDF <- as.data.frame(HRcenters)
        #plot(HRcenters, add = TRUE, col="red", pch=19, cex=1)
        
        ### All potential trap locations
        #### TRAPS
        centre.trap.rand <- spsample(access, n = 10000, "random", iter = 20)
        #centre.trap <- traps_bord[sample(nrow(traps_bord), 5000), ]
        centre.trap <- as.data.frame(centre.trap.rand)
        coordinates(centre.trap)  <-  c("x", "y")
        #plot(centre.trap, add = TRUE, col="blue", pch=19, cex=1)
        ### Number of traps
        n.trap<-length(centre.trap$x)
        #n.trap<-length(centre.trap$coords.x1)
        #n.trap
        #head(centre.trap)
        
        ### Simulate stoat home ranges in n.day from a multivariate normal distribution
        ### Covariance between x and y (values in metres); the same for all stoat (correlation between x an y = 0.5 in this example)
        cor.mat<-matrix(cor.mov, 2, 2)
        
        diag(cor.mat)<-1
        
        ### Convert to variance-covariance matrix; standard deviation in metres
        cov.mat<-cor2cov(cor.mat, sd=c(500, 500))
        
       ### Simulate stoat positions during n.days
                stoat.pos<-vector("list")

                for (i in 1:tot.stoat){

                        stoat.out<-matrix(0, ncol=3, nrow=n.day)

                        stoat.out[ ,3]<-as.numeric(i)

                        stoat.out[ , 1:2]<-mvrnorm(n.day, mu=c(HRcentersDF[i, 1], HRcentersDF[i, 2]), Sigma=cov.mat)      

                        stoat.pos[[i]]<-stoat.out   
                }

            ### Put together in a matrix

            stoat.mov<-do.call(rbind, stoat.pos)

            colnames(stoat.mov)<-c("x", "y", "ID")

            stoat.mov<-as.data.frame(stoat.mov)

            coordinates(stoat.mov)<-stoat.mov[, c('x', 'y')]

            ### Estimate 100% MCP
            mcp.stoat<-mcp(stoat.mov[, 3], percent=100)

            ### Check that the simulated home range area mkes sense (stoats: btw 100 and 500 ha); you may need to tweak the sd in line 61 to make this work
            hr.area<-as.data.frame(mcp.stoat)
            #summary(hr.area)

            ###### Loop for estimating overlap between each home range and each camera deployed and the detectability (assumes variability per individual stoat not trap)
            p.det<-runif(tot.stoat, min.det, max.det)      ### Probability of detection per night

            ### Simulate the data

            over.camera<-det.camera<-matrix(0, ncol=tot.stoat, nrow=n.trap)

            for (i in 1:tot.stoat){
                    
                    stoat.coord<-mcp.stoat@polygons[[i]]@Polygons[[1]]@coords				### Coordinates of the home range of each stoat
                    
                    over.cam1<-point.in.polygon(centre.trap$x, centre.trap$y, as.vector(stoat.coord[,1]), as.vector(stoat.coord[,2]), mode.checked=FALSE)			#### Check which cameras are within the stoat HR
                    
                    over.camera[,i]<-ifelse(over.cam1>0, 1, 0)						#### Cameras within the home range are 1 other 0
                    
                    by.catch<-rbinom(n.trap, 1, prob=1-(1-p.bycatch)^n.day)                     ### By-catch in each trap after a given number of days (0 or 1)                         
                    
                    int.tot<-1-(1-p.det[i])^n.day              ### Probability of interacting with a trap after the total number of nights
                    
                    ### Detection in a given trap (0, 1) - the probability of capture depends on three factors: (i) being within a home range;
                    ### (ii) the stoat interacting with the trap (int.tot); and (iii) the cumulative probability of bycatch in that trap
                    det.camera[, i]<-rbinom(length(over.camera[, i]), size=1, prob=over.camera[, i]*int.tot*(1-by.catch))     
                    
            }

        ### Stoat that will be caught

        stoat.camera<-apply(det.camera, 2, max)

        ### Percentage of stoat population within the range of each camera trap

        perc.catch<-round((sum(stoat.camera)/tot.stoat)*100, digits=2)

        ### Output of the function - returns the mean home range size (ha) and the percentage of stoats caught

        out.sim<-round(c(mean(hr.area$area), perc.catch), digits=2)
        
        #out.sim <- (perc.catch)
        
        return(out.sim)

}



#plot(centre.stoat$x, centre.stoat$y, pch=19, col="black", cex=1)

#plot(mcp.stoat, add=TRUE)

#points(centre.trap$x, centre.trap$y, col="blue", pch=19, cex=1)           ### Plot the traps
