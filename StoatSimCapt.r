### R script - Model to simulate home ranges and traps within the home range.
### Parameters: density.km2: stoat density per km2
### tot.area: area (km2)
### x.coord: x-coordinates of the spatial unit being analysed (in metres)
### y.coord: y-coordidates of the spatial unit being analysed (in metres)
### p.bycatch: nightly probability of by-catch
### n.trap: number of traps
### n.day: number of trapping days
### cor.mov: correlation in the movements of stoats between one day and the next
### min.det: minimum probability of capture
### max.det: maximum probability of capture

stoat.hr<-function(density.km2, tot.area, border.size, x.coord, y.coord, p.bycatch, n.trap, n.day, cor.mov, min.det, max.det){

            ### Packages required.
            require(MASS)
            require(MBESS)
            require(adehabitatHR)
            require(sp)

            ## Total number of stoats
            init.density<-density.km2/100			#### Density in individuals per ha

            tot.stoat<-init.density*tot.area

            ### Sample the centre of the home range of each individual stoat
            centre.stoat<-matrix(0, nrow=tot.stoat, ncol=2)           

            centre.stoat[ ,1]<-sample(x.coord, tot.stoat)

            centre.stoat[, 2]<-sample(y.coord, tot.stoat)

            ### Locate the traps at random in the landscape and at a minimum distance of 100 metres to try to ensure independence

            x.camera<-sample(round(x.coord, digits=-2), n.trap)

            y.camera<-sample(round(y.coord, digits=-2), n.trap)

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

                        stoat.out[ , 1:2]<-mvrnorm(n.day, mu=c(centre.stoat[i, 1], centre.stoat[i, 2]), Sigma=cov.mat)      

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

            ###### Loop for estimating overlap between each home range and each camera deployed and the detectability (assumes variability per individual stoat not trap)
            p.det<-runif(tot.stoat, min.det, max.det)      ### Probability of detection per night

            ### Simulate the data

            over.camera<-det.camera<-matrix(0, ncol=tot.stoat, nrow=n.trap)

            for (i in 1:tot.stoat){

                    stoat.coord<-mcp.stoat@polygons[[i]]@Polygons[[1]]@coords				### Coordinates of the home range of each stoat

		    over.cam1<-point.in.polygon(x.camera, y.camera, as.vector(stoat.coord[,1]), as.vector(stoat.coord[,2]), mode.checked=FALSE)			#### Check which cameras are within the stoat HR
       
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

        return(out.sim)

}

