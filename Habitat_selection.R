### Set of useful functions to analyse habitat selection by stoats using Compositional Analysis

library(lmerTest)
library(popbio)
library(logistf)
library(spdep)
library(multcomp)
library(blme)
library(adehabitatHS)

# text file with used habitats data, 
# no rows names, every row corresponds to different buffer (around trap or sighting site)
# every column corresponds to habitat type
used<-read.table("C:/Orkney/Model/Habitat selection/used.txt", header=TRUE, na.strings = "NA")
head(used)

# text file with avialble habitats data,
# number of rows equals the number of rows in used habiatats
# could be the same proportion of available habitats for all rows
# every column corresponds to habitat type
avail<-read.table("C:/Orkney/Model/Habitat selection/avail.txt", header=TRUE, na.strings = "NA")
head(avail)

#main test
# test = "parametric" or "randomisation"
# rnv - replacement for 0 values
# nrep - number of replications
# alpha - significance level
hs <- compana(used, avail, test = c("parametric"), rnv = 0.001, nrep = 1000, alpha = 0.05)
# graphical representation of results
hs$profile
# pairwise compariosn of habitats
hs$rm
# ranking of habiatats (higest values for most preferred) 
hs$rank
# global lamda test for difference between available and used
hs$test

# useful function summarising area of habitats 
# which can be exported ot Excell to calculate proportion

dat1 <- xtabs(formula=Area_hab~Habitat, data=hs)

write.table(dat1, file="C:/Orkney/Model/Habitat selection/dat1.txt", append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

sub <- subset (hs, PARISH_COD == "LH", data = hs)

summary(sub)

# useful function summarising area of habitats separately for each location (either trap or sighting site)
# which can be exported ot Excell to calculate proportion

xtabs(formula=Area_hab~TRAP_ID + Habitat,  data=hs)

dat1 <- xtabs(formula=Area_hab~ID_2 + Habitat, data=hs)

write.table(dat1, file="C:/Orkney/Model/Habitat selection/dat1.txt", append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)


