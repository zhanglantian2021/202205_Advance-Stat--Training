# construct phylogeny

#get PhyloMaker
#install.packages("devtools")
library(devtools)
#install_github("jinyizju/V.PhyloMaker")

library(V.PhyloMaker)

# input example species list
species <- read.csv("Dipterocarpaceae.specieslist.csv")
summary(species)



# generate 10 phylogenies
tree <- phylo.maker(sp.list = species, tree = GBOTB.extended, nodes = nodes.info.1, scenarios = "S2", r = 10)    # r means the number of runs, the default is r = 1. Here r is set as 10, so 10 runs will be executed and 10 phylogenies will be built.

# write out the 10 phylogenies
for (i in 1:10)
{ 
write.tree(tree$scenario.2[[i]], paste("scenario.2_run.", i, ".tre", sep=""))
}



phy <- read.tree("scenario.2_run.1.tre")
phy
plot(phy)
