library(ade4)

data("carni70")

### 2.1


carni70_phy <- newick2phylog(carni70$tre)
plot(carni70_phy)

size <- scalewt(log(carni70$tab))[,1] # scale log body size
names(size) <- row.names(carni70$tab)
symbols.phylog(carni70_phy,size)


tre <- ape::read.tree(text = carni70$tre)
adephylo::orthogram(size, tre = tre)

yrange <- scalewt(carni70$tab[,2]) # scale log geographic range size
names(yrange) <- row.names(carni70$tab)
symbols.phylog(carni70_phy,yrange)

### It can be observed from the above trees that body size and geographic range size both have most have the values
## that are too small, that is either smaller than mean. Even the values which are greater than mean are 
## in the smaller range that is very few values with 2.5 size

adephylo::orthogram(as.vector(yrange), tre = tre)

## It is obvious that the body weight and geographic range size are traits with continuous character.
## thus the usual analyzing techniques will not be used here, instead we will have to use a model
## which would analyze the continuous data.

a <- carni70$tab[1]
hist(a$size,breaks = 50)


### 2.2
library(mvSLOUCH)
carni_tree <- ape::read.tree(text = carni70$tre)

#BM_fit_1 <- BrownianMotionModel(carni_tree,carni70$tab)
