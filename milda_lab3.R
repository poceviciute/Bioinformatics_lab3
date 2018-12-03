# Question 1

###
### Chapter 3
###
#setwd("Z:/Documents/Bioinformatics_lab3")
library(msa)
x <- paste("AJ5345", 26:49, sep = "")
x <- c("Z73494", x)
sylvia.seq <- read.GenBank(x)
ape::write.dna(sylvia.seq, file ="sylvia_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)

sylvia_seqs<-ape::read.FASTA("sylvia_seqs.fasta")
#sylvia.seq <- sylvia_seqs
# Align the sequences:
sylvia.clus <- clustal(sylvia.seq)
aligment_sylvia <- msa("sylvia_seqs.fasta", type="dna")
sylvia.seq.ali <- as.DNAbin(aligment_sylvia)
rownames(sylvia.seq.ali) <- labels(sylvia.seq)

identical(sylvia.clus[x, ], sylvia.seq.ali[x, ])

taxa.sylvia <- attr(sylvia.seq, "species")
names(taxa.sylvia) <- names(sylvia.seq)
rm(sylvia.seq)
taxa.sylvia[1] <- "Sylvia_atricapilla"
taxa.sylvia[24] <- "Sylvia_abyssinica"

#setwd("D:/LiU/732A51/Bioinformatics_lab3")

sylvia.eco <- read.table("sylvia_data.txt")
str(sylvia.eco)
rownames(sylvia.eco)

save(sylvia.clus, taxa.sylvia, sylvia.eco,
     file = "sylvia.RData")

###
### Chapter 5
###

syl.K80 <- ape::dist.dna(sylvia.seq.ali, pairwise.deletion = TRUE)
syl.F84 <- ape::dist.dna(sylvia.seq.ali, model = "F84", p = TRUE)
syl.TN93 <- ape::dist.dna(sylvia.seq.ali, model = "TN93", p = TRUE)
syl.GG95 <- ape::dist.dna(sylvia.seq.ali, model = "GG95", p = TRUE)

round(cor(cbind(syl.K80, syl.F84, syl.TN93, syl.GG95)), 3)

syl.JC69 <- ape::dist.dna(sylvia.seq.ali, model = "JC69", p = TRUE)
syl.raw <- ape::dist.dna(sylvia.seq.ali, model = "raw", p = TRUE)
layout(matrix(1:2, 1))
plot(syl.JC69, syl.raw)
abline(b = 1, a = 0) # draw x = y line
plot(syl.K80, syl.JC69)
abline(b = 1, a = 0)

layout(matrix(1:3, 1))
for (i in 1:3) {
  s <- logical(3); s[i] <- TRUE
  x <- sylvia.seq.ali[, s]
  d <- dist.dna(x, p = TRUE)
  ts <- dist.dna(x, "Ts", p = TRUE)
  tv <- dist.dna(x, "Tv", p = TRUE)
  plot(ts, d, xlab = "Number of Ts or Tv", col = "blue",
       ylab = "K80 distance", xlim = range(c(ts, tv)),
       main = paste("Position", i))
  points(tv, d, col = "red")
}

y <- numeric()
for (i in 1:3) {
  s <- logical(3); s[i] <- TRUE
  y <- c(y, dist.dna(sylvia.seq.ali[, s], p = TRUE))
}
g <- gl(3, length(y) / 3)
library(lattice)
histogram(~ y | g, breaks = 20)

nj.sylvia.K80 <- nj(syl.K80)
nj.sylvia.GG95 <- nj(syl.GG95)
dist.topo(nj.sylvia.K80, nj.sylvia.GG95)

grep("Chamaea", taxa.sylvia, value = TRUE)
f <- function(xx) root(nj(dist.dna(xx, p=TRUE)), "AJ534526")
tr <- f(sylvia.seq.ali)
## same than: tr <- root(nj.sylvia.K80, "AJ534526")
nj.boot.sylvia <- boot.phylo(tr, sylvia.seq.ali, f, 200,
                             rooted = TRUE)
nj.boot.codon <- boot.phylo(tr, sylvia.seq.ali, f, 200, 3,
                            rooted = TRUE)
nj.est <- tr
nj.est$tip.label <- taxa.sylvia[tr$tip.label]
plot(nj.est, no.margin = TRUE)
nodelabels(round(nj.boot.sylvia / 200, 2), bg = "white")
add.scale.bar(length = 0.01)
write.tree(nj.est, "sylvia_nj_k80.tre")

# WHY ITS NOT WORKING - WHERE TO GET sylvia.txt_phyml_tree.txt AND sylvia.txt_phyml_stats.txt
write.dna(sylvia.seq.ali, "sylvia.txt")
phyml.sylvia <- phymltest("sylvia.txt", execname = "~/phyml")

summary(phyml.sylvia)
plot(phyml.sylvia, col = "black")
TR <- read.tree("sylvia.txt_phyml_tree.txt")
mltr.sylvia <- TR[[28]]
mltr.sylvia$tip.label <- taxa.sylvia[mltr.sylvia$tip.label]
mltr.sylvia <- root(mltr.sylvia, "Chamaea_fasciata")
plot(mltr.sylvia, no.margin = TRUE)
add.scale.bar(length = 0.01)

tr.ml <- drop.tip(mltr.sylvia, "Chamaea_fasciata")
res <- vector("list", 9)
for (L in -4:4)
  res[[L + 5]] <- chronopl(tr.ml, 10^L, 12, 16, CV = TRUE)
Lambda <- 10^(-4:4)
CV <- sapply(res, function(x) sum(attr(x, "D2")))
plot(Lambda, CV / 1e5, log = "x")

sylvia.chrono <- res[[2]]
rts <- attr(sylvia.chrono, "rates")
summary(rts)

par(mar = c(2, 0, 0, 0))
plot(sylvia.chrono, edge.width = 100*rts, label.offset = .15)
axisPhylo()
write.tree(sylvia.chrono, "sylvia.chrono.tre")

###
### Chapter 6
###

load("sylvia.RData")
nj.est <- read.tree("sylvia_nj_k80.tre")
nj.est <- drop.tip(nj.est, "Chamaea_fasciata")
DF <- sylvia.eco[nj.est$tip.label, ]
table(DF$geo.range, DF$mig.behav)

syl.er <- ace(DF$geo.range, nj.est, type = "d")
syl.sym <- ace(DF$geo.range, nj.est, type="d", model="SYM")
anova(syl.er, syl.sym)

mod <- matrix(0, 3, 3)
mod[2, 1] <- mod[1, 2] <- 1
mod[2, 3] <- mod[3, 2] <- 2
syl.mod <- ace(DF$geo.range, nj.est, type="d", model=mod)

sapply(list(syl.er, syl.sym, syl.mod), AIC)

Q <- syl.mod$index.matrix
diag(Q) <- 0
Q[1, 2] <- Q[2, 1] <- syl.mod$rates[1]
Q[2, 3] <- Q[3, 2] <- syl.mod$rates[2]

Q[] <- c(0, syl.mod$rates)[Q + 1]
diag(Q) <- -rowSums(Q)

P <- matexpo(0.05 * Q)
rownames(P) <- c("temp", "temptrop", "trop")
colnames(P) <- rownames(P)

co <- rep("grey", 24)
co[DF$geo.range == "temp"] <- "black"
co[DF$geo.range == "trop"] <- "white"
plot(nj.est, "c", FALSE, no.margin = TRUE, label.offset = 1)
tiplabels(pch = 22, bg = co, cex = 2, adj = 1)
nodelabels(thermo = syl.mod$lik.anc, cex = 0.8,
           piecol = c("black", "grey", "white"))

sylvia.chrono <- read.tree("sylvia.chrono.tre")
yule(sylvia.chrono)
birthdeath(sylvia.chrono)
1 - pchisq(2*(-1.034112 - -1.113822), 1)

x <- sylvia.eco[sylvia.chrono$tip.label, "geo.range"]
ANC <- ace(x, sylvia.chrono, type = "d", model = mod)
ANC$lik.anc[1:3, ]
anc <- apply(ANC$lik.anc, 1, which.max)
X <- factor(c(x, anc))
yule.cov(sylvia.chrono, ~ X)
1 / (1 + exp(-(-0.0535529)))
1 / (1 + exp(-(-0.0535529 -1.4608019)))
1 / (1 + exp(-(-0.0535529 -0.9775966)))

fsamp <- function(x) sample(length(x), size = 1, prob = x)
nrep <- 1e3
Pvls <- numeric(nrep)
for (i in 1:nrep) {
  anc <- apply(ANC$lik.anc, 1, fsamp)
  X <- factor(c(x, anc))
  Pvls[i] <- yule.cov(sylvia.chrono, ~ X)$Pval
}
hist(Pvls, freq = FALSE, main = "")
lines(density(Pvls))

# Question 2
#install.packages("ade4")
#install.packages("mvMORPH")
#install.packages("mvSLOUCH")
#install.packages("adephylo")
library(ade4)
library(mvMORPH)
library(mvSLOUCH)
library(ape)
#devtools::install_github("kopperud/slouch")
#install.packages("rlang")
#install.packages("tibble")
library(tibble)
library(rlang)

# Q2.1
data(carni70)

View(carni70)

carni70_tree <- newick2phylog(carni70$tre)


size <- scalewt(log(carni70$tab))[,1] # scaled log body size in kg
names(size) <- row.names(carni70$tab)
yrange <- scalewt(carni70$tab[,2])# scaled geographic range in km
names(yrange) <- row.names(carni70$tab)

plot_data <- data.frame(spieces = row.names(carni70$tab), size = size, range =yrange )
par(mfrow=c(1,1))

plot(plot_data$spieces, plot_data$size, xlab = "Spieces", ylab = "Size", pch = 19,
     col = rgb(0, 0, 0, 0.1))

plot(plot_data$spieces, plot_data$range, xlab = "Spieces", ylab = "Size", pch = 19,
     col = rgb(0, 0, 0, 0.1))

symbols.phylog(carni70.phy,size)

tre <- ape::read.tree(text = carni70$tre)
adephylo::orthogram(size, tre = tre)


symbols.phylog(carni70.phy,yrange)
adephylo::orthogram(as.vector(yrange), tre = tre)

if(adegraphicsLoaded()) {
  g1 <- s.label(cbind.data.frame(size, yrange), plabel.cex = 0)
  g2 <- addhist(g1)
} else {
  s.hist(cbind.data.frame(size, yrange), clabel = 0)
}


# This data set describes the phylogeny of 70 carnivora as reported by Diniz-Filho and Torres (2002). 
# It also gives the geographic range size and body size corresponding to these 70 species.
# tre: is a character string giving the phylogenetic tree in Newick format. Branch lengths 
# are expressed as divergence times (millions of years);
# tab: s a data frame with 70 species and two traits: size (body size (kg)) ; range (geographic range size (km))

# Q2.2
tre <- ape::read.tree(text = carni70$tre)
fit1 <- mvBM(tre, data = carni70$tab, model="BM1")

# library(nlme)
brownian_corr <- corBrownian(phy = tre)
# m1 <- gls(size ~ yrange, correlation=brownian_corr)
# summary(m1)

ace(size, as.phylo(carni70_tree))
ace(size, tre, method = "pic")
ace(x, tre, method = "GLS",
    corStruct = corBrownian(1, bird.orders))
# Both traits evolve as independent Brownian motions.

# The traits evolve as a correlated Brownian motion.

# Both traits evolve as independent Ornstein{Uhlenbeck processes
taxa <- carni70$tab
tre_phylo <- as.phylo(tre)
fit2 <- corphylo(X = taxa, phy = tre_phylo)

# The traits evolve as a bivariate Ornstein{Uhlenbeck process

# size evolves as a Brownian motion and range as an Ornstein{Uhlenbeck process adapting
# to it (use slouch or mvSLOUCH and be careful about column order).
