rm(list = ls())

# question 1

library(msa)
library(ape)

x <- paste("AJ5345", 26:49, sep = "")
x <- c("Z73494", x)

#read a data from GenBank
sylvia.seq <- read.GenBank(x)
sylvia.seq
#write.dna(sylvia.seq, file ="sylvia_seqs.fasta", format = "fasta")

#syliva sequences
#sylvia_seqs <- read.FASTA("sylvia_seqs.fasta")
sylvia.clus <- clustal(sylvia.seq)

library(phyloch)
#sylvia.maff <- mafft(sylvia.seq)
#identical(sylvia.clus[x, ], sylvia.maff[x, ])


taxa.sylvia <- attr(sylvia.seq, "species")
names(taxa.sylvia) <- names(sylvia.seq)

taxa.sylvia[1] <- "Sylvia_atricapilla"
taxa.sylvia[24] <- "Sylvia_abyssinica"



sylvia.eco <- read.table("sylvia_data.txt")
str(sylvia.eco)
rownames(sylvia.eco)

#save sylvia cluster
save(sylvia.clus, taxa.sylvia, 
     sylvia.eco,
     file = "sylvia.RData")


###
### Chapter 5
###

### my added code 
write.dna(sylvia.seq, file ="sylvia_seqs.fasta", format = "fasta")
sylvia_aligment <- msa("sylvia_seqs.fasta", type="dna")
sylvia.seq.ali <- as.DNAbin(sylvia_aligment)
rownames(sylvia.seq.ali) <- names(sylvia.seq)

#####

syl.K80 <- dist.dna(sylvia.seq.ali, pairwise.deletion = TRUE)
syl.F84 <- dist.dna(sylvia.seq.ali, model = "F84", p = TRUE)
syl.TN93 <- dist.dna(sylvia.seq.ali, model = "TN93", p = TRUE)
syl.GG95 <- dist.dna(sylvia.seq.ali, model = "GG95", p = TRUE)
syl.JC69 <- dist.dna(sylvia.seq.ali, model = "JC69", p = TRUE)
syl.raw <- dist.dna(sylvia.seq.ali, model = "raw", p = TRUE)
# 
# layout(matrix(1:2, 1))
# plot(syl.JC69, syl.raw)
# abline(b = 1, a = 0) # draw x = y line
# plot(syl.K80, syl.JC69)
# abline(b = 1, a = 0)




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

layout(matrix(1:1, 1))

y <- numeric()
for (i in 1:3) {
  s <- logical(3); s[i] <- TRUE
  y <- c(y, dist.dna(sylvia.seq.ali[, s], p = TRUE))
}
y
#generate a factor level
g <- gl(3, length(y) / 3)
library(lattice)
histogram(~ y | g, breaks = 20)


nj.sylvia.K80 <- nj(syl.K80)
nj.sylvia.GG95 <- nj(syl.GG95)
dist.topo(nj.sylvia.K80, nj.sylvia.GG95)

#finding a Chamaea Species from Taxa Syliva list
grep("Chamaea", taxa.sylvia, value = TRUE)

f <- function(xx) root(nj(dist.dna(xx, p=TRUE)), "AJ534526") ## rerooting nj tree wrt AJ534526
tr <- f(sylvia.seq.ali)
## same than: tr <- root(nj.sylvia.K80, "AJ534526")
nj.boot.sylvia <- boot.phylo(tr, sylvia.seq.ali, f, 200,
                             rooted = TRUE)
nj.boot.codon <- boot.phylo(tr, sylvia.seq.ali, f, 200, 3,
                            rooted = TRUE)
nj.est <- tr
#updaing tip lable of nj estimated to taxa.syliva tips labels
nj.est$tip.label <- taxa.sylvia[tr$tip.label]
plot(nj.est, no.margin = TRUE)
nodelabels(round(nj.boot.sylvia / 200, 2), bg = "white")
add.scale.bar(length = 0.01)
write.tree(nj.est, "sylvia_nj_k80.tre")


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
# 
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

#question 1 part 2 code

load("sylvia.RData")
nj.est <- read.tree("sylvia_nj_k80.tre")
#droping a tip "Chamaea_fasciata"
nj.est <- drop.tip(nj.est, "Chamaea_fasciata")
DF <- sylvia.eco[nj.est$tip.label, ]
table(DF$geo.range, DF$mig.behav)

syl.er <- ace(DF$geo.range, nj.est, type = "discrete",model = "ER")
syl.sym <- ace(DF$geo.range, nj.est, type="discrete", model="SYM")
syl.ard <- ace(DF$geo.range, nj.est, type="d", model="ARD")

syl.er$se
syl.sym$se
anova(syl.er, syl.sym)


#syl.mod <- ace(DF$geo.range, nj.est, type="d", model=mod)

sapply(list(syl.er, syl.sym, syl.ard), AIC)

Q <- syl.er$index.matrix
diag(Q) <- 0
Q[1, 2] <- Q[2, 1] <- syl.er$rates[1]
Q[2, 3] <- Q[3, 2] <- syl.er$rates[1]

#Q[] <- c(0, syl.mod$rates)[Q + 1]
diag(Q) <- -rowSums(Q)

P <- matexpo(0.05 * Q)
rownames(P) <- c("temp", "temptrop", "trop")
colnames(P) <- rownames(P)
P
co <- rep("grey", 24)
co[DF$geo.range == "temp"] <- "black"
co[DF$geo.range == "trop"] <- "white"
plot(nj.est, "c", FALSE, no.margin = TRUE, label.offset = 1)
tiplabels(pch = 22, bg = co, cex = 2, adj = 1)
nodelabels(thermo = syl.mod$lik.anc, cex = 0.8,
           piecol = c("black", "grey", "white"))




# sylvia.chrono <- read.tree("sylvia.chrono.tre")
# yule(sylvia.chrono)
# birthdeath(sylvia.chrono)
# 1 - pchisq(2*(-1.034112 - -1.113822), 1)
# 
# x <- sylvia.eco[sylvia.chrono$tip.label, "geo.range"]
# ANC <- ace(x, sylvia.chrono, type = "d", model = mod)
# ANC$lik.anc[1:3, ]
# anc <- apply(ANC$lik.anc, 1, which.max)
# X <- factor(c(x, anc))
# yule.cov(sylvia.chrono, ~ X)
# 1 / (1 + exp(-(-0.0535529)))
# 1 / (1 + exp(-(-0.0535529 -1.4608019)))
# 1 / (1 + exp(-(-0.0535529 -0.9775966)))
# 
# fsamp <- function(x) sample(length(x), size = 1, prob = x)
# nrep <- 1e3
# Pvls <- numeric(nrep)
# for (i in 1:nrep) {
#   anc <- apply(ANC$lik.anc, 1, fsamp)
#   X <- factor(c(x, anc))
#   Pvls[i] <- yule.cov(sylvia.chrono, ~ X)$Pval
# }
# hist(Pvls, freq = FALSE, main = "")
# lines(density(Pvls))
# 
