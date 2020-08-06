# Analysis on characteristics of plant-pollinator associations
# trying to focus the code a bit
# assumes local environment
# last used/modified jby, 2020.06.23

# setwd('~/Documents/Academic/Active_projects/pollination_networks')

require("MASS")
require("lme4")
require("nlme")
require("MuMIn")
require("piecewiseSEM")
require("fitdistrplus")

require("brms")
require("parallel")

require("xtable")

require("tidyverse")

source("../shared/Rscripts/base.R")
source("../shared/Rscripts/base_graphics.R")

#--------------------------------------------------------------------
# load data

# metadata
meta <- read.csv("data/references_all.csv", h=TRUE)

# degree distributions! New as of 2020.03.21
degrees <- read.csv("output/degree_per_plant_phy.csv", h=TRUE) %>% filter(!is.na(symmetry)) %>% left_join(meta %>% select(web, Latitude)) %>% mutate(lat.x = sqrt(abs(Latitude)), lat.sc = (lat.x-mean(lat.x, na.rm=TRUE))/sd(lat.x, na.rm=TRUE), pc1.sc = PC1/sd(PC1, na.rm=TRUE), pc2.sc = PC2/sd(PC2, na.rm=TRUE), rel.mn.n.shared=mn.n.shared/n.poll, rel.mn.n.shared.a=mn.n.shared.a/n.poll, rel.mn.n.shared.z=mn.n.shared.z/n.poll)

# check over
head(degrees)

# identify webs for with we have at least 3 of each symmetry type
deg.mn <- degrees %>% group_by(web, symmetry) %>% summarize(mnPoll = mean(n.poll), mdPoll = median(n.poll), N=length(n.poll))

mt3 <- intersect(filter(deg.mn, symmetry=="actinomorphic", N>=3)$web, filter(deg.mn, symmetry=="zygomorphic", N>=3)$web)

deg.test <- filter(degrees, web %in% mt3, !is.na(PC1)) # dataset for working, I think

head(deg.test) # coo
dim(deg.test) # 2637


#--------------------------------------------------------------------
# analysis on NEW pollinator sharing calculations

summary(degrees$n.poll)
hist(degrees$n.poll)
length(which(degrees$n.poll>100))
hist(sqrt(degrees$n.poll))

dim(degrees) # 4035!

degrees %>% group_by(symmetry) %>% summarize(mdPoll = median(n.poll), mnPoll=mean(n.poll), mdShare=median(rel.mn.n.shared, na.rm=TRUE), mnShare=mean(rel.mn.n.shared, na.rm=TRUE), mnShareA=mean(rel.mn.n.shared.a, na.rm=TRUE), mnShareZ=mean(rel.mn.n.shared.z, na.rm=TRUE))
#   symmetry      mdPoll mnPoll mdShare mnShare mnShareA mnShareZ
#   <chr>          <dbl>  <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
# 1 Actinomorphic      5   14.4   0.179   0.232    0.236    0.189
# 2 Zygomorphic        4   11.7   0.198   0.269    0.267    0.298
# ... HUH

# comparing Act and Zyg species
wilcox.test(n.poll~symmetry, data=degrees, alt="g") # p = 0.003
wilcox.test(rel.mn.n.shared~symmetry, data=degrees, alt="l") # p = 0.003
wilcox.test(rel.mn.n.shared.a~symmetry, data=degrees, alt="l") # n.s.
wilcox.test(rel.mn.n.shared.z~symmetry, data=degrees, alt="l") # p < 2.2e-16 HRM

# comparing sharing with Act and Zyg species
wilcox.test(degrees[degrees$symmetry=="Actinomorphic","mn.n.shared.a"], degrees[degrees$symmetry=="Actinomorphic","mn.n.shared.z"], paired=TRUE, alt="g") # p < 2.2e-16
wilcox.test(degrees[degrees$symmetry=="Zygomorphic","mn.n.shared.a"], degrees[degrees$symmetry=="Zygomorphic","mn.n.shared.z"], paired=TRUE, alt="l") # p < 4.8e-05 

# comparison with latitude ...
cor.test(~n.poll+abs(Latitude), data=filter(degrees, symmetry=="actinomorphic"), method="sp") # cor = 0.39, p < 2.2e-16
cor.test(~n.poll+abs(Latitude), data=filter(degrees, symmetry=="zygomorphic"), method="sp") # cor = 0.37, p < 2.2e-16

# comparison with phylospace ...
cor.test(~n.poll+PC1, data=degrees, method="sp") # cor = -0.09, 1.108e-08
t.test(PC1~symmetry, data=degrees) # n.s.
t.test(PC2~symmetry, data=degrees) # p = 0.0001919
t.test(PC3~symmetry, data=degrees) # p < 2.2e-16
t.test(PC4~symmetry, data=degrees) # p < 2.2e-16


ggplot(filter(degrees, matrix %in% mt3), aes(x=PC1, y=PC2, color=symmetry)) + geom_point() + theme(legend.position=c(0.25,0.75))
ggplot(filter(degrees, matrix %in% mt3), aes(x=PC1, y=PC3, color=symmetry)) + geom_point() + theme(legend.position=c(0.25,0.75))
ggplot(filter(degrees, matrix %in% mt3), aes(x=PC1, y=PC4, color=symmetry)) + geom_point() + theme(legend.position=c(0.25,0.75))
ggplot(filter(degrees, matrix %in% mt3), aes(x=PC2, y=PC3, color=symmetry)) + geom_point() + theme(legend.position=c(0.25,0.75))
ggplot(filter(degrees, matrix %in% mt3), aes(x=PC2, y=PC4, color=symmetry)) + geom_point() + theme(legend.position=c(0.25,0.75))
ggplot(filter(degrees, matrix %in% mt3), aes(x=PC3, y=PC4, color=symmetry)) + geom_point() + theme(legend.position=c(0.25,0.75))


dim(filter(degrees, symmetry=="actinomorphic")) # 3420
dim(filter(degrees, symmetry=="zygomorphic")) # 615

#--------------------------------------------------------------------
# okay formal model-fitting

# number of pollinators per taxon ----------------
hist(deg.test$n.poll) # hmmm
sum(deg.test$n.poll==0) # no zeroes
descdist(deg.test$n.poll, boot=100) # hmmm

npoll.M <- brm(n.poll ~ (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MS <- brm(n.poll ~ symmetry + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.ML <- brm(n.poll ~ lat.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MP <- brm(n.poll ~ pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MSL <- brm(n.poll ~ symmetry + lat.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MSxL <- brm(n.poll ~ symmetry * lat.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MSP <- brm(n.poll ~ symmetry + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MLP <- brm(n.poll ~ lat.sc + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MSLP <- brm(n.poll ~ symmetry + lat.sc + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
npoll.MSxLP <- brm(n.poll ~ symmetry * lat.sc + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)

npoll.mods <- list(npoll.M, npoll.MS, npoll.ML, npoll.MP, npoll.MSL, npoll.MSxL, npoll.MSP, npoll.MLP, npoll.MSLP, npoll.MSxLP)

save(npoll.mods, file="output/brmmods/npoll.mods.RData")

load("output/brmmods/npoll.mods.RData")

npoll.M <- npoll.mods[[1]]
npoll.MS <- npoll.mods[[2]]
npoll.ML <- npoll.mods[[3]]
npoll.MP <- npoll.mods[[4]]
npoll.MSL <- npoll.mods[[5]]
npoll.MSxL <- npoll.mods[[6]]
npoll.MSP <- npoll.mods[[7]]
npoll.MLP <- npoll.mods[[8]]
npoll.MSLP <- npoll.mods[[9]]
npoll.MSxLP <- npoll.mods[[10]]

npoll.comp <- LOO(npoll.M, npoll.MS, npoll.ML, npoll.MP, npoll.MSL, npoll.MSxL, npoll.MSP, npoll.MLP, npoll.MSLP, npoll.MSxLP) 

npoll.comp # best-fit is MSxLP

save(npoll.comp, file="output/brmmods/npoll.comp.RData")

load("output/brmmods/npoll.comp.RData")

npoll.comp

bayes_R2(npoll.MSxLP) # rsq = 0.28

summary(npoll.MSxLP)
# negative effect of zygomorphy; positive effect of lat; zygomorphy reduces lat effect


# and relative shared pollinators --------------------

descdist(deg.test$rel.mn.n.shared, boot=100) # hmmm

# Bayesian Regression Models using Stan
Rshare.M <- brm(rel.mn.n.shared ~ (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MS <- brm(rel.mn.n.shared ~ symmetry + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.ML <- brm(rel.mn.n.shared ~ lat.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MP <- brm(rel.mn.n.shared ~ pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MSL <- brm(rel.mn.n.shared ~ symmetry + lat.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MSxL <- brm(rel.mn.n.shared ~ symmetry * lat.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MSP <- brm(rel.mn.n.shared ~ symmetry + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MLP <- brm(rel.mn.n.shared ~ lat.sc + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MSLP <- brm(rel.mn.n.shared ~ symmetry + lat.sc + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)
Rshare.MSxLP <- brm(rel.mn.n.shared ~ symmetry * lat.sc + pc1.sc + pc2.sc + (1|web), data=deg.test, family=hurdle_gamma, iter=5000, cores=4)

Rshare.mods <- list(Rshare.M, Rshare.MS, Rshare.ML, Rshare.MP, Rshare.MSL, Rshare.MSxL, Rshare.MSP, Rshare.MLP, Rshare.MSLP, Rshare.MSxLP)

save(Rshare.mods, file="output/brmmods/Rshare.mods.RData")

load("output/brmmods/Rshare.mods.RData")

Rshare.M <- Rshare.mods[[1]]
Rshare.MS <- Rshare.mods[[2]]
Rshare.ML <- Rshare.mods[[3]]
Rshare.MP <- Rshare.mods[[4]]
Rshare.MSL <- Rshare.mods[[5]]
Rshare.MSxL <- Rshare.mods[[6]]
Rshare.MSP <- Rshare.mods[[7]]
Rshare.MLP <- Rshare.mods[[8]]
Rshare.MSLP <- Rshare.mods[[9]]
Rshare.MSxLP <- Rshare.mods[[10]]

Rshare.comp <- LOO(Rshare.M, Rshare.MS, Rshare.ML, Rshare.MP, Rshare.MSL, Rshare.MSxL, Rshare.MSP, Rshare.MLP, Rshare.MSLP, Rshare.MSxLP, k=10) 

Rshare.comp # best-fit is MSP; nothing with overlapping CI

save(Rshare.comp, file="output/brmmods/Rshare.comp.RData")
load("output/brmmods/Rshare.comp.RData")

summary(Rshare.MSP)
# positive effect of zygomorphy

bayes_R2(Rshare.MSP) # rsq = 0.66


