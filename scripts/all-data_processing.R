# Illustrating data set characteristics with a MAP
# from code via CJC
# assumes local environment
# last used/modified jby, 2020.03.24

# setwd('~/Documents/Academic/Active_projects/pollination_networks')

library(bipartite)
library(codependent)
library(igraph)
library(reshape)
library(stringr)

require("MASS")
require("lme4")
require("nlme")
require("MuMIn")

library(tidyverse)

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("lwgeom")
library("rgdal")

source("../shared/Rscripts/base.R")
source("../shared/Rscripts/base_graphics.R")

#--------------------------------------------------------------------
# load and parse data

# floral symmetry annotations
symm <- read.csv("data/spp_symmetry.csv", h=TRUE)

table(symm$symmetry)

spp.z <- unique(subset(symm, symmetry=="zygomorphic")$user_supplied_name)
spp.a <- unique(subset(symm, symmetry=="actinomorphic")$user_supplied_name)

intersect(spp.z, spp.a) # good to check

# list of source files from Web of Life and IWDB
files.wol <- list.files(path='data/web-of-life_2019-02-06_013943', pattern="M_PL", full=TRUE)
files.iwdb <- list.files(path='data/iwdb_20190422/matrices', pattern=".csv", full=TRUE)

# load the matrices as a list
webs.wol <- lapply(files.wol, function(x) read.csv(x, row.names=1, strip.white=TRUE))
names(webs.wol) <- gsub("data/web-of-life_.+/(.+)\\.csv", "\\1", files.wol)

webs.iwdb <- lapply(files.iwdb, function(x) read.csv(x, row.names=1, strip.white=TRUE))
names(webs.iwdb) <- gsub("[ ,]", "-", gsub("data/iwdb_.+/matrices/(.+)\\.csv", "\\1", files.iwdb))

webs <- c(webs.wol, webs.iwdb)
length(webs) # 159!

webs <- lapply(webs, function(w){ row.names(w) <- gsub("(.+) $", "\\1", row.names(w)); return(w)})


# split matrices by symmetry
webs.z <- lapply(webs, function(x) x[intersect(spp.z, rownames(x)),])
webs.a <- lapply(webs, function(x) x[intersect(spp.a, rownames(x)),])
webs.n <- lapply(webs, function(x) x[setdiff(rownames(x),c(spp.a,spp.z)),])
# issue here with trailing spaces in species names (oy) needs resolution at the source but it's turning out to be tricky ... might gain me some zyg species, tho

# what's not getting binned?
lapply(webs.n[lapply(webs.n, nrow)>0], rownames)

webs[["Robertson-1928"]] <- read.csv("data/Robertson-1928_Carlinville.csv", h=TRUE, row.names=1)


# deal with too-small sub-matrices
length(which(lapply(webs.z, nrow) >= 5)) # woah. Not a lot left!
webs.z <- webs.z[lapply(webs.z, nrow) >= 5]
webs.a <- webs.a[names(webs.z)]

# one big happy table of actual matrix descriptive stats
data.symm <- read.csv("output/sub-matrix_stats.csv", h=TRUE) # stats for symmetry-based subsets (these are CORRECTED, 2020.03.20)

data.full <- read.csv("output/whole-matrix_stats.csv", h=TRUE) # stats for whole matrices


# metadata!
meta <- read.csv("data/references_all.csv", h=TRUE, sep=",")

symm.cts <- data.frame(web=names(webs), Ntot=unlist(lapply(webs, nrow)), Nact=unlist(lapply(webs, function(x) length(which(rownames(x) %in% spp.a)))), Nzyg=unlist(lapply(webs, function(x) length(which(rownames(x) %in% spp.z))))) %>% mutate(Nnon=Ntot-(Nact+Nzyg))

data.symm <- meta %>% select(web, Latitude, Longitude) %>% inner_join(data.symm, by="web")

head(data.symm)
dim(data.symm)

data.full <- meta %>% select(web, Latitude, Longitude) %>% left_join(data.full, by="web") %>% left_join(symm.cts, by="web") %>% mutate(lat=Latitude, lon=Longitude)

head(data.full)
dim(data.full)


# degree distributions! New as of 2020.03.20
degrees <- read.csv("output/degree_per_plant_phy.csv", h=TRUE) %>% filter(!is.na(symmetry))

degrees <- degrees %>% left_join(meta[,c("web", "Latitude")]) %>% mutate(lat.x = sqrt(abs(Latitude)), lat.sc = (lat.x-mean(lat.x))/sd(lat.x), pc1.sc = PC1/sd(PC1, na.rm=TRUE), pc2.sc = PC2/sd(PC2, na.rm=TRUE), pc3.sc = PC3/sd(PC3, na.rm=TRUE), pc4.sc = PC4/sd(PC4, na.rm=TRUE), rel.mn.n.shared=mn.n.shared/n.poll) %>% select(plant, web, Latitude, lat.x, lat.sc, symmetry, PC1, PC2, PC3, PC4, pc1.sc, pc2.sc, pc3.sc, pc4.sc, n.poll, n.shared, rel.mn.n.shared)

head(degrees)
dim(degrees)



# write out
write.table(data.full, "output/full-network_stats-meta.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(data.symm, "output/sub-network_stats-meta.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(degrees, "output/plant_stats-meta.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# read back in
data.full <- read.csv("output/full-network_stats-meta.txt", h=TRUE, sep="\t")
data.symm <- read.csv("output/sub-network_stats-meta.txt", sep="\t", h=TRUE)
degrees <- read.csv("output/plant_stats-meta.txt", sep="\t", h=TRUE)

deg.mn <- degrees %>% group_by(web, symmetry) %>% summarize(mnPoll = mean(n.poll), mdPoll = median(n.poll), N=length(n.poll))

mt5 <- intersect(filter(deg.mn, symmetry=="actinomorphic", N>=5)$web, filter(deg.mn, symmetry=="zygomorphic", N>=5)$web)

deg.mt5 <- filter(degrees, web %in% mt5, !is.na(PC1))
data.symm.mt5 <- filter(data.symm, matrix %in% mt5)

dim(data.symm.mt5) # 78, meaning 39 ..?



