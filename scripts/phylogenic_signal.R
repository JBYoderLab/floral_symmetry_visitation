# Phylogenetic signal testing
# assumes local environment
# last used/modified jby, 2020.08.06

# setwd('~/Documents/Academic/Active_projects/pollination_networks')

library(bipartite)
library(codependent)
library(igraph)
library(reshape)
library(stringr)
library(treeplyr)
library(pez)
library(phytools)
library(nlme)
library(Taxonstand)

library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)

library(tidyverse)

source("../shared/Rscripts/base.R")
source("../shared/Rscripts/base_graphics.R")


#--------------------------------------------------------------------
# load and parse data

deg.phy <- read.csv("output/degree_per_plant_phy.csv", h=TRUE) # characters to map

deg.tree <- read.tree("output/working_tree.tre")

deg.phy.mn <- deg.phy %>% group_by(tip_name, symmetry) %>% summarize(mnDeg = mean(n.poll)) %>% filter(tip_name %in% deg.tree$tip.label) %>% as.data.frame(.)
filter(deg.phy.mn, duplicated(tip_name)) # whew

deg.phy.mn <- deg.phy.mn %>% as.data.frame(.) %>% mutate(symmetry = as.numeric(as.factor(symmetry))-1) # 0 = act, 1 = zyg
rownames(deg.phy.mn) <- deg.phy.mn$tip_name
deg.phy.mn <- deg.phy.mn %>% select(-tip_name)

head(deg.phy.mn)
dim(deg.phy.mn)

deg.p4d <- phylo4d(deg.tree, deg.phy.mn)
cols <- rep("blue", 2507)
cols[deg.phy.mn$symmetry=="Zygomorphic"] <- "orange"

#--------------------------------------------------------------------
# Phylogenetic signal tests

# phylogenetic signal ...
phyloSignal(p4d=deg.p4d, methods=c("all")) # this takes awhile




