# Figures and descriptive analysis
# from code via CJC
# assumes local environment
# last used/modified jby, 2020.06.26

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

data.full <- read.csv("output/full-network_stats-meta.txt", h=TRUE, sep="\t")
data.symm <- read.csv("output/sub-network_stats-meta.txt", sep="\t", h=TRUE)
degrees <- read.csv("output/plant_stats-meta.txt", sep="\t", h=TRUE)

# some data-munging
deg.mn <- degrees %>% group_by(web, symmetry) %>% summarize(mnPoll = mean(n.poll), mdPoll = median(n.poll), N=length(n.poll), mnMNshared=mean(rel.mn.n.shared), mnMNsharedA=mean(rel.mn.n.shared.a), mnMNsharedZ=mean(rel.mn.n.shared.z))

data.symm <- data.symm %>% left_join(deg.mn[,c("web", "symmetry", "mnPoll", "mnMNshared", "mnMNsharedA", "mnMNsharedZ")])

mt5 <- intersect(filter(deg.mn, symmetry=="Actinomorphic", N>=5)$web, filter(deg.mn, symmetry=="Zygomorphic", N>=5)$web)

deg.mt5 <- filter(degrees, web %in% mt5, !is.na(PC1))
data.symm.mt5 <- filter(data.symm, web %in% mt5)

dim(data.symm.mt5) # 78, meaning 39 ..?
dim(deg.mt5) # 2308

#--------------------------------------------------------------------
# color scheme

colors <- c("darkseagreen3", "sienna1", "royalblue3") # whole-matrix; zygomorphic; actinomorphic

#--------------------------------------------------------------------
# FIGURE 1

# buncha processing and transformation ...
world <- ne_download(type="land", category="physical", returnclass="sf")
class(world)
data.full.spatial <- data.full %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
data.full.spatial$Longitude = data.full$Longitude
data.full.spatial$Latitude = data.full$Latitude

# map of locations
map = ggplot(filter(world,min_zoom<2)) + geom_sf(color=NA, fill="gray70") + geom_sf(data=data.full.spatial, aes(size=Ntot), alpha=0.75, fill=colors[1], color="black", pch=21) + coord_sf(crs = "+proj=moll") + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(), panel.border=element_blank(), panel.grid = element_line(colour = 'transparent'), legend.position="none")

map

# stacked histograms of sample sizes
Ngtr <- data.full[,c("web", "Ntot", "Nact", "Nzyg", "Nnon")] %>% gather("Counting", "N", 2:5) %>% mutate(Counting = factor(Counting, c("Ntot", "Nact", "Nzyg", "Nnon")), web=factor(web, data.full$web[order(data.full$Ntot)]))

Ngtr2 <- data.full[,c("web", "Ntot", "Nact", "Nzyg", "Nnon")] %>% gather("Counting", "N", 2:5) %>% mutate(Counting = factor(Counting, c("Ntot", "Nact", "Nzyg", "Nnon")), web=factor(web, data.full$web[order(data.full$Ntot, decreasing=TRUE)]))

ss1 = ggplot(filter(Ngtr, Counting!="Ntot"), aes(y=N, x=web, fill=Counting)) + geom_col() + scale_fill_manual(values=colors[3:1], labels=c("Actinomorphic", "Zygomorphic", "n. det."), name="Floral symmetry") + coord_flip() + scale_x_discrete(labels=NULL) + labs(x="Plant-visitor network", y="Plant species") + theme_ipsum(axis_text_size=10, axis_title_size=14, plot_margin = margin(10, 10, 15, 10)) + theme(legend.position="none", panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(), axis.ticks=element_blank(), panel.border=element_blank())

ss2 = ggplot(filter(Ngtr, Counting!="Ntot", web!="Robertson-1928"), aes(y=N, x=web, fill=Counting)) + geom_col() + scale_fill_manual(values=colors[3:1], labels=c("Actinomorphic", "Zygomorphic", "n. det."), name="Floral symmetry") + coord_flip() + scale_x_discrete(labels=NULL) + labs(x=NULL, y="Plant species") + theme_ipsum(axis_text_size=10, axis_title_size=14, plot_margin = margin(10, 10, 15, 10)) + theme(legend.position=c(0.7, 0.2), panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(), legend.key.size=unit(0.3, "cm"), , legend.text = element_text(size = 10))

ss1
ss2

npoll_comp <- ggplot(degrees, aes(x=symmetry, y=n.poll, fill=symmetry, color=symmetry, group=symmetry)) + geom_point(position=position_jitter(), alpha=0.05) + geom_boxplot(alpha=0.75, color="black", width=0.5, outlier.alpha=0) + scale_y_log10() + scale_fill_manual(values=colors[3:2]) + scale_color_manual(values=colors[3:2]) + scale_x_discrete(labels=c("Act.", "Zyg.")) + labs(x="Floral symmetry", y="Visitor species") + theme_ipsum(base_size=14, axis_text_size=11, axis_title_size=14, caption_size=11, plot_margin = margin(5, 5, 10, 5), axis_title_just = "ct") + theme(legend.position="none")

npoll_comp


# new multi panel format for more dimensions ...
Rshare <- degrees %>% dplyr::select(symmetry, rel.mn.n.shared, rel.mn.n.shared.a, rel.mn.n.shared.z) %>% rename(`All~species`=rel.mn.n.shared, `Act.~species`=rel.mn.n.shared.a, `Zyg.~species`=rel.mn.n.shared.z) %>% gather("Sharing.with", "Sharing", 2:4) %>% mutate(Sharing.with=factor(Sharing.with, c("Zyg.~species", "Act.~species", "All~species")))

head(Rshare)

Rshare_multi <- ggplot(Rshare, aes(x=Sharing.with,symmetry, y=Sharing, fill=symmetry, color=symmetry, group=Sharing.with,symmetry)) + geom_boxplot(alpha=0.75, color="black", width=0.5, outlier.alpha=0) + coord_flip() + facet_wrap("symmetry", labeller=label_parsed, nrow=2) + scale_fill_manual(values=colors[3:2]) + scale_color_manual(values=colors[3:2]) + scale_x_discrete(labels=c("Zyg. species", "Act. species", "All species")) + labs(x="Sharing visitor species with", title="For species whose flowers are", y="Prop. shared visitor species") + theme_ipsum(base_size=14, axis_text_size=11, axis_title_size=13, strip_text_size=12, plot_margin = margin(7, 5, 10, 5), plot_title_size=14, plot_title_margin=5, axis_title_just = "ct") + theme(legend.position="none", panel.spacing.y=unit(0.25, "lines"))

Rshare_multi

npoll_lat <- ggplot() + geom_smooth(data=data.symm, aes(y=mnPoll, x=abs(Latitude), color=symmetry, lty=symmetry), method="lm", se=FALSE) + geom_line(data=data.symm, aes(y=mnPoll, x=abs(Latitude), group=web, color=symmetry), alpha=0.6, color="gray30", size=0.75) + geom_point(data=data.symm, aes(y=mnPoll, x=abs(Latitude), group=web, color=symmetry, shape=symmetry), alpha=0.75, size=4) + scale_color_manual(values=colors[3:2], name="Floral symmetry", labels=c("Actinomorphic", "Zygomorphic")) + scale_shape_manual(values=c(20,18), name="Floral symmetry", labels=c("Actinomorphic", "Zygomorphic")) + scale_linetype_manual(values=1:2, guide=FALSE) + scale_y_continuous(trans="log10") + theme_ipsum(axis_text_size=14, axis_title_size=14, plot_margin=margin(10,15,15,15)) + theme(legend.position=c(0.3,0.8), legend.text = element_text(size=11)) + labs(y="Mean visitor species per plant", x="Latitude (°N or S)") 

npoll_lat


# final multipanel output
{cairo_pdf("output/figures/Fig01_map-symm-npoll-Rshare.pdf", width=7, height=10)
ggdraw() + draw_plot(map, -0.07, 0.675, 0.6, 0.35) + draw_plot(ss1, 0.5, 0.7, 0.25, 0.275) + draw_plot(ss2, 0.725, 0.7, 0.25, 0.275) + draw_plot(npoll_comp, 0.02, 0.35, 0.35, 0.325) + draw_plot(Rshare_multi, 0.4, 0.35, 0.55, 0.35) + draw_plot(npoll_lat, 0.2, 0, 0.6, 0.35) + draw_plot_label(c("A", "B", "C", "D", "E"), x=c(0, 0.5, 0, 0.4, 0.18), y=c(1, 1, 0.7, 0.7, 0.35)) + 
draw_line(x=c(0.81,0.81), y=c(0.555,0.586), color="gray60", lwd=1) + draw_label("***", x=0.805, y=0.565, hjust=1, vjust=0.5) + # A->Z vs A->A
draw_line(x=c(0.872,0.872), y=c(0.42,0.451), color="gray60", lwd=1) + draw_label("**", x=0.867, y=0.435, hjust=1, vjust=0.5) + # Z->Z vs Z->A
draw_line(x=c(0.89,0.89), y=c(0.42,0.556), color="gray60", lwd=1) + draw_label("***", x=0.882, y=0.485, hjust=1, vjust=0.5) + # Z->Z vs A->Z
draw_line(x=c(0.91,0.91), y=c(0.481,0.616), color="gray60", lwd=1) + draw_label("*", x=0.905, y=0.54, hjust=1, vjust=0.5) + # Z->all vs A->all
draw_line(x=c(0.17,0.29), y=c(0.62,0.62), color="gray60", lwd=1) + draw_label("*", x=0.23, y=0.62, hjust=0.5, vjust=0) # visitor counts

}
dev.off()


#--------------------------------------------------------------------
# FIGURE 2

# test for differences in median ----------------
data.symm %>% group_by(symmetry) %>% summarize(medPollShare = median(mnMNshared), medNpoll = median(mnPoll))

wilcox.test(mnPoll~symmetry, data=data.symm.mt5, paired=TRUE, alt="g") # A > Z, p = 0.0002
wilcox.test(mnMNshared~symmetry, data=data.symm.mt5, paired=TRUE, alt="l") # A > Z, p = 0.06

wilcox.test(connectance~symmetry, data=data.symm.mt5, paired=TRUE, alt="l") # Z > A, p = 7.9e-05
wilcox.test(web.asymmetry~symmetry, data=data.symm.mt5, paired=TRUE, alt="l") # Z > A, p = 0.0007
wilcox.test(links.per.species~symmetry, data=data.symm.mt5, paired=TRUE, alt="g") # A > Z, p = 3.1e-10
wilcox.test(nestedness~symmetry, data=data.symm.mt5, paired=TRUE, alt="l") # Z > A, 3.1e-07
wilcox.test(log(host.sharing)~symmetry, data=data.symm.mt5, paired=TRUE) # n.s., p = 0.36
wilcox.test(modularity~symmetry, data=data.symm.mt5, paired=TRUE, alt="l") # Z > A, p = 0.0004

wilcox.test(coex.plant~symmetry, data=data.symm.mt5, paired=TRUE, alt="g") # A > Z, p = 3.7e-09
wilcox.test(coex.poll~symmetry, data=data.symm.mt5, paired=TRUE, alt="g") # A > Z, p = 3.1e-08

pairdat <- data.symm.mt5 %>% mutate(mnPoll.x = log10(mnPoll)) %>% dplyr::select(web, symmetry, mnPoll.x, mnMNshared, connectance, web.asymmetry, coex.plant, coex.poll) %>% rename(`log[10](Visitor~species)`=mnPoll.x, `Visitor~species~sharing`=mnMNshared, `Connectance`=connectance, `Web~asymmetry`=web.asymmetry, `Plant~coex.~robustness`=coex.plant, `Visitor~coex.~robustness`=coex.poll) %>% gather("var","val",3:8)

pairdat$var <- factor(pairdat$var, c("log[10](Visitor~species)", "Visitor~species~sharing", "Connectance", "Web~asymmetry", "Plant~coex.~robustness", "Visitor~coex.~robustness"))

siglabs <- data.frame(var = factor(c("log[10](Visitor~species)", "Visitor~species~sharing", "Connectance", "Web~asymmetry", "Plant~coex.~robustness", "Visitor~coex.~robustness")), lab=c("*", NA, "**", "*", "***", "***"), x=1.5, y=c(2, NA, 0.75, 0.9, 25, 9))

pairmeds <- pairdat %>% group_by(symmetry, var) %>% summarize(med=median(val))

pairplot <- ggplot() + geom_line(data=pairdat, aes(y=val, x=symmetry, group=web), alpha=0.4, color="gray30", size=0.75) + geom_point(data=pairdat, aes(y=val, x=symmetry, group=web, color=symmetry, shape=symmetry), alpha=0.65, size=4) + geom_text(data=siglabs, aes(x=x, y=y, label=lab), size=5) + geom_point(data=pairmeds, aes(x=symmetry, y=med), pch="—", color="black", size=6) + facet_wrap("var", scale="free", labeller=label_parsed, nrow=4) + scale_color_manual(values=colors[3:2]) + scale_shape_manual(values=c(20,18)) + scale_x_discrete(labels=c("Act.", "Zyg.")) + labs(y=NULL, x="Floral symmetry", caption=expression(paste("*")~Paired~Wilcoxon~test~p<0.001~paste("; **")~p<10^-4~paste("; ***")~p<10^-5)) + theme_ipsum(axis_text_size=10, axis_title_size=13, axis_title_just="ct", strip_text_size=11, plot_margin=margin(10,15,15,15)) + theme(legend.position="none", panel.spacing.y=unit(0.5, "line"), panel.spacing.x=unit(0.75, "line"))

pairplot

# correlations with PropZ ----------------
data.z <- filter(data.full, PropZ>0)

descdist(data.z$PropZ)
descdist(log10(data.z$PropZ)) # boom, normalized
descdist(log10(data.z$mnRshared[!is.na(data.z$mnRshared)])) # noice
descdist(log10(data.z$connectance)) # noice
descdist(data.z$web.asymmetry)
descdist(data.z$nestedness)
descdist(data.z$modularity)
descdist(log10(data.z$coex.plant))
descdist(data.z$coex.poll)

cor.test(~log10(mnPoll)+PropZ, data=data.z) # cor = -0.21, p = 0.05
cor.test(~mnRshared+PropZ, data=data.z) # cor = 0.39, p = 0.0002
cor.test(~connectance+PropZ, data=data.z) # cor = 0.46, p = 4.6e-06
cor.test(~web.asymmetry+PropZ, data=data.z) # cor = -0.30, p = 0.004
cor.test(~nestedness+PropZ, data=data.z) # cor = 0.23, p = 0.03
cor.test(~modularity+PropZ, data=data.z) # n.s.
cor.test(~coex.plant+PropZ, data=data.z) # n.s.
cor.test(~coex.poll+PropZ, data=data.z) # cor = 0.32, p = 0.003

cor.test(~log10(mnPoll)+log10(PropZ), data=data.z) # n.s.
cor.test(~log10(mnRshared)+log10(PropZ), data=data.z) # cor = 0.22, p = 0.04
cor.test(~log10(connectance)+log10(PropZ), data=data.z) # cor = 0.21, p = 0.04
cor.test(~web.asymmetry+log10(PropZ), data=data.z) # cor = -0.26, p = 0.01
cor.test(~nestedness+log10(PropZ), data=data.z) # n.s.
cor.test(~modularity+log10(PropZ), data=data.z) # n.s.
cor.test(~coex.plant+log10(PropZ), data=data.z) # n.s.
cor.test(~coex.poll+log10(PropZ), data=data.z) # cor = 0.27, p = 0.01


cordat <- data.z %>% mutate(PropZ.x = log10(PropZ), mnPoll.x = log10(mnPoll), mRs.x = log10(mnRshared), conn.x = log10(connectance)) %>% dplyr::select(web, PropZ, mnPoll.x, mRs.x, conn.x, web.asymmetry, coex.plant, coex.poll) %>% rename(`log[10](Visitor~species)`=mnPoll.x, `log[10](Visitor~species~sharing)`=mRs.x, `log[10](Connectance)`=conn.x, `Web~asymmetry`=web.asymmetry, `Plant~coex.~robustness`=coex.plant, `Visitor~coex.~robustness`=coex.poll) %>% gather("var","val",3:8)

cordat$var <- factor(cordat$var, c("log[10](Visitor~species)", "log[10](Visitor~species~sharing)", "log[10](Connectance)", "Web~asymmetry", "Plant~coex.~robustness", "Visitor~coex.~robustness" ))

sigcors <- c("log[10](Visitor~species~sharing)", "log[10](Connectance)", "Web~asymmetry", "Visitor~coex.~robustness")

#corlabs <- data.frame(var=factor(levels(cordat$var), levels(cordat$var)), cor=c(NA, "cor = 0.22\np = 0.04", "cor = 0.21\np = 0.04", "cor = -0.26\np = 0.01", NA, "cor = 0.27\np = 0.01"), y = c(NA, 0, 0, -0.5, NA, 12), x = c(NA, -1.25, -1.25, -1.25, NA, -1.25))

corplots <- ggplot() + geom_smooth(data=filter(cordat, var%in%sigcors), aes(x=PropZ, y=val), method="lm", color="white", fill="gray80") + geom_point(data=cordat, aes(x=PropZ, y=val), color=colors[1], alpha=0.75) + facet_wrap("var", scale="free_y", labeller=label_parsed, nrow=4) + theme_ipsum(axis_text_size=10, axis_title_size=13, axis_title_just="ct", strip_text_size=12, plot_margin=margin(10,15,15,15)) + labs(y=NULL, x="Proportion zygomorphic") + theme(legend.position="none", panel.spacing.y=unit(0.5, "line"), panel.spacing.x=unit(0.5, "line"))

# figure 02
{cairo_pdf("output/figures/Fig02_pairs-cors.pdf", width=9, height=7)

ggdraw() + draw_plot(pairplot, 0, 0, 0.45, 1) + draw_plot(corplots, 0.45, 0, 0.55, 1) + draw_plot_label(c("A", "B"), x=c(0, 0.45), y=1) 

}
dev.off()



