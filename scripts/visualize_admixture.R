#!/usr/bin/env Rscript
# This script is for mostly explorative purposes for some reason the clumping does not work perfectly
library(pophelper)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(data.table)
library("optparse")

option_list = list(
  make_option(c("-i", "--inputdir"), type="character", default="results/",
              help="input directory with admixture files", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="tishkoff_sgdp_merged_ld_pruned",
              help="prefix of admixture files", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default="data/tishkoff_sgdp_merged_ld_pruned.fam",
              help="metadata file with individual pop mapping", metavar="character"),
  make_option(c("-c", "--popcoords"), type="character", default="data/individual_geo_coords.tab",
              help="File with population coordinates pop, latitude, longitude", metavar='character'),
  make_option(c("-o", "--out"), type="character", default="results/admixture_plot",
              help="output file name prefix", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

fam <- fread(file=opt$metadata, select=c(1, 2))
names(fam) <- c("FID", "IID")

coords <- fread(file=opt$popcoords)
names(coords) <- c("FID", "IID", "Latitude", "Longitude")
merged <- merge(fam, coords, by.x='IID', by.y="IID")
merged_sorted <- merged[order(merged$Latitude, merged$Longitude, decreasing=TRUE),]
merged_sorted <- merged_sorted[!duplicated(merged_sorted$FID.x),]
names(merged_sorted) <- c("IID", "pop", "pop1", "Latitude", "Longitude")

subsetgrp <- merged_sorted$pop
subsetgrp[subsetgrp == "BantuTswana"] <- "Bantu Tswana"
subsetgrp[subsetgrp == "BantuKenya"] <- "Bantu Kenya"
subsetgrp[subsetgrp == "Yemenite_Jew"] <- "Yemenite Jewish"
subsetgrp[subsetgrp == "Iraqi_Jew"] <- "Iraqi Jewish"
subsetgrp[subsetgrp == "Sabue"] <- "Chabu"
subsetgrp[subsetgrp == "Khomani_San"] <- "Khomani San"

fam <- fread(file=opt$metadata, select=c(1))
names(fam) <- c("K=16")
fam[fam == "BantuTswana"] <- "Bantu Tswana"
fam[fam == "BantuKenya"] <- "Bantu Kenya"
fam[fam == "Yemenite_Jew"] <- "Yemenite Jewish"
fam[fam == "Iraqi_Jew"] <- "Iraqi Jewish"
fam[fam == "Sabue"] <- "Chabu"
fam[fam == "Khomani_San"] <- "Khomani San"

sfiles <- list.files(path=opt$inputdir, pattern=paste(opt$prefix, ".*.Q", sep=""), full.names=T)
qlist <- sortQ(readQ(files=sfiles))
# qlist <- qlist[c(12, 13)]
# print((qlist[[1]]))
# aligned_qlist <- alignK(qlist)

colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2",
          "#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368",
          "#A4A4A4","#610B5E")


# # coords <- coords[!duplicated(coords$Latitude),]
plotQ(alignK(qlist),imgoutput="join", grplab=fam, grplabangle=90, width=15, grplabjust=1,
      grplabspacer=0, grplabpos=0.95,  ordergrp=F, linealpha=0, showdiv=T, pointalpha=0, showgrplab=T, panelratio=c(1,1),
      splab=paste0("K=",sapply(qlist,ncol)), clustercol=colors, sortind='label', subsetgrp=subsetgrp,
      indlabheight=-.03, splabangle=0, imgtype="pdf", outputfilename=opt$out, exportpath=getwd(),showticks=TRUE,
      theme='theme_bw', indlabhjust=1,)
