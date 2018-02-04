#!/usr/bin/env Rscript
###################################################################
###################################################################
#######							   ########	
#######	            HISTOGRAM OF BACKBONE H-BONDS          ########
#######						           ########
###################################################################
###################################################################
# Description: This script graphs a multiple plot figure composed #
# by three histograms of the counts of backbone Hbonds in three	  #
# different point mutants. A distance value between Histograms is #
# calculated to quantify their differences			  #
# Note: This script requires the library "Histogram Tools" 	  #
# if you don't have it, you can always install it by typing 	  #
# <install.packages("HistogramTools", dependencies=TRUE)>	  #
# Author: Diego Granados					  #
# Date: January 29, 2018					  #
# R-Version: 3.2.3						  #	
###################################################################
# Loading libraries and reading datafiles
library('HistogramTools')
LQ_A<-read.table('total_bb_L117Q_lobule_A.dat')
LR_A<-read.table('total_bb_L117R_lobule_A.dat')
LK_A<-read.table('total_bb_L117K_lobule_A.dat')
LA_A<-read.table('total_bb_2LAO_lobule_A.dat')
# Definition of the breaks for the histograms, change this variable 
#  for representing your data correctly
vec<-seq(10,60,by=2)
# Define histograms that we are comparing
HLA<-hist(LA_A$V2,breaks=vec,plot=FALSE)
HLK<-hist(LK_A$V2,breaks=vec,plot=FALSE)
HLR<-hist(LR_A$V2,breaks=vec,plot=FALSE)
HLQ<-hist(LQ_A$V2,breaks=vec,plot=FALSE)
# Distance between histograms, we always use the same reference
DLA_LK<-1-intersect.dist(HLA,HLK)
DLA_LR<-1-intersect.dist(HLA,HLR)
DLA_LQ<-1-intersect.dist(HLA,HLQ)
############### PLOTTING THE FIGURE ######################
# setting up a high resolution device
tiff("FileName.tiff", height = 12, width = 28, units = 'cm',compression = "lzw", res = 300)
#Layout for the figure
par(mfrow=c(1,3))
# Expand margins, internal and outer 
par(mar=c(2.5,2.5,0.7,1.0))
par(oma=c(4.5,2.9,4.0,0))
#### Plotting the histograms ####
hist(LA_A$V2,col=rgb(1,0,0,0.4),breaks=vec,xlab='Freq',main='LAO VS L117K',cex.axis=1.5)
hist(LK_A$V2,col=rgb(0,0,1,0.4),add=T,breaks=vec) 
# Add value of intersection of histograms #
eq<-bquote(bold("H =" ~ .(DLA_LK)))
text(50,7500,eq,cex=1.7,col='red')
hist(LA_A$V2,col=rgb(1,0,0,0.4),breaks=vec,main='LAO VS L117R',,cex.axis=1.5)
hist(LR_A$V2,col=rgb(0,0,1,0.4),add=T,breaks=vec)
# Add value of intersection of histograms #
eq<-bquote(bold("H =" ~ .(DLA_LR)))
text(50,7500,eq,cex=1.7,col='red')
hist(LA_A$V2,col=rgb(1,0,0,0.4),breaks=vec,main='LAO VS L117Q',cex.axis=1.5)
hist(LQ_A$V2,col=rgb(0,0,1,0.4),add=T,breaks=vec)
# Add value of intersection of histograms #
eq<-bquote(bold("H =" ~ .(DLA_LQ)))
text(50,7500,eq,cex=1.7,col='red')
#Comon title for the multiple histogram figure
mtext("Cantidad total de puentes de hidrógeno en la cadena principal (lóbulo A)", outer = TRUE, cex = 1.5,side=3,line=2)
#Shared axes labels
mtext("Frecuencia", outer = TRUE, cex = 1.2,side=2,line=1,las=0)
mtext("Número de puentes de hidrógeno de cadena principal", outer = TRUE, cex = 1.2,side=1,line=1,las=1)
#Shared legend 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("LAO","Mutante","Intersección de histogramas"),col=c(rgb(1,0,0,0.4),rgb(0,0,1,0.4),rgb(0.5,0,1,0.8)),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,15,15), cex = 1.5)
dev.off()
