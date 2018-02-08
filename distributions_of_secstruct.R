#!/usr/bin/env Rscript
###################################################################
###################################################################
#######							   ########	
#######	            DISTRIBUTION OF SECCONDARY             ########
#######			   STRUCTURE			   ########
#######						           ########
###################################################################
###################################################################
# Description:  
# #
# #
# #
# #
# #
# #
# Author: Diego Granados					  #
# Date: February 05, 2018					  #
# R-Version: 3.2.3						  #	
###################################################################
# Loading datafiles and libraries
library('overlapping')
ss_2LAO<-read.table('secstruct_2LAO_all.dat')
ss_L117K<-read.table('secstruct_L117K_all.dat')
ss_L117R<-read.table('secstruct_L117R_all.dat')
ss_L117Q<-read.table('secstruct_L117Q_all.dat')
# This function count the percentage of a  particular asignation 
# of seccondary structure during the simulation. The three arguments 
# are the information of seccondary structure (data), 
# the stride value of the asignation we want to measure (asig)
# and the length of the protein (l) 
count_secstruct<-function(data,asig,l){
	vec <- vector()
		for (i in 1:length(data)){
		p<-(length(which(data[,i]==asig))/l)*100
		vec<-c(vec,p)
		}
	return(vec)
}
# We generate the probability densities for all the proteins
# This repetitive code must be improved with a series of apply function
# Alpha Helix

D_LA_H<-density(count_secstruct(ss_2LAO,asig="H",l=238))
D_LK_H<-density(count_secstruct(ss_L117K,asig="H",l=238))
D_LR_H<-density(count_secstruct(ss_L117R,asig="H",l=238))
D_LQ_H<-density(count_secstruct(ss_L117Q,asig="H",l=238))
#Beta Sheet
D_LA_E<-density(count_secstruct(ss_2LAO,asig="E",l=238))
D_LK_E<-density(count_secstruct(ss_L117K,asig="E",l=238))
D_LR_E<-density(count_secstruct(ss_L117R,asig="E",l=238))
D_LQ_E<-density(count_secstruct(ss_L117Q,asig="E",l=238))
# Values of distance between density plots #
# Alpha Helix
D_LA_LK_H<-overlap(list(count_secstruct(ss_2LAO,asig="H",l=238),count_secstruct(ss_L117K,asig="H",l=238)),nbins = 1000, plot = FALSE, partial.plot = FALSE)
D_LA_LR_H<-overlap(list(count_secstruct(ss_2LAO,asig="H",l=238),count_secstruct(ss_L117R,asig="H",l=238)),nbins = 1000, plot = FALSE, partial.plot =FALSE)
D_LA_LQ_H<-overlap(list(count_secstruct(ss_2LAO,asig="H",l=238),count_secstruct(ss_L117Q,asig="H",l=238)),nbins = 1000, plot = FALSE, partial.plot =FALSE)
# Beta sheet
D_LA_LK_E<-overlap(list(count_secstruct(ss_2LAO,asig="E",l=238),count_secstruct(ss_L117K,asig="E",l=238)),nbins = 1000, plot = FALSE, partial.plot = FALSE)
D_LA_LR_E<-overlap(list(count_secstruct(ss_2LAO,asig="E",l=238),count_secstruct(ss_L117R,asig="E",l=238)),nbins = 1000, plot = FALSE, partial.plot =FALSE)
D_LA_LQ_E<-overlap(list(count_secstruct(ss_2LAO,asig="E",l=238),count_secstruct(ss_L117Q,asig="E",l=238)),nbins = 1000, plot = FALSE, partial.plot =FALSE)
###########################################################
########### PLOTTING THE FIGURES ##########################
###########################################################
###########################################################
#### ALPHA HELIX #######
tiff("distro_alpha.tiff", height = 12, width = 28, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,3))
par(mar=c(2.5,2.5,1.9,1.5))
par(oma=c(4.5,2.9,4.0,0))
plot(D_LA_H,col=rgb(1,0,0,0.6),main='LAO VS L117K',cex.axis=1.6)
polygon(D_LA_H,col=rgb(1,0,0,0.6),lwd=2)
lines(D_LK_H,col=rgb(0,0,1,0.6))
polygon(D_LK_H,col=rgb(0,0,1,0.6),lwd=2)
eq<-bquote(bold("H =" ~ .(round(D_LA_LK_H$OV,2))))
text(26,0.2,eq,cex=1.7,col='red')
plot(D_LA_H,col=rgb(1,0,0,0.6),main='LAO VS L117R',cex.axis=1.6)
polygon(D_LA_H,col=rgb(1,0,0,0.6),lwd=2)
lines(D_LR_H,col=rgb(0,0,1,0.6))
polygon(D_LR_H,col=rgb(0,0,1,0.6),lwd=2)
eq<-bquote(bold("H =" ~ .(round(D_LA_LR_H$OV,2))))
text(26,0.2,eq,cex=1.7,col='red')
plot(D_LA_H,col=rgb(1,0,0,0.6),main='LAO VS L117Q',cex.axis=1.6)
polygon(D_LA_H,col=rgb(1,0,0,0.6),lwd=2)
lines(D_LQ_H,col=rgb(0,0,1,0.6))
polygon(D_LQ_H,col=rgb(0,0,1,0.6),lwd=2)
eq<-bquote(bold("H =" ~ .(round(D_LA_LQ_H$OV,2))))
text(26,0.2,eq,cex=1.7,col='red')
#Comon title
mtext("Distribuciones del porcentaje de residuos en conformación alfa", outer = TRUE, cex = 1.5,side=3,line=2)
#Shared axes labels
mtext("Densidad de probabilidad", outer = TRUE, cex = 1.2,side=2,line=1,las=0)
mtext("Porcentaje de residuos en conformación de hélice alfa", outer = TRUE, cex = 1.2,side=1,line=1,las=1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("LAO","Mutante"),col=c(rgb(1,0,0,0.6),rgb(0,0,1,0.6)),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,15), cex = 1.5)
dev.off()
#### ALPHA HELIX #######
tiff("distro_beta.tiff", height = 12, width = 28, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,3))
par(mar=c(2.5,2.5,1.9,1.5))
par(oma=c(4.5,2.9,4.0,0))
plot(D_LA_E,col=rgb(1,0,0,0.6),main='LAO VS L117K',cex.axis=1.6)
polygon(D_LA_E,col=rgb(1,0,0,0.6),lwd=2)
lines(D_LK_E,col=rgb(0,0,1,0.6))
polygon(D_LK_E,col=rgb(0,0,1,0.6),lwd=2)
eq<-bquote(bold("H =" ~ .(round(D_LA_LK_E$OV,2))))
text(18,0.15,eq,cex=1.7,col='red')
plot(D_LA_E,col=rgb(1,0,0,0.6),main='LAO VS L117R',cex.axis=1.6)
polygon(D_LA_E,col=rgb(1,0,0,0.6),lwd=2)
lines(D_LR_E,col=rgb(0,0,1,0.6))
polygon(D_LR_E,col=rgb(0,0,1,0.6),lwd=2)
eq<-bquote(bold("H =" ~ .(round(D_LA_LR_E$OV,2))))
text(18,0.15,eq,cex=1.7,col='red')
plot(D_LA_E,col=rgb(1,0,0,0.6),main='LAO VS L117Q',cex.axis=1.6)
polygon(D_LA_E,col=rgb(1,0,0,0.6),lwd=2)
lines(D_LQ_E,col=rgb(0,0,1,0.6))
polygon(D_LQ_E,col=rgb(0,0,1,0.6),lwd=2)
eq<-bquote(bold("H =" ~ .(round(D_LA_LQ_E$OV,2))))
text(18,0.15,eq,cex=1.7,col='red')
#Comon title
mtext("Distribuciones del porcentaje de residuos en conformación alfa", outer = TRUE, cex = 1.5,side=3,line=2)
#Shared axes labels
mtext("Densidad de probabilidad", outer = TRUE, cex = 1.2,side=2,line=1,las=0)
mtext("Porcentaje de residuos en conformación de hoja beta", outer = TRUE, cex = 1.2,side=1,line=1,las=1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("LAO","Mutante"),col=c(rgb(1,0,0,0.6),rgb(0,0,1,0.6)),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,15), cex = 1.5)
dev.off()




