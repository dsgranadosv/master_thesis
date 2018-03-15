#!/usr/bin/env Rscript
#################################################################
#  MULTIPLE FIGURE OF RESIDUES				    #####
#  WITH DIFFERENCES IN THEIR RAMACHANDRAN PLOTS		    #####
#################################################################
# Libraries
library('plyr')
# Reading of data file and adding one column 
pp_2LAO<-as.matrix(read.table('../pp_2LAO.dat'));pp_2LAO<-pp_2LAO[,-1]
pp_L117K<-as.matrix(read.table('../pp_L117K.dat'));pp_L117K<-pp_L117K[,-1] 
p_matrix<-as.matrix(read.table('../m_pvalues.tmp'))
index<-seq(2,237,by=1)
p_matrix<-cbind(index,p_matrix)
#### Definition of functions ####
ramachandran<-function(var1,var2,color,main_title,agg){
par(mar=c(5,4.3,4.1,2))
plot(x=var1,y=var2,pch=19,col=color,xlab='Phi Angle',ylab='Psi Angle',xlim=c(-180,180),ylim=c(-180,180),main=main_title,cex.main=2,cex.lab=1.8,cex.axis=1.6)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
}
plot_the_different<-function(mat,pp_ref,pp_comp){
	inds<-mat[,1]
	phi<-sapply(inds,function(x) ((x*2)-3))
	psi<-sapply(phi,function(x) x+1)
	pp_inds<-cbind(phi,psi)
	renglones<-(length(pp_inds[,1]))/5
	renglones<-round_any(renglones,1,f=ceiling)
	tiff("../dif_rama.tiff", height = 18, width = 36.132, units = 'cm',compression = "lzw", res = 300)
	par(mfrow=c(renglones,5))
	par(oma=c(3.5,3,5.2,0))
	for (i in 1:length(pp_inds[,1])){	
		ramachandran(pp_comp[,(pp_inds[i,1])],pp_comp[,(pp_inds[i,2])],rgb(0,0,1,0.3),paste(mat[i,1]))
		par(new=TRUE)
		ramachandran(pp_ref[,(pp_inds[i,1])],pp_ref[,(pp_inds[i,2])],rgb(1,0,0,0.3),paste(mat[i,1]))
	}
}
# Conditional print of some rows of the p-value matrix based in multiple conditions
dif_only_L117K<-p_matrix[p_matrix[,2] > 0.5 & p_matrix[,3] < 0.5 & p_matrix[,4] < 0.5 , ]
############# Calling the device and additional information to the original plot ##################
plot_the_different(dif_only_L117K,pp_2LAO,pp_L117K)
# Common title and legend 
mtext("Residuos con diferencias significativas en los gráficos de \nRamachandran exclusivos de L117K", outer = TRUE, cex = 2,side=3,line=0)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("LAO","L117K"),col=c(rgb(1,0,0,0.6),rgb(0,0,1,0.6)),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,15), cex = 2.5)
dev.off()
