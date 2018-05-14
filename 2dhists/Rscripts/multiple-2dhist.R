#Reading libraries
library('MASS')
library('RColorBrewer')
library('sinkr')
###### Reading data files and converting them to matrices  ########
LAO<-as.matrix(read.table('../data/at-angl_2LAO_open_all-dis.dat'))
L117K<-as.matrix(read.table('../data/at-angl_L117K_open_all-dis.dat'))
L117Q<-as.matrix(read.table('../data/at-angl_L117Q_open_all-dis.dat'))
L117R<-as.matrix(read.table('../data/at-angl_L117R_open_all-dis.dat'))
# Function for building a 2DHist
make_2dhist<-function(mat){
	a<-kde2d(mat[,1],mat[,2],n = 20,lims=c(90,125,40,95))
return(a)
}
#Apply function 
histograms<-sapply(list(LAO,L117K,L117Q,L117R),make_2dhist)
#Color palette for plotting
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(20)
##### Generating the multiple plot figure #######
# Device for saving the figure
tiff("2dhists_and_boxplots.tiff", height = 19.3, width = 38, units = 'cm',compression = "lzw", res = 300)
# Plot array
patron<-c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,
	  1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,
          1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,
          1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,
          6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,0,
          6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,0, 
          6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,0,
          6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,0)
layout(matrix(patron, 8, 17, byrow = TRUE))
# Custom breaks for data
b=seq(0,0.0055,by=0.0055/length(r))
#Margins for the superior panels
par(mar=c(2,3.5,4,0.5))
par(oma=c(2,5.2,5,4))
# Adding 2dhists to the plot
# The axis function is for controlling the position and labels 
at_x=seq(90,130,by=10)
at_y=seq(40,95,by=10)
image(histograms[,1],breaks=b,col=r,cex.main=3,main='LAO',axes=F)
axis(1,lwd.ticks=2,at=at_x,labels=T,cex.axis=2)
axis(2,lwd.ticks=2,at=at_y,labels=T,cex.axis=2)
image(histograms[,2],breaks=b,col=r,cex.main=3,main='L117K',axes=F)
axis(1,lwd.ticks=2,at=at_x,labels=T,cex.axis=2)
axis(2,lwd.ticks=2,at=at_y,labels=F)
image(histograms[,4],breaks=b,col=r,cex.main=3,cex.axis=2,main='L117R',axes=F)
axis(1,lwd.ticks=2,at=at_x,labels=T,cex.axis=2)
axis(2,lwd.ticks=2,at=at_y,labels=F)
image(histograms[,3],breaks=b,col=r,cex.main=3,main='L117Q',axes=F)
axis(1,lwd.ticks=2,at=at_x,labels=T,cex.axis=2)
axis(2,lwd.ticks=2,at=at_y,labels=F)
# Color scale of densities for 2dhistograms
imageScale(histograms[,1]$z, col=r, axis.pos=4,cex.axis=1.8)
# Common title and common axis for the 2d histograms
title('Exploración conformacional de LAO y mutantes puntuales',outer=TRUE,cex.main=3.5)
mtext("Ángulo de torsión(°)",side=2,cex=1.8,outer=T,at=0.72,line=1)
mtext("Ángulo de apertura(°)",side=1,cex=1.8,outer=T,at=0.48,line=-23.5)
# Dataset definition for boxplots
a<-cbind(LAO[,1],L117K[,1],L117R[,1],L117Q[,1])
c<-cbind(LAO[,2],L117K[,2],L117R[,2],L117Q[,2])
# New margins
par(mar=c(4,4.5,6.5,1.0))
# The boxplots
boxplot(a,col='chocolate1',outline=FALSE,lwd=2.2,cex.axis=2,cex.lab=2,xlab='Mutantes',ylab='Ángulo de apertura',names=c("LAO","L117K","L117R","L117Q"))
# Title of boxplot #1
title("Distribuciones de apertura",line=1.5,cex.main=2)
boxplot(c,col='cornflowerblue',outline=FALSE,lwd=2.2,cex.axis=2,cex.lab=2,xlab='Mutantes',ylab='Ángulo de torsión',names=c("LAO","L117K","L117R","L117Q"))
# Title of boxplot #2
title("Distribuciones de torsión",line=1.5,cex.main=2)
dev.off()
