#!/usr/bin/env Rscript
# Loading libraries and datafiles
library('MASS'); library('Peacock.test'); 
pp_2LAO<-as.matrix(read.table('pp_2LAO.dat'));pp_2LAO<-pp_2LAO[,-1]
pp_modelo<-read.table('pp_modelo.dat'); pp_modelo<-pp_modelo[,-1]
compare_dihedrals_mutants<-function(mut1,mut2){
	vec<-vector()	
	for (i in seq(1,length(pp_2LAO[1,]),by=2)){	
		m<-cbind(mut1[,i],mut1[,(i+1)])
		n<-cbind(mut2[,i],mut2[,(i+1)])	
		p<-peacock2(m,n)	
		vec<-c(vec,p)
		} 
return(vec)
}
# Diferencias
di_k<-compare_dihedrals_mutants(pp_2LAO,pp_modelo)
write(di_k, file = "comp_2LAO_modelo.dat",sep = " ",ncolumns=1)
