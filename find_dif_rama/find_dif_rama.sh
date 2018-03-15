#!/bin/bash
#####################################################################
####		FINDING RESIDUES WITH DIFFERENCES  	        #####
####		IN RAMACHANDRAN PLOTS	                      	#####
####						              	#####
####  Description: This shell script runs a series of 		#####
####  Rscripts. First, it performs a KS-Test of the bidimen-	#####
####  sional distribution associated to the ramachandran plots 	#####
####  of the mutants'residues. It takes 2LAO distributions as	#####
####  references. Later it finds the statiscally different	#####
####  distributions and plot them in a multiple figure 		#####
####  Author: Diego Granados					#####	
####  								#####
#####################################################################
##### Running the KS-TEST for every residue for all the mutants #####
for i in L117K L117R L117Q
do
sed "s/modelo/$i/g" Rscripts/ks_test.R > "$i".tmp
done
parallel Rscript ::: *.tmp
#### Paste for constructing a matrix of pvalues for every residue for all the mutants ####
paste comp_2LAO_L117K.dat comp_2LAO_L117R.dat comp_2LAO_L117Q.dat > m_pvalues.tmp
#### Generate the multiple plot figure ####
Rscript Rscripts/construct_figure.R 
###########################
# Cleaning temporal files #
###########################
rm *.tmp 
