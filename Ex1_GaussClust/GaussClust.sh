#!/bin/sh

##########################################################################################
#                             GaussClust v0.1.0, January 2016                              #
#   SHELL SCRIPT FOR AUTOMATING GAUSSIAN MIXTURE MODELING ON CONTINUOUS, DISCRETE, OR    #
#   COMBINED GENETIC/MORPHOLOGICAL/ECOLOGICAL DATA, FOR SPECIES DELIMITATION AND         #
#   CLASSIFICATION                                                                       #
#   Copyright (c)2016 Justin C. Bagley, Universidade de Brasília, Brasília, DF, Brazil.  #
#   See the README and license files on GitHub (http://github.com/justincbagley) for     #
#   further information. Last update: January 3, 2016. For questions, please email       #
#   jcbagley@unb.br.                                                                     #
##########################################################################################

echo "
##########################################################################################
#                             GaussClust v0.1.0, January 2016                              #
##########################################################################################
"
echo "INFO      | $(date) | STEP #1: SETUP. SETTING OPTIONS AND PATH VARIABLE... "

############ SCRIPT OPTIONS
## OPTION DEFAULTS ##
NUM_NMDS_DIMS=4
CALL_UNSUPERGMM=1
RANGE_NBCLUST=0
CALL_DISCRIMINANT=1
MY_PROBS_MATRIX=probs.txt
CALL_BGMM=0

## PARSE THE OPTIONS ##
while getopts 'k:u:r:d:b:p:c:' opt ; do
  case $opt in
    k) NUM_NMDS_DIMS=$OPTARG ;;
    u) CALL_UNSUPERGMM=$OPTARG ;;
    r) RANGE_NBCLUST=$OPTARG ;;
    d) CALL_DISCRIMINANT=$OPTARG ;;
    b) CALL_BGMM=$OPTARG ;;
    p) MY_PROBS_MATRIX=$OPTARG ;;
    c) NUM_COMPONENTS=$OPTARG ;;
  esac
done

## SCRIPT USAGE ##
##--Check for mandatory positional parameters and echo usage, then wait for commands...
shift $((OPTIND-1))
if [ $# -lt 1 ]; then
  echo "
Usage: $0 [options] inputFile
  "
  echo "Options: -k nmdsDimensions (specify number of dimensions, k, to retain during NMDS \
on Gower distances) | -u unsuperGMM (0=no unsupervised GMM is carried out; 1=conduct unsupervised GMM \
using 'Rmixmod' R package, for comparative or individual purposes) | -r rangeNumClusters (optional \
numeric listing of a range, x:y, of the number of clusters to be modeled over during unsupervised GMM \
in Rmixmod) | -d GMMDiscrim (0=(semi-)supervised, GMM-based discriminant analysis is not carried out \
with 'mixmodelLearn' in Rmixmod; 1=conduct discriminant analysis) | -b beliefBasedMM (0=no mixture \
modeling is carried out using the 'bgmm' R package; 1=calls 'supervised' GMM analysis, 2=calls \
'semisupervised' GMM analysis, and 3=calls both supervised and semisupervised analyses in bgmm) | -p \
probsMatrix (specify matrix of plausibilities, or weights of the prior probabilities for labeled \
observations, for bgmm; also sets belief matrix) | -c numComponents (specify number of components \
(e.g. Gaussian components) or 'clusters' to assign individuals to during regular GMM (single value, \
rather than a range; see -r above) or bgmm modeling)

The -k flag sets the number of k dimensions to be retained during NMDS, which affects both
regular Gaussian mixture modeling and also the different models that are implemented in
the bgmm R package. Like file name, there is no default value; however, k=4 is recommended
by the authors based on discussion in Edwards and Knowles (Proc. Roy. Soc. B. 2014) and 
Hausdorf and Hennig (Syst. Biol. 2014). (By contrast, k=2 would be normal for most other
ecological data, but may not contain sufficient information for interspecific datasets.)

The -u flag calls the unsupervised Gaussian mixture modeling method implemented in the 
'mixmodCluster' function of the Rmixmod R package. See the Rmixmod R site and documentation
for additional information on this package (available at: 
https://cran.r-project.org/web/packages/Rmixmod/index.html). Set this flag to '0' to
skip this analysis.

The -r flag gives the user the ability to conduct unsupervised modeling (called using -u 
above) over a range of nbCluster values. In the case that rangeNumClusters is specified 
(e.g. as '5:20'), Rmixmod will calculate unsupervised GMMs over this range and select the 
best model using the Bayesian information criterion (BIC). If a range of values is not 
specified for -r, then a GMM analysis in Rmixmod will use the number of components/clusters 
specified using the -c flag (see below).

The -d flag calls the supervised or semi-supervised discriminant analysis method implemented 
in the 'mixmodLearn' and 'mixmodPredict' functions of Rmixmod. The discriminant analysis is
based on GMMs and is conducted in a two-step (A, Learning; B, Prediction) procedure, which 
estimates a discriminant function from known labeled data and uses it to predict (classify) 
unknown samples that correspondto the same knowns, i.e. species or clusters. Set this flag 
to '0' to skip this analysis.

The -b flag allows users to request two belief-based Gaussian mixture modeling options 
available in the 'bgmm' R package. The two currently supported models are specified in 
different functions by passing the script a value of '1', which calls the 'supervised' 
function for supervised GMM analysis, or '2', which calls the 'semisupervised' function 
for semisupervised GMM analysis. You can also call both of these functions by passing a 
value of '3' to this option. See the bgmm R site and documentation for more information 
(available at: https://cran.r-project.org/web/packages/bgmm/index.html). Set this flag 
to '0' to skip this analysis.

The -p flag specifies the filename of the bgmm 'B' matrix file in the working dir.

The -c flag specifies the number of components or 'clusters' that will be modeled during
regular GMM or bgmm modeling (except see other option available using -r flag above). This 
corresponds to 'k' or the number of columns in 'B', based on definitions in the bgmm 
documentation.

Input file: Script expects as inputFile a single plain-text data table with a header and 
several columns of information followed by columns containing single-type or mixed data
(e.g. categorical, discrete, or continuous data for different morphological characters 
measured) for the sample. The first column will be named 'samples' and typically contain 
sample IDs/codes for each individual (preferably with a species-specific abbreviation 
followed by a museum voucher number or individual code). The second column is headed as 
'type' and specifies whether each individual ID is 'known' or 'unknown'. The third column 
contains labels (e.g. four-letter codes) for each known individual (e.g. by species), and
'NA' values for samples of unknown type, assigning individuals to species. The example 
input file contains a header with four-letter codes for each datacolumn, but users can 
make the names a little longer if needed.
"

  exit 1
fi


## Make input file a mandatory parameter, and set the path variable to the current dir:
	MY_INPUT_FILE="$1"
	MY_PATH=`pwd -P`


##--FIX issues with echoing shell text containing dollar signs to R:
	MY_POINTS_VAR=$(echo "\$points")	## Make points variable with '$points' text for Rscript...
	MY_TYPE_VAR=$(echo "\$type")		## Make type variable with '$type' text for Rscript...
	MY_SAMP_VAR=$(echo "\$samples")		## Same as above but for '$samples'...
	MY_SPECIES_VAR=$(echo "\$species")	## Same as above but for '$species'...
	MY_STRESS_VAR=$(echo "\$stress")	## Same as above but for '$stress'...
	MY_NMDS1_VAR=$(echo "\$nmds_1")
	MY_NMDS2_VAR=$(echo "\$nmds_2")
	MY_NMDS3_VAR=$(echo "\$nmds_3")
	MY_NMDS4_VAR=$(echo "\$nmds_4")
	MY_TIJ_VAR=$(echo "\$tij")


############ MAKE R SCRIPT
echo "INFO      | $(date) | STEP #2: MAKE GAUSSIAN CLUSTERING R SCRIPT CONTAINING ENVIRONMENTAL VARIABLES AND ANALYSIS CODE... "

echo "
#!/usr/bin/env Rscript

#################################### GaussClust.R ########################################

############ CONDUCT SETUP, READ IN AND PLOT THE DATA
setwd('$MY_PATH')
# 
##--Load needed library, R code, or package stuff. Install packages if not present.
packages <- c('bgmm', 'Rmixmod', 'StatMatch', 'MASS', 'ggfortify', 'vegan')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(bgmm)
library(Rmixmod)
library(StatMatch)
library(MASS)
library(ggfortify)
library(vegan)

##--Read in the data:
mydata_names <- read.table('$MY_INPUT_FILE', h=T)
str(mydata_names)
mydata <- mydata_names[,-c(1:3)]
str(mydata)

##--Graph of pairwise data plots when there are small numbers of characters (else basic R 
##--function and windows cannot handle the plots resulting from calling 'plot' function):
if( dim(mydata_names)[2] < 10 ){
print('Making pairwise plots of the data... ')
pdf('pairwise_data_plots.pdf')
plot(mydata_names)
dev.off()} else {print('Cannot do pairwise plots of the data (too much data)... ')
}


############ ANALYSIS

############ I. SCALE / STANDARDIZE DATA USING NMDS 
##--Estimate Gower distances from the original data. Here, we are starting from a single 
##--data frame, and the code was originally written with a single morphological data matrix
##--in mind. However, multiple datasets from different data types could be used. Either way,
##--Gower distances are ideal because they allow us to get an estimate of the similarities
##--and dissimilarities between all individuals/samples, while allowing for missing data 
##--('NA' observations) in the dataset.
mydata_gower <- gower.dist(mydata)

##--Use ggfortify to visualize the Gower distnces:
pdf('gower_dist_gg_autoplot.pdf')
autoplot(mydata_gower)
dev.off()

##--Conduct NMDS using simple 'nmds' function from labdsv pkg, RETAINING k DIMENSIONS,
##--with the goal being to get NMDS points for use downstream in Rmixmod and/or bgmm:
# mydata_gower_nmds <- nmds(mydata_gower, k=$NUM_NMDS_DIMS)
# pdf('gower_nmds_plot.pdf')
# plot(mydata_gower_nmds)
# dev.off()
##--Use ggfortify to visualize the NMDS:
# pdf('gower_gg_nmds_plot.pdf')
# autoplot(isoMDS(mydata_gower, k=$NUM_NMDS_DIMS), colour = 'orange', size = 1, shape = 3)
# dev.off()

##--Conduct NMDS using 'metaMDS' function in vegan package, retaining k dimensions:
mydata_gower_metaMDS <- metaMDS(mydata_gower, k=$NUM_NMDS_DIMS)
summary(mydata_gower_metaMDS)
#
nmds_stress <- mydata_gower_metaMDS$MY_STRESS_VAR * 100
nmds_stress
#
mydata_gower_metaMDS$MY_POINTS_VAR
metaMDS_points <- as.data.frame(mydata_gower_metaMDS$MY_POINTS_VAR)
metaMDS_points
row.names(metaMDS_points) <- mydata_names[,1]
names(metaMDS_points) <- c('nmds_1', 'nmds_2', 'nmds_3', 'nmds_4')
#
pdf('gower_nmds_plot.pdf')
plot(metaMDS_points)
dev.off()
pdf('gower_nmds_plot_1vs2.pdf')
plot(metaMDS_points$MY_NMDS1_VAR, metaMDS_points$MY_NMDS2_VAR, xlab='NMDS dim 1', ylab='NMDS dim 2')
text(-0.12,0.12, round(c(mydata_gower_metaMDS$MY_STRESS_VAR * 100), digits=2), col='red')
dev.off()
pdf('gower_nmds_plot_1vs3.pdf')
plot(metaMDS_points$MY_NMDS1_VAR, metaMDS_points$MY_NMDS3_VAR, xlab='NMDS dim 1', ylab='NMDS dim 3')
text(-0.12,0.12, round(c(mydata_gower_metaMDS$MY_STRESS_VAR * 100), digits=2), col='red')
dev.off()
pdf('gower_nmds_plot_2vs3.pdf')
plot(metaMDS_points$MY_NMDS2_VAR, metaMDS_points$MY_NMDS3_VAR, xlab='NMDS dim 2', ylab='NMDS dim 3')
text(-0.12,0.12, round(c(mydata_gower_metaMDS$MY_STRESS_VAR * 100), digits=2), col='red')
dev.off()


##--Save each dimension of values retained from NMDS into a separate variable, and then
##--in a data frame (extension 'df'). NOTE: This and other code is weak in being written 
##--to assume that the user will retain four NMDS dimensions (-k 4).
nmds_1 <- metaMDS_points[,1]
nmds_2 <- metaMDS_points[,2]
nmds_3 <- metaMDS_points[,3]
nmds_4 <- metaMDS_points[,4]
nmds_dims_df = data.frame(nmds_1, nmds_2, nmds_3, nmds_4)


############ II. PREP AND CHECK DATA FOR GMM ANALYSES
##--Make data frame containing the individual sample names as well as columns of points
##--from NMDS dimensions retained during STEP I above. Also save the new data frame(s) to
##--file(s) in the working dir; we may want it handy in case we need it later...
sample_names <- mydata_names[,1]
type <- mydata_names[,2]
species <- mydata_names[,3]
mydata_names_df <- data.frame(sample_names, type, species, nmds_1, nmds_2, nmds_3, nmds_4)
write.table(mydata_names_df, file='mydata_names_df.txt')

##--Subset the NMDS points by 'known' and 'unknown' individuals. We also stop to write the
##--resulting new data frames back to file in working dir--in case of subsequent checks: 
attach(mydata_names_df)
known_0 <- mydata_names_df[ which(mydata_names_df$MY_TYPE_VAR=='known'), ]
detach(mydata_names_df)
row.names(known_0) <- known_0[,1]
known_0
write.table(known_0, file='known_0.txt')
str(known_0)
#
knowns <- known_0[,-c(1:3)]
knowns
dim(knowns)[1]
write.table(knowns, file='knowns.txt')
str(knowns)
#
known_labels <- subset(mydata_names_df$MY_SPECIES_VAR, type=='known')
length(known_labels)
known_labels
#
attach(mydata_names_df)
unknown_0 <- mydata_names_df[ which(mydata_names_df$MY_TYPE_VAR=='unknown'), ]
detach(mydata_names_df)
row.names(unknown_0) <- unknown_0[,1]
unknown_0 <- unknown_0[,-c(1:3)]
write.table(unknown_0, file='unknown_0.txt')
str(unknown_0)
#
unknown_labels <- subset(mydata_names_df$MY_SPECIES_VAR, type=='unknown')
length(unknown_labels)
unknown_labels

##--Remove names from data frame of unknowns and place in var 'X' (bgmm unknowns var). Also 
##--read in the belief matrix (B) containing prior probabilities for the knowns, with 0.95 
##--probability for 'known' labeled individuals and all other cells receiving probs of
##--0.05/k (where k is number of components or clusters, and number of columns in B). 
##--***IMPORTANT***: matrix B is generated by the user prior to running this script and
##--ONLY contains individuals classified as 'knowns'.
row.names(unknown_0) <- c(1:dim(unknown_0)[1])
unknowns <- unknown_0
X <- unknowns
X
dim(X)
#
B <- read.table('$MY_PROBS_MATRIX', header=TRUE, sep='\t')
names(B) <- c(0:$NUM_COMPONENTS)
row.names(B) <- B[,1]
B <- B[,-c(1)]
B
dim(B)



############ III. CONDUCT UNSUPERVISED GMM CLUSTERING ANALYSIS IN RMIXMOD
if( $CALL_UNSUPERGMM == '0' ){print('Skipping unsupervised GMM analysis... ')} else {if($CALL_UNSUPERGMM == '1' ){
date()
print('Conducting unsupervised GMM analysis of the data using Rmixmod... ')
#
if( $RANGE_NBCLUST == '0'){print('Running unsupervised GMM using a single cluster number... ')
mydata_gower_gmm <-mixmodCluster(nmds_dims_df, nbCluster=$NUM_COMPONENTS)
summary(mydata_gower_gmm)
pdf('gower_gmm_result.pdf')  ## SAVE THIS PLOT!
plot(mydata_gower_gmm)
dev.off()
mydata_gower_gmm['partition']
} else {print('Running unsupervised GMMs over a range of values for nbCluster, then selecting the best model using BIC... ')
mydata_gower_gmm <-mixmodCluster(nmds_dims_df, nbCluster=$RANGE_NBCLUST)
summary(mydata_gower_gmm)
pdf('gower_gmm_result.pdf')  ## SAVE THESE PLOTS--THEY'RE AWESOME!!!
plot(mydata_gower_gmm)
plot(mydata_gower_gmm, c(1, 2))
plot(mydata_gower_gmm, c(1, 3))
plot(mydata_gower_gmm, c(1, 4))
plot(mydata_gower_gmm, c(2, 3))
plot(mydata_gower_gmm, c(2, 4))
plot(mydata_gower_gmm, c(3, 4))
dev.off()
mydata_gower_gmm['partition']}
	}
}

####### IV. (SEMI-)SUPERVISED GMM-BASED DISCRIMINANT ANALYSIS IN RMIXMOD:
if( $CALL_DISCRIMINANT == '0' ){print('Skipping GMM-based discriminant analysis in Rmixmod... ')} else {if($CALL_DISCRIMINANT == '1' ){
## Analysis:
##    A. Learning:
mydata_known_learn <- mixmodLearn(as.data.frame(knowns), as.factor(known_labels), nbCVBlocks = 10)
mydata_known_learn['bestResult']
graphics.off()
quartz()
pdf('superRmixmod_learn_results.pdf')   ## SAVE THESE PLOTS--THEY'RE AWESOME!!!
plot(mydata_known_learn)
plot(mydata_known_learn, c(1, 2))
plot(mydata_known_learn, c(1, 3))
plot(mydata_known_learn, c(1, 4))
plot(mydata_known_learn, c(2, 3))
plot(mydata_known_learn, c(2, 4))
plot(mydata_known_learn, c(3, 4))
dev.off()
#
##    B. Prediction:
## My prior experience with this (supervised prediction on lizard morphological dataset)
## suggests that prediction success when going from a set of knowns (1:1 or partial coverage)
## to unknowns is usually not so good (~20%). Nevertheless, here goes:
mydata_unknown_prediction <- mixmodPredict(data = X, classificationRule = mydata_known_learn['bestResult'])
summary(mydata_unknown_prediction)
mean(as.integer(unknown_labels) == mydata_unknown_prediction['partition'])
	}
}

####### V. (SEMI-)SUPERVISED BELIEF-BASED GMM ANALYSIS IN BGMM:
if( $CALL_BGMM == '0' ){print('Skipping belief-based GMM analysis in bgmm... ')} else if($CALL_BGMM == '1' ){
supervisedModel <- supervised(as.data.frame(knowns), class = as.factor(known_labels))
supervisedModel
# pdf('bgmm_supervised_result.pdf')
# plot(supervisedModel)
# dev.off()
} else if($CALL_BGMM == '2' ){
semisupervisedModel <- semisupervised(as.data.frame(X), as.data.frame(knowns), class = as.factor(known_labels), k = $NUM_COMPONENTS, P = B)
$semisupervisedModel
pdf('bgmm_semisupervised_result.pdf')
plot(semisupervisedModel)
dev.off()
z <- as.data.frame(semisupervisedModel$MY_TIJ_VAR)
write.table(z, file='bgmm_semisupervised_posteriorProbs.txt', sep='\t')} else if($CALL_BGMM == '3' ){
supervisedModel <- supervised(as.data.frame(knowns), class = as.factor(known_labels))
supervisedModel
# pdf('bgmm_supervised_result.pdf')
# plot(supervisedModel)
# dev.off()
semisupervisedModel <- semisupervised(as.data.frame(X), as.data.frame(knowns), class = as.factor(known_labels), k = $NUM_COMPONENTS, P = B)
$semisupervisedModel
pdf('bgmm_semisupervised_result.pdf')
plot(semisupervisedModel)
dev.off()
z <- as.data.frame(semisupervisedModel$MY_TIJ_VAR)
write.table(z, file='bgmm_semisupervised_posteriorProbs.txt', sep='\t')} else {print('Sorry, the belief, soft, and unsupervised routines in bgmm are not yet supported in GaussClust... ')}





######################################### END ############################################
" > GaussClust.r



############ FINAL STEPS:
echo "INFO      | $(date) | STEP #3: RUN THE R SCRIPT. "
	R CMD BATCH ./GaussClust.R

echo "INFO      | $(date) | STEP #4: CLEAN UP THE WORKSPACE. "
##--Cleanup:
echo "INFO      | $(date) |          Moving R output files to new folder named 'R'... "
	mkdir R
	mv ./*.pdf ./*.Rout ./R/
	if [ "$CALL_BGMM" -gt "0" ]; then
		mv ./bgmm_semisupervised_posteriorProbs.txt ./R/
	fi

## Next: some questions-based flow control for the cleanup...
## Rscript:
	read -p "FLOW      | $(date) |          Would you like to keep the Rscript output by GaussClust? (y/n) : " DEL_SCRIPT
	if [ "$DEL_SCRIPT" != "y" ]; then
		rm ./GaussClust.r
	else
		mv ./GaussClust.r ./R/
	fi

## Text files:
	read -p "FLOW      | $(date) |          Would you like to keep text files output by GaussClust? (y/n) : " TO_KEEP
	if [ "$TO_KEEP" = "y" ]; then
		echo "INFO      | $(date) |          Moving text files to new folder named 'txt'... "
		mkdir txt
		mv ./known_0.txt ./knowns.txt ./unknown_0.txt ./mydata_names_df.txt ./txt/
	else
		rm ./known_0.txt ./knowns.txt ./unknown_0.txt ./mydata_names_df.txt
	fi


echo "INFO      | $(date) | Done conducting Gaussian clustering and related analyses using GaussClust."
echo "INFO      | $(date) | Bye.
"
#
#
#
######################################### END ############################################

exit 0