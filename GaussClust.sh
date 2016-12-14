#!/bin/sh

##########################################################################################
#                             GaussClust v1.0, December 2016                             #
#   SHELL SCRIPT FOR AUTOMATING GAUSSIAN MIXTURE MODELING ON CONTINUOUS, DISCRETE, OR    #
#   COMBINED GENETIC/MORPHOLOGICAL/ECOLOGICAL DATA, FOR SPECIES DELIMITATION AND/OR      #
#   CLASSIFICATION                                                                       #
#   Copyright (c)2016 Justin C. Bagley, Universidade de Brasília, Brasília, DF, Brazil.  #
#   See the README and license files on GitHub (http://github.com/justincbagley) for     #
#   further information. Last update: December 14, 2016. For questions, please email     #
#   jcbagley@unb.br.                                                                     #
##########################################################################################

echo "
##########################################################################################
#                             GaussClust v1.0, December 2016                             #
##########################################################################################
"
echo "INFO      | $(date) | STEP #1: SETUP. SETTING OPTIONS AND PATH VARIABLE... "

############ SCRIPT OPTIONS
## OPTION DEFAULTS ##
NUM_NMDS_DIMS=4
MY_PROBS_MATRIX=probs.txt
CALL_REG_GMM=1
CALL_SUPERMIXMOD=1
CALL_BGMM=0

## PARSE THE OPTIONS ##
while getopts 'd:r:n:s:b:p:c:' opt ; do
  case $opt in
    d) NUM_NMDS_DIMS=$OPTARG ;;
    r) CALL_REG_GMM=$OPTARG ;;
    n) NUM_RANGE_GMM_NBCLUSTERS=$OPTARG ;;
    s) CALL_SUPERMIXMOD=$OPTARG ;;
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
  echo "Options: -d numNMDSDimensions (specify number of dimensions, k, to retain during NMDS \
on Gower distances) | -r regGMM (0=no unsupervised GMM is carried out; 1=conduct unsupervised GMM \
using 'Rmixmod' R pacakge, for comparative or individual purposes) | -n numGMMClusters (optional \
numeric listing of a range, x:y, of the number of clusters to be modeled over during unsupervised GMM \
in Rmixmod) | -s superGMM (0=no (semi-)supervised GMM is carried out in Rmixmod; 1=conduct (semi-)\
supervised GMM in Rmixmod) | -b beliefBasedMM (0=no mixture modeling is carried out using the 'bgmm' \
R package; the following other options conduct different kinds of mixture modeling using separate \
functions available in bgmm: belief, soft, semisupervised, supervised) | -c numComponents (specify \
number of components (e.g. Gaussian components) or 'clusters' to assign individuals to during regular \
GMM (single value, rather than a range; see -n above) or bgmm modeling) 

The -d flag sets the number of k dimensions to be retained during NMDS, which affects both
regular Gaussian mixture modeling and also the different models that are implemented in
the bgmm R package. Like file name, there is no default value; however, k=4 is recommended
by the authors based on discussion in Edwards and Knowles (Proc. Roy. Soc. B. 2014) and 
Hausdorf and Hennig (Syst. Biol. 2014). (By contrast, k=2 would be normal for most other
ecological data, but may not contain sufficient information for interspecific datasets.)

The -r flag calls the unsupervised Gaussian mixture modeling method implemented in the 
'mixmodCluster' function of the Rmixmod R package. See the Rmixmod R site and documentation
for additional information on this package (available at: 
https://cran.r-project.org/web/packages/Rmixmod/index.html). Set this flag to '0' to
skip this analysis.

The -n flag is *optional* and gives the user the ability to conduct unsupervised modeling
(called using -r above) over a range of nbCluster values. In the case that a numGMMClusters 
range is specified (e.g. '5:20'), Rmixmod will calculate unsupervised GMMs over this range 
and select the best model using the Bayesian information criterion (BIC). If a range of 
values is not specified for -n, then a GMM analysis in Rmixmod will use the number of 
components/clusters specified using the -c flag (see below).

The -s flag calls the (semi-)supervised Gaussian mixture modeling method implemented in
the 'mixmodLearn' and 'mixmodPredict' functions of Rmixmod. Set this flag to '0' to
skip this analysis.

The -b flag allows users to request the Gaussian mixture modeling or belief-based mixture
modeling options available in the 'bgmm' R package. You may call four different models,
specified in different functions in bgmm, by passing the script the function names 'belief', 
'soft', 'semisupervised', or 'supervised'. See the bgmm R site and documentation for 
more information on this method (available at: 
https://cran.r-project.org/web/packages/bgmm/index.html).

The -p flag specifies the filename of the bgmm 'B' matrix file in the working dir.

The -c flag specifies the number of components or 'clusters' that will be modeled during
regular GMM or bgmm modeling (except see other option available using -n flag above). This 
corresponds to 'k' or the number of columns in 'B', based on definitions in the bgmm 
documentation.

Input file: Script expects as inputFile a single plain-text data table with a header and 
several columns of information followed by columns containing categorical or discrete data
(e.g. for different morphological characters measured) for the sample. The first column 
will be named 'samples' and contain code names or IDs for each individual (preferably with 
a species-specific abbreviation followed by a museum voucher number or individual code). 
The second column is headed as 'type' and specifies whether each individual ID is 'known'
or 'unknown'. The third column contains integer values corresponding to codes/numbers (1 
to k, where k is total number of components/clusters) assigning each individual to a 
species (usually, 1 species = 1 cluster). The example input file contains a header with 
four-letter codes for each column, but users can make the names a little longer if needed.
"

  exit 1
fi


## Make input file a mandatory parameter:
MY_INPUT_FILE="$1"

MY_PATH=`pwd -P`
# MY_PATH=$(pwd)

##--FIX some issues with echoing shell text to R below:
## Make points variable with '$points' text for Rscript...
MY_POINTS_VAR=$(echo "\$points")
## Make type variable with '$type' text for Rscript...
MY_TYPE_VAR=$(echo "\$type")
## Same as above but for '$samples'...
MY_SAMP_VAR=$(echo "\$samples")
## Same as above but for '$species'...
MY_SPECIES_VAR=$(echo "\$species")


############ MAKE R SCRIPT
echo "INFO      | $(date) | STEP #2: MAKE GAUSSIAN CLUSTERING R SCRIPT CONTAINING ENVIRONMENTAL VARIABLES AND ANALYSIS CODE... "

echo "
#!/usr/bin/env Rscript

#################################### GaussClust.R ########################################

############ SETUP
setwd('$MY_PATH')
# 
##--Load needed library, R code, or package stuff. Install packages if not present.
##--opt: source('GaussClust.R', chdir = TRUE)
packages <- c('bgmm', 'stringi', 'Rmixmod', 'StatMatch', 'MASS', 'cluster', 'ggfortify',
'labdsv')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(bgmm)
library(stringi)
library(Rmixmod)
library(StatMatch)
library(MASS)
library(cluster)
library(ggfortify)
library(labdsv)

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
##--Estimate Gower distances from the original morphological data matrix w/o names: 
mydata_gower <- gower.dist(mydata)

##--Use ggfortify to visualize the gower distnces:
pdf('gower_dist_gg_autoplot.pdf')
autoplot(mydata_gower)
dev.off()

##--Conduct NMDS using simple 'nmds' function from labdsv pkg, RETAINING k DIMENSIONS,
##--with the goal being to get NMDS points for use downstream in Rmixmod and/or bgmm:
mydata_gower_nmds <- nmds(mydata_gower, k=$NUM_NMDS_DIMS)
pdf('gower_nmds_plot.pdf')
plot(mydata_gower_nmds)
dev.off()

##--Use ggfortify to visualize the NMDS:
pdf('gower_gg_nmds_plot.pdf')
autoplot(isoMDS(mydata_gower, k=$NUM_NMDS_DIMS), colour = 'orange', size = 1, shape = 3)
dev.off()

##--Save each dimension of values retained from NMDS into a separate variable, and then
##--in a data frame (extension 'df'):
nmds_1 <- mydata_gower_nmds$MY_POINTS_VAR[,1]
nmds_2 <- mydata_gower_nmds$MY_POINTS_VAR[,2]
nmds_3 <- mydata_gower_nmds$MY_POINTS_VAR[,3]
nmds_4 <- mydata_gower_nmds$MY_POINTS_VAR[,4]
nmds_dims_df = data.frame(nmds_1, nmds_2, nmds_3, nmds_4)


############ II. PREP AND CHECK DATA FOR GMM ANALYSES
##--Make data frame containing the individual sample names as well as columns of points
##--from NMDS dimensions retained during STEP I above. Also save the new data frame(s) to
##--file(s) in the working dir; we may want it handy in case we need it later...
sample_names <- mydata_names[,1]
type <- mydata_names[,2]
mydata_names_4bgmm_df <- data.frame(sample_names, type, nmds_1, nmds_2, nmds_3, nmds_4)
write.table(mydata_names_4bgmm_df, file='mydata_names_4bgmm_df.txt')

##--Subset the NMDS points by 'known' and 'unknown' individuals, for bgmm. Also write the
##--resulting new data frames back to file in working dir--in case of subsequent checks: 
attach(mydata_names_4bgmm_df)
known_0 <- mydata_names_4bgmm_df[ which(mydata_names_4bgmm_df$MY_TYPE_VAR=='known'), ]
detach(mydata_names_4bgmm_df)
row.names(known_0) <- known_0[,1]
known_0 <- known_0[,-c(1:2)]
known_0
dim(known_0)[1]
write.table(known_0, file='known_0.txt')
str(known_0)
knowns <- known_0
#
attach(mydata_names_4bgmm_df)
unknown_0 <- mydata_names_4bgmm_df[ which(mydata_names_4bgmm_df$MY_TYPE_VAR=='unknown'), ]
detach(mydata_names_4bgmm_df)
row.names(unknown_0) <- unknown_0[,1]
unknown_0 <- unknown_0[,-c(1:2)]
write.table(unknown_0, file='unknown_0.txt')
str(unknown_0)

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
#
##--Here, let's make a couple of arbitrary belief matrices, to analyze and compare:
p1=0.95
B2 <- (matrix(p1, nrow = 206, ncol = 15))
row.names(B2) <- row.names(knowns)
p2=0.066666666666667
B3 <- (matrix(p2, nrow = 206, ncol = 15))
row.names(B3) <- row.names(knowns)



############ III. CONDUCT UNSUPERVISED GMM ANALYSIS IN RMIXMOD
if( $CALL_REG_GMM == '0' ){print('Skipping unsupervised GMM analysis... ')} else {if($CALL_REG_GMM == '1' ){
date()
print('Conducting unsupervised GMM analysis of the data using Rmixmod... ')
mydata_gower_gmm <-mixmodCluster(nmds_dims_df, nbCluster=$NUM_RANGE_GMM_NBCLUSTERS)
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
mydata_gower_gmm['partition']} else { date();
print('WARNING: Something went wrong deciding whether or not to call unsupervised GMM analysis... ')
	}
}


####### IV. (SEMI-)SUPERVISED GMM ANALYSIS IN RMIXMOD:
if( $CALL_SUPERMIXMOD == '0' ){print('Skipping semi-supervised or supervised GMM analysis in Rmixmod... ')} else {if($CALL_SUPERMIXMOD == '1' ){
## (1):
mydata_names_clusters <- read.table('mydata_names_clusters.txt', h=T)
str(mydata_names_clusters)
## (2):
known_1 <- mydata_names_clusters[ which(mydata_names_4bgmm_df$MY_TYPE_VAR=='known'), ]
## (3):
##    (3A):
mydata_known_learn <- mixmodLearn(knowns, known_1$MY_SPECIES_VAR, nbCVBlocks = nrow(mydata_names_4bgmm_df))
mydata_known_learn['bestResult']
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
##    (3B):
## My prior experience with this (supervised prediction on lizard morphological data)
## suggests that prediction success when going from a set of knowns (1:1 or partial coverage)
## to unknowns is usually not so good. Nevertheless, here goes:
unknown_1 <- mydata_names_clusters[ which(mydata_names_4bgmm_df$MY_TYPE_VAR=='unknown'), ]
mydata_unknown_prediction <- mixmodPredict(data = X, classificationRule = mydata_known_learn['bestResult'])
summary(mydata_unknown_prediction)
mean(as.integer(unknown_1$MY_SPECIES_VAR) == mydata_unknown_prediction['partition'])
	}
}

####### V. BGMM MODELING STEPS (IN PROGRESS):
modelSoft1 <- soft(X=X, knowns=knowns, P=B2)
pdf('bgmm_modelSoft1_result.pdf')  ## SAVE THIS PLOT!
plot(modelSoft1)
dev.off()

modelSoft2 <- soft(X=X, knowns=knowns, P=B3)
pdf('bgmm_modelSoft2_result.pdf')  ## SAVE THIS PLOT!
plot(modelSoft2)
dev.off()

quartz()
pdf('bgmm_modelSoft1_vs_modelSoft2_2x1.pdf')  ## SAVE THIS PLOT!
par(mfrow=c(2,1))
plot(modelSoft1)
plot(modelSoft2)
dev.off()


### ADD SEMISUPERVISED ANALYSIS IN BGMM HERE ONCE YOU GET IT WORKING.


######################################### END ############################################
" > GaussClust.r



############ FINAL STEPS:
echo "INFO      | $(date) | STEP #3: RUN THE R SCRIPT. "
R CMD BATCH ./GaussClust.R

echo "INFO      | $(date) | STEP #4: CLEAN UP THE WORKSPACE. "
##--Cleanup:
mkdir R
mv ./bgmm_modelSoft1_vs_modelSoft2_2x1.pdf ./bgmm_modelSoft2_result.pdf \
./bgmm_modelSoft1_result.pdf ./superRmixmod_learn_results.pdf ./gower_gmm_result.pdf \
./gower_gg_nmds_plot.pdf ./gower_nmds_plot.pdf ./gower_dist_gg_autoplot.pdf ./*.Rout ./R/

#
## Next-some questions-based flow control for the cleanup:
read -p "FLOW      | $(date) |          Would you like to keep the Rscript output by GaussClust? (y/n) : " DEL_SCRIPT
if [ "$DEL_SCRIPT" != "y" ]; then
	rm ./GaussClust.r
else
	mv ./GaussClust.r ./R/
fi

read -p "FLOW      | $(date) |          Would you like to keep text files output by GaussClust? (y/n) : " TO_KEEP
if [ "$TO_KEEP" = "y" ]; then
	mkdir txt
	mv ./known_0.txt ./unknown_0.txt ./mydata_names_4bgmm_df.txt ./txt/
else
	rm ./known_0.txt ./unknown_0.txt ./mydata_names_4bgmm_df.txt
fi


echo "INFO      | $(date) | Done conducting Gaussian clustering and related analyses using GaussClust."
echo "INFO      | $(date) | Bye.
"
#
#
#
######################################### END ############################################

exit 0
