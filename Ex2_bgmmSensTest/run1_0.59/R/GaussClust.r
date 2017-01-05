
#!/usr/bin/env Rscript

#################################### GaussClust.R ########################################

############ CONDUCT SETUP, READ IN AND PLOT THE DATA
setwd('/Users/justinbagley/Documents/GaussClust/Ex2_bgmmSensTest/run1_0.59')
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
mydata_names <- read.table('Enyalius_35.txt', h=T)
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

##--Conduct NMDS using 'metaMDS' function in vegan package, retaining k dimensions:
mydata_gower_metaMDS <- metaMDS(mydata_gower, k=4)
summary(mydata_gower_metaMDS)
#
nmds_stress <- mydata_gower_metaMDS$stress * 100
nmds_stress
#
mydata_gower_metaMDS$points
metaMDS_points <- as.data.frame(mydata_gower_metaMDS$points)
metaMDS_points
row.names(metaMDS_points) <- mydata_names[,1]
names(metaMDS_points) <- c('nmds_1', 'nmds_2', 'nmds_3', 'nmds_4')
#
pdf('gower_nmds_plot.pdf')
plot(metaMDS_points)
dev.off()
pdf('gower_nmds_plot_1vs2.pdf')
plot(metaMDS_points$nmds_1, metaMDS_points$nmds_2, xlab='NMDS dim 1', ylab='NMDS dim 2')
text(-0.12,0.12, round(c(mydata_gower_metaMDS$stress * 100), digits=2), col='red')
dev.off()
pdf('gower_nmds_plot_1vs3.pdf')
plot(metaMDS_points$nmds_1, metaMDS_points$nmds_3, xlab='NMDS dim 1', ylab='NMDS dim 3')
text(-0.12,0.12, round(c(mydata_gower_metaMDS$stress * 100), digits=2), col='red')
dev.off()
pdf('gower_nmds_plot_2vs3.pdf')
plot(metaMDS_points$nmds_2, metaMDS_points$nmds_3, xlab='NMDS dim 2', ylab='NMDS dim 3')
text(-0.12,0.12, round(c(mydata_gower_metaMDS$stress * 100), digits=2), col='red')
dev.off()


##--Save each dimension of values retained from NMDS into a separate variable, and then
##--in a data frame (extension 'df'):
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
known_0 <- mydata_names_df[ which(mydata_names_df$type=='known'), ]
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
known_labels <- subset(mydata_names_df$species, type=='known')
length(known_labels)
known_labels
#
attach(mydata_names_df)
unknown_0 <- mydata_names_df[ which(mydata_names_df$type=='unknown'), ]
detach(mydata_names_df)
row.names(unknown_0) <- unknown_0[,1]
unknown_0 <- unknown_0[,-c(1:3)]
write.table(unknown_0, file='unknown_0.txt')
str(unknown_0)
#
unknown_labels <- subset(mydata_names_df$species, type=='unknown')
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
B <- read.table('./run1_0.59.txt', header=TRUE, sep='	')
names(B) <- c(0:2)
row.names(B) <- B[,1]
B <- B[,-c(1)]
B
dim(B)



############ III. CONDUCT UNSUPERVISED GMM CLUSTERING ANALYSIS IN RMIXMOD
if( 1 == '0' ){print('Skipping unsupervised GMM analysis... ')} else {if(1 == '1' ){
date()
print('Conducting unsupervised GMM analysis of the data using Rmixmod... ')
#
if( 0 == '0'){print('Running unsupervised GMM using a single cluster number... ')
mydata_gower_gmm <-mixmodCluster(nmds_dims_df, nbCluster=2)
summary(mydata_gower_gmm)
pdf('gower_gmm_result.pdf')  ## SAVE THIS PLOT!
plot(mydata_gower_gmm)
dev.off()
mydata_gower_gmm['partition']
} else {print('Running unsupervised GMMs over a range of values for nbCluster, then selecting the best model using BIC... ')
mydata_gower_gmm <-mixmodCluster(nmds_dims_df, nbCluster=0)
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
if( 0 == '0' ){print('Skipping GMM-based discriminant analysis in Rmixmod... ')} else {if(0 == '1' ){
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
if( 3 == '0' ){print('Skipping belief-based GMM analysis in bgmm... ')} else if(3 == '1' ){
supervisedModel <- supervised(as.data.frame(knowns), class = as.factor(known_labels))
supervisedModel
# pdf('bgmm_supervised_result.pdf')
# plot(supervisedModel)
# dev.off()
} else if(3 == '2' ){
semisupervisedModel <- semisupervised(as.data.frame(X), as.data.frame(knowns), class = as.factor(known_labels), k = 2, P = B)

pdf('bgmm_semisupervised_result.pdf')
plot(semisupervisedModel)
dev.off()
z <- as.data.frame(semisupervisedModel$tij)
write.table(z, file='bgmm_semisupervised_posteriorProbs.txt', sep='	')} else if(3 == '3' ){
supervisedModel <- supervised(as.data.frame(knowns), class = as.factor(known_labels))
supervisedModel
# pdf('bgmm_supervised_result.pdf')
# plot(supervisedModel)
# dev.off()
semisupervisedModel <- semisupervised(as.data.frame(X), as.data.frame(knowns), class = as.factor(known_labels), k = 2, P = B)

pdf('bgmm_semisupervised_result.pdf')
plot(semisupervisedModel)
dev.off()
z <- as.data.frame(semisupervisedModel$tij)
write.table(z, file='bgmm_semisupervised_posteriorProbs.txt', sep='	')} else {print('Sorry, the belief, soft, and unsupervised routines in bgmm are not yet supported in GaussClust... ')}




######################################### END ############################################

