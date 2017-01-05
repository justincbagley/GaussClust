
#!/usr/bin/env Rscript

################################### compareProbs.R #######################################

############ I. SETUP
setwd('/Users/justinbagley/Documents/GaussClust/Ex2_bgmmSensTest')
#
##--Load needed library, R code, or package stuff. Install packages if not present.
packages <- c('tools', 'ggfortify', 'grid', 'gplots')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(tools)
library(ggfortify)
library(grid)
library(gplots)


############ II. READ AND PREP tij MATRICES FOR PLOTTING
##--Create vector list of filenames for the posterior probability matrices from each run:
tijFileNames <- list.files(getwd(), pattern='^tij')
#
##--The following function will rotate the matrices into correct orientation during image
##--plotting; the idea for this rotation method is from the following site: https://blog.snap.uaf.edu/2012/06/08/matrix-rotation-for-image-and-contour-plots-in-r/
f <- function(m) t(m)[,nrow(m):1]
#
##--Loop through the filenames to read each file in and assign it to the appropriately 
##--named object:
for(i in 1:length(tijFileNames)){
#	x <- tij_file.names[i]
	basename <- file_path_sans_ext(tijFileNames[i])
	oname = paste(basename)	
	assign(basename, as.matrix(read.table(tijFileNames[i], header=TRUE, sep='	')))
}

############ III. PLOT THE MATRICES IN A SINGLE GRAPH TO PERMIT BROAD VISUAL COMPARISON OF RESULTS ACROSS RANGE OF TEST VALUES
##--Plot the matrices in the quartz window with par-specified parameters (given above):
tijObjList <- list(tij_0.5, tij_0.59, tij_0.68, tij_0.77, tij_0.86, tij_0.95)
pdf('plots.pdf')

##--Use par and col to specify the plot layout and color characteristics:
par(mfrow=c(3,2), mai = c(0.5, 0.1, 0.1, 0.1))
col<- colorRampPalette(c('ivory2', 'yellow', 'red'))(20)

for(i in tijObjList){
	image(f(i), col=col)
}
dev.off()



##--Some logical testing to see if any matrices are the same (likely, they are NOT):
identical(tijObjList[1][[1]][1:31,],tijObjList[2][[1]][1:31,]);
identical(tijObjList[1][[1]][1:31,],tijObjList[3][[1]][1:31,]);
identical(tijObjList[1][[1]][1:31,],tijObjList[4][[1]][1:31,]);
identical(tijObjList[1][[1]][1:31,],tijObjList[5][[1]][1:31,]);
identical(tijObjList[1][[1]][1:31,],tijObjList[6][[1]][1:31,]);
identical(tijObjList[2][[1]][1:31,],tijObjList[3][[1]][1:31,]);
identical(tijObjList[2][[1]][1:31,],tijObjList[4][[1]][1:31,]);
identical(tijObjList[2][[1]][1:31,],tijObjList[5][[1]][1:31,]);
identical(tijObjList[2][[1]][1:31,],tijObjList[6][[1]][1:31,]);
identical(tijObjList[3][[1]][1:31,],tijObjList[4][[1]][1:31,]);
identical(tijObjList[3][[1]][1:31,],tijObjList[5][[1]][1:31,]);
identical(tijObjList[3][[1]][1:31,],tijObjList[6][[1]][1:31,]);
identical(tijObjList[4][[1]][1:31,],tijObjList[5][[1]][1:31,]);
identical(tijObjList[4][[1]][1:31,],tijObjList[6][[1]][1:31,]);
identical(tijObjList[5][[1]][1:31,],tijObjList[6][[1]][1:31,]);


##--Some mathematical/arithmetic testing:
##--Calculate the differences between the posterior probability matrices (tij's, placed in
##--'tijObjList' above), which gives one difference matrix for each pair of tij matrices. 
##--Then sum across all difference matrices to get the cumulative difference across them.
##--Compare this to the sum of all the original tij matrices (which since each probs matrix 
##--row sums to one is simply equal to the number of rows minus the header, in each probs 
##--matrix). Two tests might be whether 1) the differences between any pair of matrices is 
##--more/less than five percent of the sum of entries in the 2 matrices, or 2) whether the
##--sum of differences across all difference matrices is more/less than five percent of the
##--sum of entries across all difference matrices.
#
##--First, put all individual tij matrices, calculated from probs matrices with different 
##--on/off-diagonal test values, into separate objects, using a for loop:
for(i in 1:6) {
    name <- paste('m', i, sep = '')
    assign(name, tijObjList[i][[1]])
}
##--Next, calculate a difference matrix for each pair of tij matrices:
dif_m12 <- m1 - m2; 
dif_m13 <- m1 - m3; 
dif_m14 <- m1 - m4; 
dif_m15 <- m1 - m5; 
dif_m16 <- m1 - m6; 
dif_m23 <- m2 - m3; 
dif_m24 <- m2 - m4; 
dif_m25 <- m2 - m5; 
dif_m26 <- m2 - m6; 
dif_m34 <- m3 - m4; 
dif_m35 <- m3 - m5; 
dif_m36 <- m3 - m6; 
dif_m45 <- m4 - m5; 
dif_m46 <- m4 - m6; 
dif_m56 <- m5 - m6; 
#
##--Next, get the sum of all differences across the pairwise diff matrices, as well as the
##--sum of all tij matrices themselves. Then compare these two to do the 'test' #2 mentioned
##--above under 'Some mathematical/arithmetic testing', by simple division:
sum_difm <- sum(c(dif_m12, dif_m13, dif_m14, dif_m15, dif_m16, dif_m23, dif_m24, dif_m25, dif_m26, dif_m34, dif_m35, dif_m36, dif_m45, dif_m46, dif_m56))
sum_difm
sum_tijm <- sum(c(m1, m2, m3, m4, m5, m6))
sum_tijm
sum_difm/sum_tijm	## Is this value > < or = to 0.05?
#
##--Now, source the ./diffProTests.r code generated in the final section of the sh code
##--before this Rscript, in order to calculate the difference proportion of total (from 
##--original two tij matrices) for each pairwise difference matrix. 
source('diffProTests.r')
#
diffpro_list <- as.list(apropos('diffpro'))
#
for(i in 3:length(diffpro_list)){
	cat('---', diffpro_list[i][[1]][1],'---
');
	print(get(diffpro_list[[i]][1]))
}
diffpros <- abs(c(diffpro_m12, diffpro_m13, diffpro_m14, diffpro_m15, diffpro_m16, diffpro_m23, diffpro_m24, diffpro_m25, diffpro_m26, diffpro_m34, diffpro_m35, diffpro_m36, diffpro_m45, diffpro_m46, diffpro_m56))

##--Print min, max, average, and range of difference proportions across diffpros values.
##--Ideally, all of these would be minimal, e.g. <0.05 (arbitrary value).
min(diffpros)
max(diffpros)
mean(diffpros)
range(diffpros)

##--Now, plot out the difference proportions of total across all pairs, and save the 
##--resulting line graph as a PDF file. We can do this twice, once with the diff pros
##--sorted in ascending order, and the other with them in their original order:
pdf('diffpro_sort.pdf')
plot(sort(diffpros, decreasing=FALSE), type='o', col='blue')
dev.off()
#
pdf('diffpro_unsort.pdf')
plot(diffpros, type='o', col='blue')
dev.off()
#
pdf('diffpro_unsort_ylim0.1.pdf')
plot(diffpros, type='o', col='blue', xlim=c(0, 15), ylim=c(0, 0.1))
dev.off()
#
pdf('diffpro_unsort_ylim1.pdf')
plot(diffpros, type='o', col='blue', xlim=c(0, 15), ylim=c(0, 1))
dev.off()
#


######################################### END ############################################

