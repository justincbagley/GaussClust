#!/bin/sh

##########################################################################################
# File: bgmmSensTest.sh                                                                  #
  VERSION="v0.1.0"                                                                       #
# Author: Justin C. Bagley                                                               #
# Date: created by Justin Bagley on Wed Jan 4 09:09:34 2017 -0600                        #
# Last update: March 3, 2019                                                             #
# Copyright (c) 2017-2019 Justin C. Bagley. All rights reserved.                         #
# Please report bugs to <bagleyj@umsl.edu>.                                              #
#                                                                                        #
# Description:                                                                           #
#  SHELL SCRIPT FOR CONDUCTING SIMPLE SENSITIVITY TESTING TO EXPLORE BELIEF-BASED        #
#  GAUSSIAN MIXTURE MODEL ROBUSTNESS TO VARYING PRIORS ON KNOWN OBSERVATIONS (P MATRIX)  #
#                                                                                        #
##########################################################################################

if [[ "$1" == "-v" ]] || [[ "$1" == "--version" ]]; then
	echo "$(basename $0) $VERSION";
	exit
fi

echo "
##########################################################################################
#                           bgmmSensTest v0.1.0, December 2016                           #
##########################################################################################
"

######################################## START ###########################################
echo "INFO      | $(date) | Starting bgmmSensTest analysis... "
echo "INFO      | $(date) | STEP #1: SETUP AND USER INPUT. "
###### Set paths and filetypes as different variables:
##--Set new path/dir environmental variable to the current directory:
#	MY_PATH=`pwd -P`
	MY_PATH="$(pwd -P)"
echo "INFO      | $(date) |          Setting working directory to: $MY_PATH "
	CR=$(printf '\r')
	calc () {
	   	bc -l <<< "$@"
	}

echo "INFO      | $(date) |          Reading in sensitivity test conditions and file names from configuration file. "
###### Use grep and awk to read in information from the user about the desired sensitivity
##--test conditions, by pulling the required information from the 'bgmm_sens_test.cfg'
##--configuration file supplied by the user, present in current working directory.	
	OPT_STRING_PART="$(grep -n "opt_string_part" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	NUM_COMPONENTS="$(grep -n "num_clusters" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	MY_DATA_FILE="$(grep -n "data_file" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	MY_PROBS_FILE="$(grep -n "probs_file" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	TARGET_VALUE="$(grep -n "target_diag_value" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	NUM_TEST_VALUES="$(grep -n "num_test_vals" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	TEST_MIN="$(grep -n "test_val_min" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	TEST_MAX="$(grep -n "test_val_max" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";
	TEST_INC="$(grep -n "test_val_inc" ./bgmm_sens_test.cfg | awk -F"=" '{print $NF}')";


echo "INFO      | $(date) | STEP #2: MAKE DIFFERENT INPUT FILES FOR THE TEST, VARYING TEST STAT ACROSS A RANGE OF VALUES. "
###### Make functions that make new input files, each of which has the target value changed 
##--to one test value in the range made up by the minimum, increment, and maximum values
##--given by the user above. The functions thus also have to change the corresponding off-
##--diagonal values, to the appropriate 1-y/k value (where y is the on-diagonal value). We 
##--use a shell one-liner (using grep, sed, and head) to pull the off-diagonal values out 
##--of the input file. Because we are dealing with probs matrices composed of probability 
##--values, which sum to one thus can never individually be greater than 1, we are working 
##--with float values, so our function must deal with float values throughout. We name 
##--this function "makeInputFiles".
makeInputFiles () {
count=0

		INPUT_NCOL=$(cat $MY_PROBS_FILE | awk '{print NF}' | head -n1);
		NCOL_DROPONE=$(calc $INPUT_NCOL - 1);

		AFTER_DECI="$(echo $TARGET_VALUE | sed 's/0\.//g')";

		OFF_DIAG_VALUES=$(grep -n ''$TARGET_VALUE'\t' $MY_PROBS_FILE | \
		sed 's/^[0-9\:A-Za-z\_\-]*[	]*//g; s/0\.${AFTER_DECI} //g' | head | \
		sed 's/[0-9\.]*   //g' | head -n1);
		
		NEW_OFF_DIAG_VALUES=$(calc $(calc 1 - $TARGET_VALUE) / $NCOL_DROPONE);
		
	(
		for i in $(seq $TEST_MIN $TEST_INC $TEST_MAX); do
			IFIX="$(echo $i | sed 's/0\.//g')";
			sed 's/0\.${AFTER_DECI}	/0\.'$IFIX'	/g' $MY_PROBS_FILE | sed 's/$OFF_DIAG_VALUES/$NEW_OFF_DIAG_VALUES/g' > ./run"$count"_$i.txt
			MY_RUN_TXT_FILE=./run*.txt;		## There will only be one of these, which we will use and then move into a folder so that it does not conflict with the next txt file generated by this loop.

					basename=$(echo $MY_RUN_TXT_FILE | sed 's/\.\///g; s/\.txt//g');
			
			mkdir $basename;
			mv $MY_RUN_TXT_FILE ./$basename/;

		count=$((count+1))
       
		done
	)
}

##--Don't forget to run the function!
makeInputFiles


echo "INFO      | $(date) | STEP #3: RUN ALGORITHM ACROSS TEST INPUT FILES, SAVE & COLLATE RESULTS. "
###### Move a copy of the data file into each test run folder created using the loop in 
##--STEP #2 above, and then run the GaussClust 'lite' script (in wkdir) from within each 
##--folder. We will specify options so that we run belief-based GMM using the 'bgmm' R 
##--package, using the same exact parameters, on each input file, which are brought in
##--from the OPT_STRING_PART environmental variable (captured from cfg file). We'll also 
##--make sure the output has appropriate names, and we will extract and format, as needed, 
##--any desired information from results of each run.
	(
		for j in ./*/; do
			cp ../GaussClust_lite.sh $MY_DATA_FILE $j; 
			cd $j;
			MY_LOCAL_PROBS_FILE=./run*.txt;
			./GaussClust_lite.sh $OPT_STRING_PART -p $MY_LOCAL_PROBS_FILE -c $NUM_COMPONENTS $MY_DATA_FILE;
			cd ..;
		done
	)

##--Organize posterior probability matrices resulting from each run on a separate input
##--file created in STEP #2 above, by placing one copy named as 'tij' (posterior probs
##--matrix name; see bgmm R package documentation) followed by the test value (e.g. "tij_0.5")
##--into the working directory.
	(
		for k in ./*/; do
			FOLDERNAME=$(echo $k | sed 's/\.\/run[0-9]*\_//g; s/\/$//g');
			cp "$k"R/*_posteriorProbs.txt ./tij_$FOLDERNAME.txt;
		done
	)


###### Prep some values as environmental variables for use in Rscript below:
##--Make environmental variable containing posterior probability (PP) matrix names:
	MY_PP_MATRICES=./tij_*.txt;
	FUTURE_TIJ_OBJ_NAMES=$(echo $MY_PP_MATRICES | sed 's/\.\///g; s/.txt//g; s/\ /\,\ /g'); 	## This is used in the Rscript below.
#
##--Account for possible odd or even value for $NUM_TEST_VALUES
	EVEN_CHECK=$(( $NUM_TEST_VALUES % 2 ))
	if [ $EVEN_CHECK -eq 0 ]; then
		HALF_NUM_TEST_VALUES=$(calc $NUM_TEST_VALUES/2);
	else
		HALF_NUM_TEST_VALUES=$(calc $(calc $NUM_TEST_VALUES +1)/2);
	fi


###### Prep some text and code for math/arithmetic work in Rscript below:
##--Make files with all pairs of combinations of the sequence from 1 to the total number of 
##--tij matrices, which is already $NUM_TEST_VALUES.
	set -- $(seq $NUM_TEST_VALUES)
	for a; do     shift;     for b; do         printf "%s %s\n" "$a" "$b";     done; done > pairs.txt
	set -- $(seq $NUM_TEST_VALUES)
	for a; do     shift;     for b; do         printf "%s %s %s %s\n" "$a" "$b" "$a" "$b";     done; done > pairs_x2.txt

##--Make code for subtracting matrices to get pairwise difference matrices:
	sed 's/\(^[0-9]*\)\ /\1/' pairs_x2.txt > pairs1.tmp;
	sed 's/^/dif\_m/g' pairs1.tmp > pairs2.tmp;
	sed 's/\ \([0-9]*\)\ \([0-9*]\)/\ \<\-\ m\1\ \-\ m\2\;\ /' pairs2.tmp > pairs3.tmp;
	MAKE_PAIR_DIFF_MATRICES=$(cat ./pairs3.tmp);	## THIS goes into the Rscript section (STEP #4 below) as a one-liner, using only the name of the environmental variable. All text stored within this variable will be echoed into the script by the shell, and then the script is run at the end.

##--Make code for summing all difference matrices:
	sed 's/\([0-9]*\)\ /dif_m\1/; s/$/\ /' ./pairs.txt | sed 's/\,\ $'$CR'^$//' > pairs4.tmp;
	MY_PAIRS=$(cat ./pairs4.tmp);
	SUM_ALL_DIFF_MATRICES=$(echo $MY_PAIRS | sed 's/\ /\,\ /g');	## THIS goes into the Rscript section as 'sum(c($SUM_ALL_DIFF_MATRICES))' to sum all difference matrices, for overall comparison.

##--Make code for summing original tij matrices:
	MATRIX_NUMBERS=$(echo $(seq 6) | sed 's/\([0-9*]\)/m\1/g');
	SUM_TIJ_MATRICES=$(echo $MATRIX_NUMBERS | sed 's/\ /\,\ /g');	## THIS goes into the Rscript section as 'sum(c($SUM_TIJ_MATRICES))'

##--Make code for calculating difference proportion of total (from original two tij matrices)
##--for each pairwise difference matrix. We will make the R code by using sed to modify the 
##--pairs of combinations present in the "pairs.txt" file. This operation is more complex 
##--than the simple sed commands in the preceding lines, and though I'd like to use comment
##--symbols (pound signs) in the code the shell doesn't handle those so well. So we will 
##--write the code to an Rscript, which we will then source from within the Rscript below,
##--using the 'source' function. First, let's write the code into a txt file:
	sed 's/\([0-9]*\)\ \([0-9]*\)/'$CR'x\ \<\-\ sum(dif_m\1\2)'$CR'y\ \<\-\ sum(m\1\,\ m\2)'$CR'diffpro_m\1\2\ \<\-\ x\/y'$CR'diffpro_m\1\2'$CR''$CR'/g' ./pairs.txt > ./diffProTests.txt ;

##--Now we have modified the pairs into a set of R code in a local txt file, so we need to 
##--make an appropriate header for the Rscript (I will just echo it) and then concatenate 
##--the header and the txt code into a '.r' file that R can call:
echo "
#!/usr/bin/env Rscript
" > Rtop.tmp

	cat ./Rtop.tmp ./diffProTests.txt > ./diffProTests.r ;

#
##--Make code for logical testing of whether tij matrices are identical in their portions
##--containing known observations, or not, etc. We will make this code based on the total
##--number of rows with known observations in the :
	cd ./run0*/txt/
	NUM_KNOWN_ROWS=$(grep "\"known" ./mydata_names_df.txt | wc -l)
	cd ..; cd ..;

	##--Make txt file with the logical testing code and place in variable and in Rscript:
	sed 's/\([0-9]*\)\ \([0-9]*\)/identical(tijObjList\[\1\]\[\[1\]\]\[1\:'"$NUM_KNOWN_ROWS"'\,\]\,\ tijObjList\[\2\]\[\[1\]\]\[1\:'"$NUM_KNOWN_ROWS"'\,\])/g; s/\ //g' ./pairs.txt > ./logicalTests.txt ;
	IDENT_MATRIX_TESTS=$(cat ./logicalTests.txt | sed 's/)/)\;/g'); 	## Place logical testing code into environmental variable.
	cat ./Rtop.tmp ./logicalTests.txt > ./logicalTests.r ;			## Make final Rscript for this code, which you will source in the compareProbs script below.
	rm ./Rtop.tmp;

##--Make code for combining all difference proportion of total ('diffpro_') calculations
##--in order to plot them:
	COMBINE_DIFF_PROPS=$(echo $MY_PAIRS | sed 's/\ /\,\ /g; s/dif/diffpro/g');	## THIS goes into the Rscript section as 'sum(c($SUM_ALL_DIFF_MATRICES))' to sum all difference matrices, for overall comparison.





echo "INFO      | $(date) | STEP #4: MAKE R SCRIPT CONTAINING ENVIRONMENTAL VARIABLES AND ANALYSIS CODE FOR VISUALIZING AND COMPARING RESULTS... "
############ MAKE R SCRIPT
echo "
#!/usr/bin/env Rscript

################################### compareProbs.R #######################################

############ I. SETUP
setwd('$MY_PATH')
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
	assign(basename, as.matrix(read.table(tijFileNames[i], header=TRUE, sep='\t')))
}

############ III. PLOT THE MATRICES IN A SINGLE GRAPH TO PERMIT BROAD VISUAL COMPARISON OF RESULTS ACROSS RANGE OF TEST VALUES
##--Plot the matrices in the quartz window with par-specified parameters (given above):
tijObjList <- list($FUTURE_TIJ_OBJ_NAMES)
pdf('plots.pdf')

##--Use par and col to specify the plot layout and color characteristics:
par(mfrow=c(3,2), mai = c(0.5, 0.1, 0.1, 0.1))
col<- colorRampPalette(c('ivory2', 'yellow', 'red'))(20)

for(i in tijObjList){
	image(f(i), col=col)
}
dev.off()



##--Some logical testing to see if any matrices are the same (likely, they are NOT):
$IDENT_MATRIX_TESTS


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
$MAKE_PAIR_DIFF_MATRICES
#
##--Next, get the sum of all differences across the pairwise diff matrices, as well as the
##--sum of all tij matrices themselves. Then compare these two to do the 'test' #2 mentioned
##--above under 'Some mathematical/arithmetic testing', by simple division:
sum_difm <- sum(c($SUM_ALL_DIFF_MATRICES))
sum_difm
sum_tijm <- sum(c($SUM_TIJ_MATRICES))
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
	cat('---', diffpro_list[i][[1]][1],'---\n');
	print(get(diffpro_list[[i]][1]))
}
diffpros <- abs(c($COMBINE_DIFF_PROPS))

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
" > comparePPMats.r


############ FINAL STEPS:
echo "INFO      | $(date) | STEP #5: RUN THE R SCRIPT AND CONDUCT CLEANUP. "
	R CMD BATCH ./comparePPMats.R ;

##--Cleanup remaining unnecessary files:
	rm ./tij*.txt;
	rm ./*.tmp;
	rm ./pairs*.txt;
	rm ./diffProTests.txt ./logicalTests.txt;

echo "INFO      | $(date) | Done conducting sensitivity test(s) examining the effect of varying the 'prior' probabilities of known \
observations in the beliefs matrix supplied to bgmm, using bgmmSensTest."
echo "INFO      | $(date) | Bye.
"
#
#
#
######################################### END ############################################

exit 0
