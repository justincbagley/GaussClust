# GaussClust
Clustering using Gaussian mixture modeling for species delimitation and classification

## LICENSE

All code within the GaussClust repository is available "AS IS" under a generous GNU license. See the [LICENSE](LICENSE) file for more information.

## CITATION

If you use scripts from this repository as part of your published research, I require that you cite the repository as follows (also see DOI information below): 
  
- Bagley, J.C. 2016. GaussClust. GitHub repository, Available at: http://github.com/justincbagley/GaussClust.

Alternatively, please provide the following link to this software repository in your manuscript:

- https://github.com/justincbagley/GaussClust

## DOI

The DOI for GaussClust, via [Zenodo](https://zenodo.org), is coming soon.

## INTRODUCTION

*"Clustering and discriminant analysis (or classification) methods are among the most important
techniques in multivariate statistical learning." - Lebret et al. (2015)*

This repository focuses on ways to use Gaussian Mixture Models (GMMs) to conduct clustering analyses that address problems in systematics, particularly in the delimitation of species, and classification of individuals to species, using univariate or multivariate data. Recent papers discuss the promise of such methods for species delimitation, with or without multiple data-types (e.g. genetic, morphological, ecological data), and with or without accounting for noise during clustering (Hausdorf & Hennig 2010; Edwards & Knowles 2014). The current pre-release version of the repository focuses on GaussClust, a shell script that works with R to conduct various GMM analyses. 

As noted by Lebret et al. (2015), two foci of multivariate approaches related to clustering are (1) clustering proper, which aims to group observations (e.g. individuals) into groups or 'clusters' that are more similar to one another than to other clusters, and (2) classification methods where a discriminant function is used to assign new data to groups that are known a priori. GaussClust primarily focuses on the former, but the latter also sneaks into clustering analyses and therefore is also included. As in the case of the author's other software on GitHub (e.g. [PIrANHA](https://github.com/justincbagley/PIrANHA)), GaussClust is fully command line-based and is available as open-source software according to the license. 

## GETTING STARTED

### Dependencies
![alt tag](https://cloud.githubusercontent.com/assets/10584087/21243724/b12a94d2-c2df-11e6-9d20-5cf06877ad94.png)
- The [R](https://www.r-project.org) software environment (available from download mirrors from The Comprehensive R Archive Network (CRAN) such as https://cloud.r-project.org), as well as several R packages including:
- [bgmm](https://cran.r-project.org/web/packages/bgmm/index.html)
- [Rmixmod](https://cran.r-project.org/web/packages/Rmixmod/index.html)
- [StatMatch](https://cran.r-project.org/web/packages/StatMatch/index.html)
- [MASS](https://cran.r-project.org/web/packages/MASS/index.html)
- [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html)
- [labdsv](https://cran.r-project.org/web/packages/labdsv/index.html)
[More info coming soon.]

# USAGE
````
## Assuming GaussClust.sh is present in the current working directory, display usage for the script by simply entering its name at the command line:
$ ./GaussClust.sh

##########################################################################################
#                             GaussClust v1.0, December 2016                             #
##########################################################################################

INFO      | Thu Dec 15 23:50:28 CST 2016 | STEP #1: SETUP. SETTING OPTIONS AND PATH VARIABLE... 

Usage: ./GaussClust.sh [options] inputFile
  
Options: -k nmdsDimensions (specify number of dimensions, k, to retain during NMDS on Gower distances) | -u unsuperGMM (0=no unsupervised GMM is carried out; 1=conduct unsupervised GMM using 'Rmixmod' R pacakge, for comparative or individual purposes) | -r rangeNumClusters (optional numeric listing of a range, x:y, of the number of clusters to be modeled over during unsupervised GMM in Rmixmod) | -d GMMDiscrim (0=(semi-)supervised, GMM-based discriminant analysis is not carried out with 'mixmodelLearn' in Rmixmod; 1=conduct discriminant analysis) | -b beliefBasedMM (0=no mixture modeling is carried out using the 'bgmm' R package; 1=calls 'supervised' GMM analysis, 2=calls 'semisupervised' GMM analysis, and 3=calls both supervised and semisupervised analyses in bgmm) | -p probsMatrix (specify matrix of plausibilities, or weights of the prior probabilities for labeled observations, for bgmm; also sets belief matrix) | -c numComponents (specify number of components (e.g. Gaussian components) or 'clusters' to assign individuals to during regular GMM (single value, rather than a range; see -r above) or bgmm modeling)

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
regular GMM or bgmm modeling (except see other option available using -n flag above). This 
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

````

## Real-world example \#1 with screen output:
Below is an example of the usage for GaussClust (shown above), using the example data files ("probs_35.txt" and "Enyalius_35.txt") located in the "examples" folder of the distro. In this simplified example, we use GaussClust to conduct only unsupervised Gaussian clustering on 6 clusters, on 4 dimensions of data retained from NMDS. All screen (Terminal) output from the analysis is shown. When you run GaussClust, 'INFO' and date printed to screen along with info for each of the broader steps of the script (details from R go to an *.Rout file). Question-response lines output to screen are marked 'FLOW'. As you can see, the code runs with no error messages.
````
$ ./GaussClust.sh -k 4 -u 1 -r 0 -d 0 -b 0 -p ./probs_35.txt -c 2 ./Enyalius_35.txt

##########################################################################################
#                             GaussClust v1.0, December 2016                             #
##########################################################################################

INFO      | Thu Dec 15 23:24:00 CST 2016 | STEP #1: SETUP. SETTING OPTIONS AND PATH VARIABLE... 
INFO      | Thu Dec 15 23:24:00 CST 2016 | STEP #2: MAKE GAUSSIAN CLUSTERING R SCRIPT CONTAINING ENVIRONMENTAL VARIABLES AND ANALYSIS CODE... 
INFO      | Thu Dec 15 23:24:00 CST 2016 | STEP #3: RUN THE R SCRIPT. 
INFO      | Thu Dec 15 23:24:03 CST 2016 | STEP #4: CLEAN UP THE WORKSPACE. 
INFO      | Thu Dec 15 23:24:03 CST 2016 |          Moving R ouput files to new folder named 'R'... 
FLOW      | Thu Dec 15 23:24:03 CST 2016 |          Would you like to keep the Rscript output by GaussClust? (y/n) : y
FLOW      | Thu Dec 15 23:24:06 CST 2016 |          Would you like to keep text files output by GaussClust? (y/n) : y
INFO      | Thu Dec 15 23:24:07 CST 2016 | Done conducting Gaussian clustering and related analyses using GaussClust.
INFO      | Thu Dec 15 23:24:07 CST 2016 | Bye.

````

## Real-world example \#2:
Here, (1) we specify to keep 4 NMDS dimensions (-k 4); (2) conduct a regular unsupervised GMM analysis in Rmixmod, using multiple models across 5-20 clusters (-u 1 -r 5:20), which are compared to identify the best model using BIC; (3) call supervised discriminant analysis in Rmixmod (-d 1); (4) attempt supervised and semisupervised analysis in bgmm (-b 3 for 'both'); and (4) specify that analyses use or require 13 clusters (as needed; -c 13). Output is not redirected (e.g. by placing something like '> output.txt' at the end, so all output from the script is printed to screen (except for steps conducted in R).
````
./GaussClust.sh -k 4 -u 1 -r 5:20 -d 1 -b 3 -p B_206.txt -c 13 ./mydata.txt
````

## REFERENCES
- Edwards DL, Knowles LL (2014) Species detection and individual assignment in species delimitation: can integrative data increase efficacy? Proceedings of the Royal Society B, 281, 20132765. 
- Hausdorf B, Hennig C (2010). Species delimitation using dominant and codominant multilocus markers. Systematic Biology, 59, 491-503.
- Lebret R, Iovleff S, Langrognet F, Biernacki C, Celeux G, Govaert G (2015) Rmixmod: the R package of the model-based unsupervised, supervised, and semi-supervised classification Mixmod Library. Journal of Statistical Software, 67(6). doi:10.18637/jss.v067.i06

## TODO
- Solve two input file problem. **DONE!** :white_check_mark:
- Make script do more with bgmm, including semisupervised analysis using belief probs matrix. **DONE!** :white_check_mark:
- Change Usage section to include code for working with example files.  **Partially DONE!**
- Explore different methods for BIC calculation (e.g. in mclust) during unsupervised clustering.
- Explain defaults better.
- Would be nice if not giving a flag meant the corresponding analysis was not run (equivalent to passing '0' for some of the options)

December 15, 2016
Justin C. Bagley, Tuscaloosa, AL, USA
