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

INFO      | Wed Dec 14 09:03:11 CST 2016 | STEP #1: SETUP. SETTING OPTIONS AND PATH VARIABLE... 

Usage: ./GaussClust.sh [options] inputFile
  
Options: -d numNMDSDimensions (specify number of dimensions, k, to retain during NMDS on Gower distances) | -r regGMM (0=no unsupervised GMM is carried out; 1=conduct unsupervised GMM using 'Rmixmod' R pacakge, for comparative or individual purposes) | -n numGMMClusters (optional numeric listing of a range, x:y, of the number of clusters to be modeled over during unsupervised GMM in Rmixmod) | -s superGMM (0=no (semi-)supervised GMM is carried out in Rmixmod; 1=conduct (semi-)supervised GMM in Rmixmod) | -b beliefBasedMM (0=no mixture modeling is carried out using the 'bgmm' R package; the following other options conduct different kinds of mixture modeling using separate functions available in bgmm: belief, soft, semisupervised, supervised) | -c numComponents (specify number of components (e.g. Gaussian components) or 'clusters' to assign individuals to during regular GMM (single value, rather than a range; see -n above) or bgmm modeling) 

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
````

## Real-world example #1:
Here, (1) we specify to keep 4 NMDS dimensions; (2) conduct a regular unsupervised GMM analysis in Rmixmod, using multiple models across 5-20 clusters, which are compared to identify the best model using BIC; (3) call supervised GMM analysis in Rmixmod; (4) attempt semisupervised analysis in bgmm; and (4) specify that analyses (where needed) specify 15 clusters. Output is not redirected (e.g. '> output.txt' at the end, so all output from the script (but NOT from parts conducted in R) are output to screen.
````
./GaussClust.sh -d 4 -r 1 -n 5:20 -s 1 -b semisupervised -p B_206.txt -c 15 ./mydata_names.txt
````

## REFERENCES
- Edwards DL, Knowles LL (2014) Species detection and individual assignment in species delimitation: can integrative data increase efficacy? Proceedings of the Royal Society B, 281, 20132765. 
- Hausdorf B, Hennig C (2010). Species delimitation using dominant and codominant multilocus markers. Systematic Biology, 59, 491-503.
- Lebret R, Iovleff S, Langrognet F, Biernacki C, Celeux G, Govaert G (2015) Rmixmod: the R package of the model-based unsupervised, supervised, and semi-supervised classification Mixmod Library. Journal of Statistical Software, 67(6). doi:10.18637/jss.v067.i06

## TODO
- Fix two input file problem (get code down to single file).
- Make script do more with bgmm, including semisupervised analysis.
- Change Usage section to include code for working with example files.

December 14, 2016
Justin C. Bagley, Tuscaloosa, AL, USA
