# GaussClust
Clustering using Gaussian mixture modeling for species delimitation and classification

# Usage
````
## Assuming GaussClust.sh is present in the current working directory, display usage for the script by simply entering its name at the command line:
./GaussClust.sh
````

## Real-world example #1:
Here, (1) we specify to keep 4 NMDS dimensions; (2) conduct a regular unsupervised GMM analysis in Rmixmod, using multiple models across 5-20 clusters, which are compared to identify the best model using BIC; (3) call supervised GMM analysis in Rmixmod; (4) attempt semisupervised analysis in bgmm; and (4) specify that analyses (where needed) specify 15 clusters. Output is not redirected (e.g. '> output.txt' at the end, so all output from the script (but NOT from parts conducted in R) are output to screen.
````
./GaussClust.sh -d 4 -r 1 -n 5:20 -s 1 -b semisupervised -p B_206.txt -c 15 ./mydata_names.txt
````
