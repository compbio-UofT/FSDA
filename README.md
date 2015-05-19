# FSDA
##

## Pre-process and reference files creation
The files in the build_ref and preprocess directories are used for creating the Mixture of Gaussians parameters and other reference files.

## Prediction
The script predict_region.sh in the directory prediction can be used to predict the copy number of a region. The test sample and the reference sample locations and the test sample's admixture ratio are set inside the script. The script usage is as follows:

./predict_region.sh [chr] [begin] [end] [ks_threshold] [neighbours_threshold] [admixture] [cov] [cnv]

The last three arguments are for simulation purposes. The valid values fo the "cnv" argument are "nor", "dup".

example:
./predict_region.sh chr3 10000000 11000000 0.002 5 0.13 40 dup

