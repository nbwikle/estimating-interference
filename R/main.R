### main.R
### author: Nathan Wikle.
###

### Recreates the analysis found in our manuscript, "Causal health impacts of 
###   power plant emission controls under modeled and uncertain physical process 
###   interference." This script serves to 1) cleans/processes the data, 2) fit the
###   process and outcome models, 3) perform the simulation study, and 4) reproduce
###   relevant results and figures found in the manuscript. However, if run 
###   sequentially, it will take a very long time (over 3110 hours! on a 2.3 GHz
###   Quad-Core Intel Core i7 processor). It is recommended that Sections 5
###   and 6 be run in parallel.

### !! CAUTION: THIS R SCRIPT WILL TAKE A LONG TIME !!! ###
###
###  -> It is strongly recommended that Sections 5 (outcome analysis)
###     and 6 (Simulation Study) be run in parallel. The script has been    
###     set up to run in parallel on POSIX systems by default, with 
###     'workers = detectCores()'.                   
###

### Storage Requirements
###
###  -> Note: This analysis will require ~32.5 GB of storage. 
###     The majority of this will be from the initial raw data
###     download (3.4 GB) and the results download (28.9 GB) in
###     Section 3.

##################################################################
### 1. Required packages (version number)
##################################################################

library(abind)          # 1.4.5
library(BART)           # 2.9.3
library(CBPS)           # 0.23
library(cobalt)         # 4.3.2
library(copula)         # 1.1.0
library(data.table)     # 1.14.2
library(data.tree)      # 1.0.0
library(dbarts)         # 0.9.22
library(dplyr)          # 1.0.8
library(fst)            # 0.9.8
library(ggplot2)        # 3.3.5
library(GIGrvg)         # 0.7
library(glmnet)         # 4.1.4
library(here)           # 1.0.1
library(inborutils)     # 0.3.0
library(maps)           # 3.4.0
library(MASS)           # 7.3.55
library(MatchIt)        # 4.4.0
library(Matrix)         # 1.4.0
library(matrixStats)    # 0.61.0
library(mvtnorm)        # 1.1.3
library(ncdf4)          # 1.19
library(nlme)           # 3.1.155
library(nnet)           # 7.3.17
library(npreg)          # 1.0.9
library(numDeriv)       # 2016.8.1.1
library(parallel)       # 4.1.3
library(raster)         # 3.5.15
library(rcartocolor)    # 2.0.0
library(Rcereal)        # 1.2.1.1
library(RColorBrewer)   # 1.1.2
library(Rcpp)           # 1.0.9
library(RcppArmadillo)  # 0.11.2.0.0
library(RcppDist)       # 0.1.1
library(rstan)          # 2.21.3
library(rwc)            # 1.11
library(sf)             # 1.0.7
library(sp)             # 1.4.6
library(StanHeaders)    # 2.21.0.7
library(stars)          # 0.5.5
library(stringr)        # 1.4.0
library(survival)       # 3.2.13
library(tidycensus)     # 1.2
library(tidyr)          # 1.2.0
library(truncnorm)      # 1.0.8
library(units)          # 0.8.0
library(USAboundaries)  # 0.4.0
library(WeightIt)       # 0.13.1

### install 'countbart' from source

# install
install.packages(
  here::here("src", "countbart_0.1.tar.gz"), 
  repos = NULL, 
  type = "source"
)

# load package
library(countbart)


##################################################################
### 2. Create subdirectories 
##################################################################

# create a new folder for data
dir.create(here::here("data"), showWarnings = FALSE)
# create a new folder for output
dir.create(here::here("output"), showWarnings = FALSE)


##################################################################
### 3. Download data and results
##################################################################

# data used in analysis, size = 3.4 GB, doi = 10.5281/zenodo.8015586
download_zenodo(
  doi = "10.5281/zenodo.8015586",
  path = here::here("data"), parallel = FALSE, quiet = FALSE
)

# main results from analysis, size = 28.9 GB, doi = 10.5281/zenodo.8015752
download_zenodo(
  doi = "10.5281/zenodo.8015752",
  path = here::here("output"), parallel = FALSE, quiet = FALSE
)


##################################################################
### 4. Sulfate analysis 
###     -> Total time: ~110 hours 
##################################################################

# i. so4 functions
source(here::here("R", "sulfate-analysis", "so4-functions.R"))

# ii. create facility data from AMPD emissions database  (~40 seconds)
source(here::here("R", "sulfate-analysis", "make-facility-data.R"))

# iii. create SO4 raster data  (~30 seconds)
source(here::here("R", "sulfate-analysis", "so4-data.R"))

# iv. sample from sulfate model posterior with MCMC  (~110 hours)
#       !!! CAUTION: TAKES ~110 hours !!!
source(here::here("R", "sulfate-analysis", "so4-model.R"))


##################################################################
### 5. Health outcome analysis 
###     -> Total time: ~48 hours / # workers
##################################################################

# i. outcome functions
source(here::here("R", "outcome-analysis", "outcome-functions.R"))
Rcpp::sourceCpp(here::here("src", "pois-reg.cpp"))

## ii. query US Census database
       
# Note: This code downloads US Census data from the 2016 American 
#   Community Survey (ACS) using an API key provided by the US Census
#   Bureau. You can request an API key here: 
#       https://api.census.gov/data/key_signup.html
#   Once requested, paste the key into the appropriate variable at the 
#   top of 'census-data.R'.
#   
#   For convenience, the output of this script has been included in the 
#   initial data download ("./data/Census_2016_TxZCTA.RDS") and does not
#   need to be repeated.
# source(here::here("R", "outcome-analysis", "census-data.R"))

### iii. process/clean the covariate and outcome data
#     -> time: ~ 1hr 10min
source(here::here("R", "outcome-analysis", "smoking-data.R"))
source(here::here("R", "outcome-analysis", "climate-data.R")) 
source(here::here("R", "outcome-analysis", "facility-data.R"))
source(here::here("R", "outcome-analysis", "outcome-data.R"))
source(here::here("R", "outcome-analysis", "analysis-setup.R")) 
source(here::here("R", "outcome-analysis", "bc-data.R"))

### iv. asthma analysis 

# poisson regression, plug-in inference 
#   -> time: ~5 minutes
source(here::here("R", "outcome-analysis", "asthma-pois-plugin.R"))

# poisson regression, with uncertainty propagation
#   -> time: ~6.25 hours / # workers (default: workers = detectCores())
source(here::here("R", "outcome-analysis", "asthma-pois-cut.R"))

# log-linear BART regression, plug-in inference
#   -> time: ~50 minutes
source(here::here("R", "outcome-analysis", "asthma-bart-plugin.R"))

# log-linear BART regression, with uncertainty propagation
#   -> time: ~41.5 hours / # workers (default: workers = detectCores())
source(here::here("R", "outcome-analysis", "asthma-bart-cut.R"))


##################################################################
### 6. Simulation study 
###     -> Total time: ~2952 hours / # workers
##################################################################

### i. initial simulation study set-up

# source functions
source(here::here("R", "sim-study", "sim-study-functions.R"))
Rcpp::sourceCpp(here::here("src", "pois-reg.cpp"))

# compile 'stan' LM code and create copula for simulation study
source(here::here("R", "sim-study", "sim-study-setup.R"))

### ii. Continuous Outcome Models (CM)

# CM1 with BART outcome model
#   -> time: ~222 hrs / # workers
source(here::here("R", "sim-study", "CM1-BART.R"))
# CM1 with Bayesian linear regression outcome model
#   -> time: ~42 hrs / # workers
source(here::here("R", "sim-study", "CM1-lm.R"))

# CM2 with BART outcome model
#   -> time: ~230 hrs / # workers
source(here::here("R", "sim-study", "CM2-BART.R"))
# CM2 with Bayesian linear regression outcome model
#   -> time: ~42 hrs / # workers
source(here::here("R", "sim-study", "CM2-lm.R"))

# CM3 with BART outcome model
#   -> time: ~250 hrs / # workers
source(here::here("R", "sim-study", "CM3-BART.R"))
# CM3 with Bayesian linear regression outcome model
#   -> time: ~17 hrs / # workers
source(here::here("R", "sim-study", "CM3-lm.R"))

### iii. Poisson Outcome Models (PM)

# PM1 with log-linear BART outcome model
#   -> time: ~600 hrs / # workers
source(here::here("R", "sim-study", "PM1-BART.R"))
# PM1 with Bayesian Poisson regression outcome model
#   -> time: ~90 hrs / # workers
source(here::here("R", "sim-study", "PM1-pois.R"))

# PM2 with log-linear BART outcome model
#   -> time: ~600 hrs / # workers
source(here::here("R", "sim-study", "PM2-BART.R"))
# PM2 with Bayesian Poisson regression outcome model
#   -> time: ~90 hrs / # workers
source(here::here("R", "sim-study", "PM2-pois.R"))

# PM2 with log-linear BART outcome model
#   -> time: ~600 hrs / # workers
source(here::here("R", "sim-study", "PM3-BART.R"))
# PM2 with Bayesian Poisson regression outcome model
#   -> time: ~96 hrs / # workers
source(here::here("R", "sim-study", "PM3-pois.R"))

### iv. Log-linear BART sensitivity analysis

# Sensitivity to the number of trees.
#   -> time: ~20.5 hours / # workers
source(here::here("R", "sim-study", "sensitivity-ntrees.R"))

# Sensitivity to tree depth.
#   -> time: ~52 hours / # workers
source(here::here("R", "sim-study", "sensitivity-power.R"))

### v. Process simulation study results
#   -> creates list with bias, coverage, and variance results from
#             each simulation study
#   -> recreates individual simulation study plots, as seen in the supplement
#   -> recreates sensitivity analysis plots, as seen in the supplement
#   -> time: ~15 minutes
source(here::here("R", "sim-study", "sim-study-results.R"))
source(here::here("R", "sim-study", "sensitivity-results.R"))


##################################################################
### 7. Results
###     -> Total time: ~12 minutes
##################################################################

### i. recreate figures from main text
#   -> time: ~10 minutes
source(here::here("R", "results", "main-functions.R"))
source(here::here("R", "results", "main-results.R"))

### ii. recreate results and figures from supplement
#   -> time: ~2 minutes
source(here::here("R", "results", "supplement-functions.R"))
source(here::here("R", "results", "supplement-results.R"))
