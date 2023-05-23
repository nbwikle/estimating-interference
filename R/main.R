### main.R
### author: Nathan Wikle.
###

### Recreates the analysis found in our manuscript, "Causal health impacts of 
###   power plant emission controls under modeled and uncertain physical process 
###   interference." This script serves to 1) cleans/processes the data, 2) fit the
###   process and outcome models, 3) perform the simulation study, and 4) reproduce
###   relevant results and figures found in the manuscript. However, if run 
###   sequentially, it will take a very long time (X hours on 2.3 GHz
###   Quad-Core Intel Core i7 processor). It is recommended that components of 
###   X and Y be run in parallel.

### !! CAUTION: THIS R SCRIPT WILL TAKE A LONG TIME !!! ###

### 1. Required packages (version number) ###

library(colortools)       # 0.1.5
library(data.table)       # 1.13.2
library(dplyr)            # 1.0.2
library(forcats)          # 0.5.0
library(geosphere)        # 1.5.10
library(ggplot2)          # 3.3.2
library(ggridges)         # 0.5.2
library(here)             # 1.0.1
library(hrbrthemes)       # 0.8.0
library(inborutils)       # 0.1.0.9086
library(INLA)             # 19.9.3
library(maps)             # 3.3.0
library(Matrix)           # 1.2.18
library(matrixStats)      # 0.57.0
library(mvnfast)          # 0.2.5.1
library(ncdf4)            # 1.17
library(raster)           # 3.4.5
library(rgdal)            # 1.5.18
library(rwc)              # 1.11
library(sp)               # 1.4.4
library(stringr)          # 1.4.0
library(TruncatedNormal)  # 2.2
library(viridis)          # 0.5.1

### 2. Create subdirectories ###

# create a new folder for data
dir.create(here::here("data"), showWarnings = FALSE)
# create a new folder for output
dir.create(here::here("output"), showWarnings = FALSE)

### 3. Download data ###

# data used in analysis, size = X GB, doi = ZENODO.DOI
download_zenodo(
  doi = "ZENODO.DOI",
  path = here::here("data"), parallel = FALSE, quiet = FALSE
)
