# Source Code and Supplementary Material

*Authors: Nathan B. Wikle and Corwin M. Zigler*

This repository contains source code and additional supplementary materials from our manuscript, "Causal health impacts of power plant emission controls under modeled and uncertain physical process interference."  Supplementary materials from the manuscript, including a detailed simulation study, can be found within `supp-material.pdf`. Relevant data (doi = `10.5281/zenodo.8015586`) and results (doi = `10.5281/zenodo.8015752`) have been archived for reproducibility.

The following instructions provide details on how to run the source code underlying the analysis. 

## Requirements

The code has been tested with R version 4.1.3, "One Push-Up."  The following **R packages** must be installed before the code will run successfully **(version in parentheses)**:

- [`abind`](https://cran.r-project.org/web/packages/abind/index.html)  (1.4.5)
- [`BART`](https://cran.r-project.org/web/packages/BART/index.html)  (2.9.3) 
- [`CBPS`](https://cran.r-project.org/web/packages/CBPS/index.html)  (0.23)
- [`cobalt`](https://cran.r-project.org/web/packages/cobalt/index.html)  (4.3.2)
- [`copula`](https://cran.r-project.org/web/packages/copula/index.html)  (1.1.0)
- [`data.table`](https://CRAN.R-project.org/package=data.table)  (1.14.2)
- [`data.tree`](https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html)  (1.0.0)
- [`dbarts`](https://cran.r-project.org/web/packages/dbarts/index.html)  (0.9.22)
- [`dplyr`](https://CRAN.R-project.org/package=dplyr)  (1.0.8)
- [`fst`](https://cran.r-project.org/web/packages/fst/index.html)  (0.9.8)
- [`ggplot2`](https://CRAN.R-project.org/package=ggplot2)  (3.3.5)
- [`GIGrvg`](https://cran.r-project.org/web/packages/GIGrvg/index.html)  (0.7)
- [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html)  (4.1.4)
- [`here`](https://CRAN.R-project.org/package=here)  (1.0.1)
- [`inborutils`](https://github.com/inbo/inborutils)  (0.3.0)
- [`maps`](https://cran.r-project.org/web/packages/maps/index.html)  (3.4.0)
- [`MASS`](https://cran.r-project.org/web/packages/MASS/index.html)  (7.3.55)
- [`MatchIt`](https://cran.r-project.org/web/packages/MatchIt/index.html)  (4.4.0)
- [`Matrix`](https://cran.r-project.org/web/packages/Matrix/index.html)  (1.4.0)
- [`matrixStats`](https://cran.rstudio.com/web/packages/matrixStats/index.html)  (0.61.0)
- [`mvtnorm`](https://cran.r-project.org/web/packages/mvtnorm/index.html)  (1.1.3)
- [`ncdf4`](https://cran.r-project.org/web/packages/ncdf4/index.html)  (1.19)
- [`nlme`](https://cran.r-project.org/web/packages/nlme/index.html)  (3.1.155)
- [`nnet`](https://cran.r-project.org/web/packages/nnet/index.html)  (7.3.17)
- [`npreg`](https://cran.r-project.org/web/packages/npreg/index.html)  (1.0.9)
- [`numDeriv`](https://cran.r-project.org/web/packages/numDeriv/index.html)  (2016.8.1.1)
- [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)  (4.1.3)
- [`raster`](https://CRAN.R-project.org/package=raster)  (3.5.15)
- [`rcartocolor`](https://cran.r-project.org/web/packages/rcartocolor/rcartocolor.pdf)  (2.0.0)
- [`Rcereal`](https://cran.r-project.org/web/packages/Rcereal/index.html)  (1.2.1.1)
- [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html)  (1.1.2)
- [`Rcpp`](https://cran.r-project.org/web/packages/Rcpp/index.html)  (1.0.9)
- [`RcppArmadillo`](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)  (0.11.2.0.0)
- [`RcppDist`](https://cran.r-project.org/web/packages/RcppDist/index.html)  (0.1.1)
- [`rstan`](https://cran.r-project.org/web/packages/rstan/index.html)  (2.21.3)
- [`rwc`](https://cran.r-project.org/web/packages/rwc/index.html)  (1.11)
- [`sf`](https://cran.r-project.org/web/packages/sf/index.html)  (1.0.7)
- [`sp`](https://CRAN.R-project.org/package=sp)  (1.4.6)
- [`StanHeaders`](https://cran.r-project.org/web/packages/StanHeaders/index.html)  (2.21.0.7)
- [`stars`](https://cran.r-project.org/web/packages/stars/index.html)  (0.5.5)
- [`stringr`](https://CRAN.R-project.org/package=stringr)  (1.4.0)
- [`survival`](https://cran.r-project.org/web/packages/survival/index.html)  (3.2.13)
- [`tidycensus`](https://cran.r-project.org/web/packages/tidycensus/index.html)  (1.2)
- [`tidyr`](https://cran.r-project.org/web/packages/tidyr/index.html)  (1.2.0)
- [`truncnorm`](https://cran.r-project.org/web/packages/truncnorm/index.html)  (1.0.8)
- [`units`](https://cran.r-project.org/web/packages/units/index.html)  (0.8.0)
- [`USAboundaries`](https://github.com/ropensci/USAboundaries)  (0.4.0)
- [`WeightIt`](https://cran.r-project.org/web/packages/WeightIt/index.html)  (0.13.1)

### `countbart` package install

In addition to the above packages, the `countbart` package must be installed from source. This package implements the log-linear BART models used in our paper. The package tarbal has been included as `./src/countbart_0.1.tar.gz`, and the install is performed in Section 1 of the script `./R/main.R`.

## Data

Publicly available data have been archived for download at the beginning of the analysis (doi = 10.5281/zenodo.8015586). They require 3.4 GB of space within the subdirectory. Unfortunately, due to privacy considerations, the health outcome data (pediatric asthma ED visits and Medicare all-cause mortality in Texas) cannot be publicly distributed. Instead, we have included a synthetic dataset, `synth-ped-asthma-data`, as a placeholder. 

**Eleven main data sources are downloaded:**

1. USA coal-fired power plant emissions data
- `AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv`
- Coal-fired power plan emissions data from the US Environmental Protection Agency's [Air Markets Program Data (AMPD)](https://ampd.epa.gov/ampd/). The AMPD provides access to current and historical monthly emissions data on electricity generating units (EGUs), collected as part of EPA's emissions trading programs. A comprehensive list of AMPD column definitions is included for reference (see the downloaded `AMPD_column_definitions.csv` file).

2. US Census 2016 American Community Survey (ACS) data
- `Census_2016_TxZCTA.RDS`
- Demographic data were downloaded from the [US Census ACS](https://www.census.gov/programs-surveys/acs) using the R package [`tidycensus`](https://walker-data.com/tidycensus/). The R script which generates this data (`./R/outcome-analysis/census-data.R`) has been included for reference.

3. Daymet Annual Climate Summaries
- `daymet_v4_prcp_annttl_na_2016.nc`
- `daymet_v4_tmax_annavg_na_2016.nc`
- `daymet_v4_tmin_annavg_na_2016.nc`
- `daymet_v4_vp_annavg_na_2016.nc`
- High-resolution climate data (annual total precipitation, annual average daily maximum 2-meter air temperature, annual average daily minimum 2-meter air temperature, and annual average daily water vapor pressure), downloaded from ORNL's [Daymet Version 4](https://daac.ornl.gov/DAYMET/guides/Daymet_V4_Annual_Climatology.html).

4. 2016 mean SO4 and black carbon concentrations
- `GWRwSPEC_SO4_NA_201601_201612.nc`
- `GWRwSPEC_BC_NA_201601_201612.nc`
- Annual mean sulfate and black carbon concentrations were obtained from Randall Martin's Athmospheric Composition Analysis Group's [North American Regional Estimates (V4.NA.03) dataset](https://sites.wustl.edu/acag/datasets/surface-pm2-5/#V4.NA.03).

5. HyADS coal-attributed PM2.5 concentrations
- `HyADS_grids_pm25_byunit_2016.fst`
- `HyADS_grids_pm25_total_2016.fst`
- Estimated 2016 coal-attributed PM2.5 concentrations from the Hysplit average dispersion (HyADS) model [(Henneman et al. (2019)](https://doi.org/10.1097/EDE.0000000000001024), a deterministic, reduced-complexity model of pollution transport.

6. Mexico coal-fired power plant emissions data
- `Mexico_2016_point_interpolated_02mar2018_v0.csv`
- 2016 coal-fired power plant SO2 emissions data, as estimated by the US EPA's [National Emissions Inventory platform, 2016v1](https://www.epa.gov/air-emissions-modeling/2014-2016-version-7-air-emissions-modeling-platforms).

7. North American Regional Reanalysis meteorological data
- `rhum.2m.mon.mean.nc`
- `uwnd.10m.mon.mean.nc`
- `vwnd.10m.mon.mean.nc`
- Meteorological data (relative humidity, 10m wind velocities) from the NOAA Physical Science Laboratory [NCEP North American Regional Reanalysis (NARR)](https://psl.noaa.gov/data/gridded/data.narr.monolevel.html) database.

8. Cigarette smoking data
- `smokedatwithfips_1996-2012.csv`
- Data on cigarette smoking prevalance in US counties from 1996-2012, obtained from [Dwyer-Lindgren et al. (2014)](https://doi.org/10.1186/1478-7954-12-5).

9. Synthetic pediatric asthma data
- `synth-ped-asthma-data.csv`
- *Simulated* pediatric asthma data to match the format, but not the observations, from the [Texas Health Care Information Collection (THCIC), Texas DSHS](https://www.dshs.texas.gov/texas-health-care-information-collection).

10. Texas state shape file
- `texas-state-sf.RDS`
- A shapefile of the state of Texas obtained from [US Census](https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html)

11. US ZIPs-to-county data crosswalk
- `tx-zip-to-county.csv`
- A data crosswalk matching Texas ZIP codes to counties, obtained from the [Missouri Census Data Center](mcdc.missouri.edu/applications/geocorr2014.html).

## Instructions

Before running any code, make sure the required R packages have been installed. Open and run the `main.R` file, found in the `./R/` folder.  **Note that this script will take a very long time to run sequentially**. On a 2.3 GHz Quad-Core Intel Core i7 processor with 32 GB of memory, the script would take approximately 3110 hours to run. It is **highly recommended** that Sections 5 and 6 be run in parallel; by default, they will run in parallel on POSIX systems, with the number of workers = `detectCores()`. 

The script has been divided into seven sequential steps, which perform the following:

### Step One: 

- Loads required packages into R.
- Install `countbart` from source.
- Time: < 5 minutes.

### Step Two: 

- Creates `./data/` and `./output/` subdirectories, which will hold the underlying data and analysis output, respectively.
- Time: < 0.01 seconds.

### Step Three:

- Downloads the raw data sources used in the analysis. These data are publicly available, and have been archived (doi = 10.5281/zenodo.8015586) for reproducibility. These data require 3.4 GB of space.
- Downloads the main results from the analyses and simulation studies; they have been archived at Zenodo (doi = 10.5281/zenodo.8015752). The results require 28.9 GB of space.
- Time: < 1 hour.
- Output size: 32.3 GB.

### Step Four: 

- Performs the **sulfate analysis** from the main text, including: (i) cleaning and processing the emissions, wind, and sulfate data, and (ii) sampling from the sulfate model's posterior distribution via MCMC. 
- Time: ~110 hours.
- **Note: `./R/sulfate-analysis/so4-model.R` can be skipped, as the main results have already been downloaded! This saves ~ 110 computing hours.**

### Step Five: 

- Performs the **outcome model analysis**. This includes (i) cleaning and processing the census, climate, facility, black carbon, and outcome data, and (ii) sampling from the outcome model's posterior distbitution via MCMC. In particular, four types of inference are considered: 1. log-linear BART regression *with* uncertainty in the interference struxture (cut), 2. log-linear BART regression *without* uncertainty propagation (plugin), 3. Poisson regression with cut, and 4. Poisson regression with plugin. 
- Time: ~ 48 hours / # workers
- **Note: `./R/outcome-analysis/asthma-pois-plugin.R`, `./R/outcome-analysis/asthma-pois-cut.R`, `./R/outcome-analysis/asthma-bart-plugin.R`, and `./R/outcome-analysis/asthma-bart-cut.R` can be skipped, as the main results have already been downloaded! This saves ~ 47 computing hours.**

### Step Seven: 

- Performs the **simulation studies** and **sensitivity analysis** found in the Supplement. This includes (i) creating data structures used throughout the simulation studies, (ii) six simulation studies, (iii) log-linear BART prior sensitivity analysis, and (iv) plotting the results of the simulation studies and sensitivity analysis (found in `./output/supplement/`).
- Time: ~ 2952 hours / # workers.
- **Note: The code which performs the simulation studies and sensitivity analysis can be skipped, as the main results have already been downloaded! This saves ~ 2951 computing hours.**

### Step Seven: 

- - Generates the Figures and results found in the main (`./output/main`) and supplementary (`./output/supplement`) text. 
- Time: < 15 minutes.

*Note: If the recommended steps have been skipped, the total runtime of `./main.R` is ~ 2.5 hours.*

## Output

Upon successful completion of `main.R`, the results and figures found in the main text and supplement are saved as PNGs in the `./output/main` and `./output/supplement` folders. These will require ~ 80 MB of space. 
