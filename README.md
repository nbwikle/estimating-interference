# Source Code and Supplementary Material

*Authors: Nathan B. Wikle and Corwin M. Zigler*

This repository contains source code and additional supplementary materials from our manuscript, "Causal health impacts of power plant emission controls under modeled and uncertain physical process interference."  Supplementary materials from the manuscript, including a detailed simulation study, can be found within `supp-material.pdf`. Relevant data have been archived (`doi = ZENODO`) for reproducibility.

The following instructions provide details on how to run the source code underlying the analysis. 

## Requirements

The code has been tested with R version 4.1.3, "One Push-Up."  The following **R packages** must be installed before the code will run successfully **(version in parentheses)**:

- [`colortools`](https://CRAN.R-project.org/package=colortools)  (0.1.5)
- [`data.table`](https://CRAN.R-project.org/package=data.table)  (1.13.2)
- [`dplyr`](https://CRAN.R-project.org/package=dplyr)  (1.0.2)
- [`forcats`](https://CRAN.R-project.org/package=forcats)  (0.5.0)
- [`geosphere`](https://CRAN.R-project.org/package=geosphere)  (1.5.10)
- [`ggplot2`](https://CRAN.R-project.org/package=ggplot2)  (3.3.2)
- [`ggridges`](https://CRAN.R-project.org/package=ggridges)  (0.5.2)
- [`here`](https://CRAN.R-project.org/package=here)  (1.0.1)
- [`hrbrthemes`](https://CRAN.R-project.org/package=hrbrthemes)  (0.8.0)
- [`inborutils`](https://github.com/inbo/inborutils)  (0.1.0.9086)
- [`INLA`](https://www.r-inla.org/home)  (19.9.3)
- [`maps`](https://CRAN.R-project.org/package=maps)  (3.3.0)
- [`Matrix`](https://CRAN.R-project.org/package=Matrix)  (1.2.18)
- [`matrixStats`](https://CRAN.R-project.org/package=matrixStats)  (0.57.0)
- [`mvnfast`](https://CRAN.R-project.org/package=mvnfast)  (0.2.5.1)
- [`ncdf4`](https://CRAN.R-project.org/package=ncdf4)  (1.17)
- [`raster`](https://CRAN.R-project.org/package=raster)  (3.4.5)
- [`rgdal`](https://CRAN.R-project.org/package=rgdal)  (1.5.18)
- [`rwc`](https://CRAN.R-project.org/package=rwc)  (1.11)
- [`sp`](https://CRAN.R-project.org/package=sp)  (1.4.4)
- [`stringr`](https://CRAN.R-project.org/package=stringr)  (1.4.0)
- [`TruncatedNormal`](https://CRAN.R-project.org/package=TruncatedNormal)  (2.2)
- [`viridis`](https://CRAN.R-project.org/package=viridis)  (0.5.1)

## Data

Publicly available data have been archived for download at the beginning of the analysis (doi = ZENODO). They require X GB of space within the subdirectory. Unfortunately, due to privacy considerations, the health outcome data (pediatric asthma ED visits and Medicare all-cause mortality in Texas) may not be publicly distributed. Instead, we have included synthetic data as a placeholder. There are **X main data sources:**

### Coal-fired power plant emissions

Emissions data are from the US Environmental Protection Agency's [Air Markets Program Data (AMPD)](https://ampd.epa.gov/ampd/). The AMPD provides access to current and historical monthly emissions data on electricity generating units (EGUs), collected as part of EPA's emissions trading programs. The downloaded AMPD are named `AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv`. This file contains monthly emissions totals of SO2 (and other chemical emissions), as well as characteristics of the EGUs from 1995-2017. The following columns are most relevant to our analysis:

- `Facility.ID`: unique six-digit facility identification number, also called an ORISPL, assigned by the Energy Information Administration
- `Unit.ID`: unique identifier for each unit at a facility
- `Year`: the calendar year during which activity occurred
- `Month`: the month in which activity occurred
- `Facility.Longitude`: the physical longitude of the facility
- `Facility.Latitutde`: the physical latitude of the facility
- `SO2..tons`: sulfur dioxide (SO2) emissions, in short tons
- `Has.SO2.Scrub`: denotes use of flue-gas desulfurization (FGD) emissions reduction technologies (i.e., a scrubber is in use) 

Additional variables are not relevant to this analysis, however, a comprehensive list of AMPD column definitions is included for reference (see the downloaded `AMPD_column_definitions.csv` file). 

### 2011 mean sulfate concentrations

Annual mean sulfate concentrations were obtained from the Randall Martin Atmospheric Composition Analysis Group's [North American Regional Estimates (V4.NA.03) dataset](https://sites.wustl.edu/acag/datasets/surface-pm2-5/#V4.NA.03). The data consist of raster annual mean SO4 data (micrograms per cubic meter; grid resolution = 0.01 x 0.01 degrees), as described in [van Donkelaar et al. (2019)](https://pubs.acs.org/doi/10.1021/acs.est.8b06392). The downloaded SO4 file is named `GWRwSPEC_SO4_NA_201101_201112.nc`. 

### 2010 US population density

US population density estimates were downloaded from the [USGS ScienceBase Catelog's](http://dx.doi.org/10.5066/F74J0C6M) block-level population density rasters. The 2010 dataset includes density rasters at 60-m resolution. The downloaded data file is named `pden2010_60m.tif`.

### Meteorological data

Meteorological data are from the NOAA Physical Science Laboratory [NCEP North American Regional Reanalysis (NARR)](https://psl.noaa.gov/data/gridded/data.narr.monolevel.html) database. Downloaded files include:
- `air.2m.mon.mean.nc`: monthly mean air temperature at 2 m
- `apcp.mon.mean.nc`: monthly average of daily accumulated total precipitation
- `rhum.2m.mon.mean.nc`: monthly mean relative humidity at 2 m
- `uwnd.10m.mon.mean.nc`: monthly mean U-wind at 10 m
- `vwnd.10m.mon.mean.nc`: monthly mean V-wind at 10 m
These files are read into R as raster layers.

### Data Preprocessing

Downloaded data are processed at the beginning of this analysis, using `make-facility-data.R` and `data-cleaning.R`. These two files produce a single R list object, saved as `central-usa-data.RDS`, which contains:

1. `so4`: raster of 2011 SO4 surface
2. `wind.big` and `wind.small`: rasters of wind vector elements, used to create operator matrix C
3. `em`: 2011 power plant facility emissions data
4. `X`: vector of annual SO2 emissions, matched to SO4 raster elements
5. `pop`: raster of 2010 US population density 
6. `temp`, `precip`, `rel.hum`: rasters of meteorological data

## Instructions

Before running any code, make sure the required R packages have been installed. Open and run the `main.R` file, found in the `./src/` folder.  Note that this script will take a long time to run sequentially. On a 2.3 GHz Intel Core i7 processor with 32 GB of memory, the script will take approximately 115 hours to run.  The script contains 9 sequential steps, which perform the following:

### Step One: 

- Loads required packages into R.
- Time: < 0.01 seconds.

### Step Two: 

- Creates `./data/` and `./output/` subdirectories, which will hold the underlying data and analysis output, respectively.
- Time: < 0.01 seconds.

### Step Three:

- Downloads the raw data sources used in the analysis. These data are publicly available, and have been archived (doi = 10.5281/zenodo.4072504) for reproducibility. These data require 1.7 GB of space.
- Time: ~2.5 minutes.
- Output size: 1.7 GB.

### Step Four: 

- Adds all functions from the script, `functions.R`, found in the `./src/` folder.
- Time: < 0.01 seconds.

### Step Five: 

- Loads the raw coal-fired power plant facilities data, cleans the data, and creates a data frame with relevant covariate values. After step three, the raw facility data are stored in the `./data/` folder as `AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv`. This section saves four output RDS files - `MonthlyUnitData.RDS`, `AnnualUnitData.RDS`, `MonthlyFacilityData.RDS`, and `AnnualFacilityData.RDS` - in the `./data/` folder. 
- Time: ~30 seconds.
- Output size: 13.3 MB.

### Step Six: 

- Cleans all data (including environmental covariates, the SO4 response variable, and facilities data), and stores them as a single raster. This raster object is saved as `./data/central-usa-data.RDS`. This raster contains all data needed for the remaining analysis. The relevant raw data sources can be found in the subfolder, `./data/`, created in Step Two (see Data section above for more details).
- Time: ~1 minute.
- Output size: 3.8 MB.

### Step Seven: 

- Generates posterior draws (via MCMC) from the 4 models considered in the manuscript. The samples are stored as RDS files in `./output/`. 
- Time: **CAUTION: THIS WILL TAKE A VERY LONG TIME**. To reproduce **all results in the manuscript takes ~111 hours**; to reproduce all results from the **supplementary materials takes an additional ~245 hours**. If possible, it is recommended that the individual steps 2-5 found in `./src/so4-mcmc.R` be completed in parallel. Note that steps 7-8 in `./src/so4-mcmc.R` reproduce results from the supplementary material; by default they are commented out.
- Output size: 17.5 MB (manuscript only); 69.6 MB (manuscript + supp. materials).

### Step Eight:  

- Summarizes the posterior draws with Figures, Tables, and results found in Section 4 of the manuscript. The figures will be saved as PNG files in `./output/`.  
- Time: **CAUTION: THIS MAY TAKE ~1.25 HOURS** (due to the creation of Figure 4c).
- Output size: 4.1 MB.

### Step Nine:

- Generates the plots found within `supp-materials.pdf`. These are saved as PNGs in `./output/`. 
- Time: ~10 seconds.
- Output size: 7.4 MB.

*Note: Step Nine is commented out in* `main.R`; *uncomment the source call to* `supp-plots.R` *and the last two numbered steps in* `so4-mcmc.R` *if you wish to reproduce the content in the Supplementary Material.*

## Output

Upon successful completion of `main.R`, the results are saved as PNGs in the `./output/ms-figures` (or `./output/supp-figures`) folder. Example format includes `./output/ms-figures/fig1a.png`, etc. These figures will look very similar, if not identical, to those found in the published manuscript. Differences can be explained by small changes induced by random draws from the posterior. However, all results should be qualitatively the same.
