### bc-data.R
### Nathan Wikle

### Required packages ###
# library(raster)
# library(maps)
# library(ncdf4)
# library(data.table)
# library(stringr)
# library(rgdal)
# library(Matrix)

###############################################################################
### 1. Read in BC data, create raster
###############################################################################

### Read in SO4 data; files from Randall Martin's group available here:
###    https://sites.wustl.edu/acag/datasets/surface-pm2-5/

# 2016 BC data
file.in <- here::here("data", "GWRwSPEC_BC_NA_201601_201612.nc")
black.carbon <- raster(file.in)
so4.data <- readRDS(here::here("data", "texas-2016-data.RDS")) 
bc <- crop(black.carbon, so4.data$so4)
bc <- aggregate(bc, 16)

###############################################################################
### 2. Aggregate to ZIP level
###############################################################################

# projection from raster to TX zips
proj.mat <- readRDS(here::here("data", "projection-mat.RDS")) 

# black carbon raster values
bc.vec <- values(bc)

# remove any influence from NA raster elements
proj.mat[is.na(bc.vec), ] <- 0
for (k in 1:ncol(proj.mat)) {
   column_k <- proj.mat[,k]
   proj.mat[,k] <- column_k / sum(column_k)
}

# calculate ZCTA black carbon
bc.zips <- Matrix::crossprod(proj.mat, bc.vec)

# save bc.zips
saveRDS(bc.zips, file = here::here("data", "bc-2016-zcta.RDS"))


