### climate-data.R
### Nathan Wikle

### Process the climate data. We want average daily max and min temp., 
###   total precipitation, and average relative humidity. Temp and precip 
###   data are from ORNL DAAC's Daymet V4 model 
###   (https://doi.org/10.3334/ORNLDAAC/1852).


###############################################################################
### 1. Daymet climate data (tmin, tmax, pcrp, vp)
###############################################################################

### A. load climate data

# load Daymet climate data
tmin_stars <- read_ncdf(here::here("data", "daymet_v4_tmin_annavg_na_2016.nc"), var = "tmin")
tmax_stars <- read_ncdf(here::here("data", "daymet_v4_tmax_annavg_na_2016.nc"), var = "tmax")
prcp_stars <- read_ncdf(here::here("data", "daymet_v4_prcp_annttl_na_2016.nc"), var = "prcp")
vp_stars <- read_ncdf(here::here("data", "daymet_v4_vp_annavg_na_2016.nc"), var = "vp")

# combine into single stars object
climate_stars <- c(tmin_stars, tmax_stars, prcp_stars, vp_stars)

### B. load census data

# load data
census_df <- readRDS(here::here("data", "Census_2016_TxZCTA.RDS"))
# match crs with raster data
census_matching <- st_transform(census_df, st_crs(climate_stars))

### C. aggregate raster values to zips

# calculate mean value for each zipcode
zip_means <- matrix(NA_real_, nrow = nrow(census_matching), ncol = 5)
colnames(zip_means) <- c("tmin", "tmax", "prcp", "vp", "rel_humid")
system.time(
  for (k in 1:nrow(census_matching)) {
    sf_obj <- st_transform(census_matching[k, "geometry"], st_crs(climate_stars))
    zip_means[k, 1] <- mean(climate_stars[sf_obj]$tmin, na.rm = TRUE)
    zip_means[k, 2] <- mean(climate_stars[sf_obj]$tmax, na.rm = TRUE)
    zip_means[k, 3] <- mean(climate_stars[sf_obj]$prcp, na.rm = TRUE)
    zip_means[k, 4] <- mean(climate_stars[sf_obj]$vp, na.rm = TRUE)
  }
)

###############################################################################
### 2. North American Regional Reanalysis (NARR) data (relative humidity)
###############################################################################

### A. Input data as a raster object

# year of interest
year <- 2016
# relevant raster bands
bands <- ((year - 1979) * 12) + 1:12
# input rasters into a list
raster_month <- list()
  
for (i in 1:length(bands)) {
  # grab raster for each month in the year of interest
  raster_month[[i]] <- raster(
    here::here("data", "rhum.2m.mon.mean.nc"), band = bands[i]
  )
}

# combine into a raster stack
rh_stack <- stack(raster_month)

# average over the 12 months
rh_avg <- raster::calc(rh_stack, fun = mean)

### B. Convert CRS and resolution to match Daymet data

# tmin raster
tmin_raster <- raster(
  here::here("data", "daymet_v4_tmin_annavg_na_2016.nc"), var = "tmin"
)

# project rh onto CRS and resolution of tmin
rh_avg2 <- projectRaster(
  rh_avg,
  crs = crs(tmin_raster),
  res = res(tmin_raster)
)

# crop to same extent of tmin
rh_avg_final <- crop(rh_avg2, tmin_raster)

### C. aggregate raster values to zips

# change raster to stars object
rh_stars <- st_as_stars(rh_avg_final)
names(rh_stars) <- "rh"

# calculate average rh for each zip
for (k in 1:nrow(census_matching)) {
  sf_obj <- st_transform(census_matching[k, "geometry"], st_crs(rh_stars))
  zip_means[k, 5] <- mean(rh_stars[sf_obj]$rh, na.rm = TRUE)
}


###############################################################################
### 3. Save data
###############################################################################

### A. add GEOID to climate data

# GEOIDs
GEOID <- census_matching[, "GEOID"] %>% st_drop_geometry()
# combined data
climate_data <- data.frame(GEOID, zip_means)

### B. save data

# save data as csv
write.csv(
  x = climate_data,
  file = here::here("data", "climate-data-2016-ZCTA.csv"), 
  row.names = FALSE
)


