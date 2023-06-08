### so4-data.R
### Nathan Wikle

### Creates a list containing:
###     i. raster of annual average SO4 concentrations from 2016,
###    ii. raster of annual average 10m wind velocities from 2016,
###   iii. power plant emissions data from 2016 (USA and Mexico).

#################################################################
### 1. SO4 data
#################################################################

# 2016 SO4 data
file.in <- here::here("data", "GWRwSPEC_SO4_NA_201601_201612.nc")
so4.16 <- raster(file.in)

# change to different resolution
so4.16.newres <- aggregate(so4.16, 16)

# trim to region around Texas
tx.extent <- c(-110, -88, 25, 40)
so4.tx <- raster()
extent(so4.tx) <- extent(tx.extent)
res(so4.tx) <- res(so4.16.newres)
so4.tx <- crop(so4.16.newres, so4.tx)

# bigger extent (for wind structure)
tx.big.extent <- c(-110.5, -87.5, 24.5, 40.5)
so4.tx.big <- raster()
extent(so4.tx.big) <- extent(tx.big.extent)
res(so4.tx.big) <- res(so4.16.newres)
so4.tx.big <- crop(so4.16.newres, so4.tx.big)


#################################################################
### 2. Create emissions vectors
#################################################################

### i. USA

# store monthly data in a list
month <- list()
month$fac <- readRDS(here::here("data", "MonthlyFacilityData.RDS"))
month$unit <- readRDS(here::here("data", "MonthlyUnitData.RDS"))

# clean monthly unit data
month$unit <- cleanData(month$unit)

# store annual data in a list
annual <- list()
annual$fac <- readRDS(here::here("data", "AnnualFacilityData.RDS"))
annual$unit <- readRDS(here::here("data", "AnnualUnitData.RDS"))

# rename columns
colnames(month$unit) <- unitNames()
colnames(annual$unit) <- unitNames(TRUE)

# trim emissions data by year
month.2016 <- trimYear(month, 2016)
annual.2016 <- trimYear(annual, 2016)

# trim emissions data by spatial extent
em.2016 <- trimData(annual.2016, tx.extent)

# save yearly emissions data in a list
emissions <- list(em.usa = em.2016) 

# convert emissions data to X vector (design vector used in model)
X.2016 <- createSimpleX(so4.tx, emissions$em.usa$fac)

# save emissions vectors in list
X.list <- list(X.usa = X.2016)

### ii. Mexico

# read mexico emissions data
#  -> source: The 2016v1 emissions modeling platform from 
#       the National Emissions Inventory Collaborative
#     description: (https://www.epa.gov/air-emissions-modeling/2014-2016-version-7-air-emissions-modeling-platforms)
#       
#  -> data: https://gaftp.epa.gov/Air/emismod/2016/v1/2016emissions/2016fh_inventory_oth_27sep2019.zip)
mx.data <- fread(here::here("data", "Mexico_2016_point_interpolated_02mar2018_v0.csv"))

# remove columns with no data
mx.data <- mx.data[, c(1:2, 4:7, 12:14, 16:20, 22:25, 37)]

# keep SO4 pollution (ie, poll == SO2)
mx.data <- mx.data %>%
  dplyr::filter(poll == "SO2") 

# restrict to region of interest
mx.data <- mx.data %>%
  dplyr::filter(
    (longitude > -110) &
    (longitude < -88) &
    (latitude > 25) &
    (latitude < 40)
  )

# combine units into facility-level data
mx.fac <- mx.data %>%
  group_by(facility_name) %>%
  summarise(
    country_cd = dplyr::first(country_cd),
    region_cd = dplyr::first(region_cd),
    facility_id = dplyr::first(facility_id), 
    n.units = length(unit_id),
    scc = mode(scc),
    total.SO2 = sum(ann_value),
    longitude = dplyr::first(longitude),
    latitude = dplyr::first(latitude), 
    projection_factor = mode(projection_factor)
  )

emissions$em.mx <- mx.fac

# create X vector
X.list$X.mexico <- createMexicoX(so4.tx, mx.fac)


#################################################################
### 3. Add wind velocities
#################################################################

# grab 2011 wind data
wind.big <- list()
wind.norm <- list()

uwind.raster <- yearRaster(2016, met = "u-wind")
vwind.raster <- yearRaster(2016, met = "v-wind")

# trim to raster size  
uwind.new <- trimToSO4(uwind.raster, so4.tx)
vwind.new <- trimToSO4(vwind.raster, so4.tx)
  
uwind.big <- trimToSO4(uwind.raster, so4.tx.big)
vwind.big <- trimToSO4(vwind.raster, so4.tx.big)
  
# save as "wind.big" and "wind.norm"
big <- stack(uwind.big, vwind.big)
names(big) <- c("uwind", "vwind")
  
normal <- stack(uwind.new, vwind.new)
names(normal) <- c("uwind", "vwind")
  
wind.big[[1]] <- big
wind.norm[[1]] <- normal

names(wind.big) <- c("wind.2016") 
names(wind.norm) <- c("wind.2016") 


#################################################################
### 4. save data in list
#################################################################

# save results for future use
texas.data <- list(so4 = so4.tx,
                         wind.big = wind.big,
                         wind.small = wind.norm,
                         X = X.list,
                         em = emissions)

saveRDS(texas.data, here::here("data", "texas-2016-data.RDS"))

















