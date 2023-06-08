### facility-data.R
### Nathan Wikle

###############################################################################
### 1. Load data 
###############################################################################

# load census data
census_df <- readRDS(here::here("data", "Census_2016_TxZCTA.RDS"))
# load sulfate data
so4_data <- readRDS(here::here("data", "texas-2016-data.RDS")) 
# get SO2 power plant facility data
fac_data <- so4_data$em$em.usa$fac


###############################################################################
### 2. Key-Associated Power Plants
###############################################################################

### Determine key-associated facilities

# facility coordinates
fac_locs <- fac_data[, 20:19]
names(fac_locs) <- c("long", "lat")

# convert long/lat coords to sf "points"
fac_coords <- fac_locs %>%
  st_as_sf(
    coords = c("long", "lat"),
    crs = st_crs(census_df)
  )

# calculate zcta centroids
zcta_centroids <- census_df %>% mutate(zip_centroids = st_centroid(geometry))

# key-associated power plants
key_plants <- apply(
  set_units(st_distance(zcta_centroids$zip_centroids[1:1935], fac_coords), km),
  1, which.min
)

# add to zcta_centroids df
zcta_centroids$key_plant <- as.factor(key_plants)

### combine SA and Austin facilities (these facilities are in the same location)
#   -> facilities 6181 (JT Deely) & 7097 (JK Spruce)  (60/71)
#   -> facilities 6648 (Sandow No. 4) & 52071 (Sandow Station)  (71/78)

# "new" version of facility data
fac_data.cleaned <- fac_data

# SA facilities
fac_data.cleaned[60, c(3:4, 10:16)] <- fac_data[60, c(3:4, 10:16)] +
  fac_data[71, c(3:4, 10:16)]
fac_data.cleaned[60, 5] <- 0.5
fac_data.cleaned[60, 6] <- (fac_data[60, 6] + fac_data[71, 6]) / 2
fac_data.cleaned[60, 17] <- (fac_data[60, 17] + fac_data[71, 17]) / 2

# change facility 60 to "scrubbed" status
fac_data.cleaned[60, 7] <- TRUE

# Austin facilities
fac_data.cleaned[67, c(3:4, 10:16)] <- fac_data[67, c(3:4, 10:16)] +
  fac_data[78, c(3:4, 10:16)]
fac_data.cleaned[67, 5] <- 2/3
fac_data.cleaned[67, 6] <- (fac_data[67, 6] + fac_data[78, 6]) / 2
fac_data.cleaned[67, 17] <- (fac_data[67, 17] + fac_data[78, 17]) / 2

# remove 71/78
fac_data.final <- fac_data.cleaned[-c(71,78),]

# save data
saveRDS(fac_data.final, file = here::here("data", "facility-data-2016.RDS"))

###############################################################################
### 3. Create power plant covariates
###############################################################################

# convert key_plants
key_plants.cleaned <- key_plants
key_plants.cleaned[which(key_plants == 71)] <- 60
key_plants.cleaned[which(key_plants == 78)] <- 67

# relevant facility covariates
fac_covars <- fac_data.cleaned %>% mutate(
    fac_id = FacID,
    log_heat = log(totHeatInput),
    log_optime = log(totOpTime),
    scrubbed = 1 * ScrubbedFacility,
    tot_nox = totNumNOxControls,
    pct_capacity = pctCapacity,
    pct_SnCR = pctS_n_CR, 
    phase2 = Phase2
  ) %>%
  dplyr::select(c(
    "fac_id", "scrubbed", "tot_nox", "log_heat", "log_optime",
    "pct_capacity", "pct_SnCR", "phase2"
  )) 

# create facility covariate using key-associated plant information
fac_df <- cbind(census_df %>% pull(GEOID), fac_covars[key_plants.cleaned, ])
colnames(fac_df)[1] <- 'GEOID'

# add distance to key-associated plants
all_distances <- set_units(st_distance(zcta_centroids$zip_centroids[1:1935], fac_coords), km)
fac_df$dist_to_key <- all_distances[cbind(1:1935, key_plants.cleaned[1:1935])]

# save covariate data
saveRDS(fac_df, here::here("data", "zcta-facility-data.RDS"))

