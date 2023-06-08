### analysis-setup.R
### Nathan Wikle


###############################################################################
### 1. Packages/data/functions
###############################################################################

# # source functions
# source(here::here("R", "outcome-analysis", "outcome-functions.R"))

# read in Texas SO4 and emissions data
so4.data <- readRDS(here::here("data", "texas-2016-data.RDS"))

# outcome data
outcome.data <- readRDS(here::here("data", "outcome-data.RDS"))
# remove missing health data
out.df <- na.omit(outcome.data)

# facility data
fac.df <- readRDS(here::here("data", "facility-data-2016.RDS"))



###############################################################################
### 2. Determine relative influence of each power plant
###############################################################################

# zip geometries
zip.geoms <- out.df[, "geometry"]

# so4 raster (convert NA to 0 so that raster projection is correct)
so4.raster <- so4.data$so4
values(so4.raster)[is.na(values(so4.raster))] <- 0

# caution: removes NA elements
so4.stars <- st_as_stars(so4.raster)
so4.sf <- st_as_sf(so4.stars)

# match crs
so4.sf <- so4.sf %>% st_transform(st_crs(zip.geoms))

# determine proportion of each raster element per zip
n.zips <- nrow(zip.geoms)
n.ras.elems <- length(values(so4.raster)) 
props.mat <- Matrix(data = 0, nrow = n.ras.elems, ncol = n.zips, sparse = TRUE)

# turn on spherical geometries
sf_use_s2(TRUE)

for (zip in 1:n.zips) {
  # intersection with raster elements
  int.sf <- st_intersection(so4.sf, zip.geoms[zip, ])
  # determine area
  int.sf <- int.sf %>%
    dplyr::mutate(area = st_area(int.sf))
  # determine relevant raster indices
  ras.int <- st_intersects(so4.sf, zip.geoms[zip, ], sparse = FALSE)
  ras.inds <- ras.int[, 1]
  # calculate proportion of each raster element that makes up the zip
  props.mat[ras.inds, zip] <- drop_units(int.sf$area / sum(int.sf$area))
}

# save props.mat 
saveRDS(props.mat, file = here::here("data", "projection-mat.RDS"))

### Create emissions matrix

n.fac <- nrow(so4.data$em$em.usa$fac)
n.ras.elems <- length(values(so4.data$so4)) 
X.mat <- Matrix(0, nrow = n.ras.elems, ncol = n.fac)

for (k in 1:n.fac){
  em.k <- so4.data$em$em.usa$fac
  em.k[, 12] <- rep(0, n.fac)
  em.k[k, 12] <- 1000
  X.mat[, k] <- createX(so4.data$so4, em.k)
}

# remove duplicate facilities (i.e., facilities are in the same place)
#   -> facilities 6181 (JT Deely) & 7097 (JK Spruce)  (60/71)
#   -> facilities 6648 (Sandow No. 4) & 52071 (Sandow Station)  (71/78)
em.mat <- X.mat[,-c(71,78)]

# save results
saveRDS(em.mat, file = "./data/emissions-mat.RDS")


### vectors containing key-associated and scrubber info

# vector with key-associated facilities
fac.trim <- fac.df
all.ids <- unique(fac.trim$FacID)
key.assoc <- rep(NA, nrow(out.df))

for (k in 1:nrow(out.df)){
  key.assoc[k] <- which(all.ids == out.df$fac_id[k])
}

# save data
saveRDS(key.assoc, here::here("data", "key-assoc-vec.RDS"))

# vector with scrubber status of each facility
scrubber.vec <- as.numeric(fac.trim$ScrubbedFacility)

# save data
saveRDS(scrubber.vec, here::here("data", "scrubber-vec.RDS"))

### advection-diffusion matrices

# create D matrix
dims <- dim(so4.data$so4)
inds.d <- FVMInds(dims, "insulated")
D <- diffFVM(dims, inds.d, "insulated")

# create C matrix
inds.w <- FVMInds(dims, "periodic")
C <- createC(
  big = so4.data$wind.big$wind.2016,
  small = so4.data$wind.small$wind.2016,
  dims, inds = inds.w
)

# combine in a list
advdiff.mats <- list(D = D, C = C)

# save data
saveRDS(advdiff.mats, here::here("data", "advec-diff-mats.RDS"))


