### outcome-data.R
### Nathan Wikle

### load all relevant data

# asthma data
asthma.dt <- fread(here::here("data", "synth-ped-asthma-data.csv"))
asthma.dt[, 'GEOID'] <- sapply(asthma.dt[, 'GEOID'], as.character)

# census data
census.df <- readRDS(here::here("data", "Census_2016_TxZCTA.RDS"))

# climate data
climate.dt <- fread(here::here("data", "climate-data-2016-ZCTA.csv"))
climate.dt[, 'GEOID'] <- sapply(climate.dt[, 'GEOID'], as.character)

# smoking data
smoking.dt <- fread(here::here("data", "tx-zcta-2012smokers.csv"))
colnames(smoking.dt)[1] <- "GEOID"
smoking.dt[, 'GEOID'] <- sapply(smoking.dt[, 'GEOID'], as.character)

# facility data
fac.dt <- readRDS(here::here("data", "zcta-facility-data.RDS"))

### combine all data into data frame
outcome.df <- census.df %>%
  dplyr::left_join(climate.dt, by = "GEOID") %>%
  dplyr::left_join(smoking.dt, by = "GEOID") %>%
  dplyr::left_join(fac.dt, by = "GEOID") %>%
  dplyr::left_join(asthma.dt, by = "GEOID") %>%
  dplyr::select(., -asthma_rate)

### save data
saveRDS(outcome.df, here::here("data", "outcome-data.RDS"))



