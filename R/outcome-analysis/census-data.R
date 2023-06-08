### census-data.R
### Nathan Wikle

# WARNING:
#  The Census API can be finicky. You may get the following error messages: 
#
#   "The API message returned is Sorry, the system is currently undergoing 
#     maintenance or is busy.  Please try again later..""
#  OR
#   "You have supplied an invalid or inactive API key. To obtain a valid API key, 
#     visit https://api.census.gov/data/key_signup.html. To activate your key,
#     be sure to click the link provided to you in the email from the Census Bureau 
#     that contained your key."
#
#  In either case, repeated tries will eventually lead to successful completion
#   of this code.

###############################################################################
### I. load libraries, set up census API
###############################################################################

# libraries

library(maps)
library(tidycensus)
library(units)
library(tidyr)
library(sf)
library(stringr)

# census api (user must register for a key with Census Bureau)
user.key <- "Your US Census API key goes here!"
census_api_key(user.key)

###############################################################################
### II. query census data tables
###############################################################################

# ACS 5 data tables (2016)
ACS_2016 <- load_variables(2016, "acs5", cache = TRUE)

# decennial data tables (2010)
SF_2010 <- load_variables(2010, "sf1", cache = TRUE)

### 1. Demographic data (age, race, male/female) ###

pop_groups <- ACS_2016 %>% 
  filter(str_detect(ACS_2016$name, "B01001B|B01001H|B01001I")) %>% 
  mutate(Ecodes = paste0(name,"E"))

zip_pop_16 <- get_acs(
  geography = "zcta",
  variables = c(
    pop_groups$name,
    "B01001_001", "B01001_026", "B03002_012", "B03002_013", "B03002_014", "B01002_001"
  ),
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
)

pop_est <- zip_pop_16 %>%
  mutate(
    Total_16 = B01001_001E,
    White_ped_16 = rowSums(across(pop_groups$Ecodes[c(34:37, 49:52)])),  # 0-17 y.o., M and F, white (not hispanic or latino)
    Black_ped_16 = rowSums(across(pop_groups$Ecodes[c(3:6, 18:21)])),    # 0-17 y.o., M and F, black (not hispanic or latino)
    Latinx_ped_16 = rowSums(across(pop_groups$Ecodes[c(65:68, 80:83)])), # 0-17 y.o., M and F, hispanic or latino
    White_adu_16 = rowSums(across(pop_groups$Ecodes[c(38:47, 53:62)])),  # 18+ y.o., M and F, white (not hispanic or latino)
    Black_adu_16 = rowSums(across(pop_groups$Ecodes[c(7:16, 22:31)])),   # 18+ y.o., M and F, black (not hispanic or latino) 
    Latinx_adu_16 = rowSums(across(pop_groups$Ecodes[c(69:78, 84:93)])), # 18+ y.o., M and F, hispanic or latino
    L_total_16 = B03002_012E,                       # total, hispanic or latino (all races)
    WL_total_16 = B03002_013E,                      # total, hispanic or latino (white alone)
    BL_total_16 = B03002_014E,                      # total, hispanic or latino (black alone)
    WL_share_16 = B03002_013E / B03002_012E,        # percent white of total hispanic or latino pop.
    BL_share_16 = B03002_014E / B03002_012E,        # percent black of total hispanic or latino pop.
    pct_female = B01001_026E / B01001_001E,         # percent female
    Median_age = B01002_001E) %>%                    # median age
  dplyr::select(c(
    "GEOID", "Total_16", "White_ped_16", "Black_ped_16",
    "Latinx_ped_16", "White_adu_16", "Black_adu_16",
    "Latinx_adu_16", "L_total_16", "WL_total_16", "BL_total_16",
    "BL_share_16", "WL_share_16", "pct_female", "Median_age"
  ))

zip_age_16 <- get_acs(
  geography = "zcta",
  variables = c(
    "B01001_001", "B01001_003", "B01001_004", "B01001_005", "B01001_006", 
    "B01001_027", "B01001_028", "B01001_029", "B01001_030"
  ),
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
)

age_est <- zip_age_16 %>%
  mutate(
    Ped_16 = rowSums(across(c(
    "B01001_003E", "B01001_004E", "B01001_005E", "B01001_006E", 
    "B01001_027E", "B01001_028E", "B01001_029E", "B01001_030E"
  )))) %>%          # 0-17 y.o., M and F, all races            
  dplyr::select(c(
    "GEOID", "Ped_16"
  ))


### 2. Poverty status (B17) and Household Income (B19) ###

pov_groups <- ACS_2016 %>%
  filter(str_detect(ACS_2016$name, "B17001")) %>%
  mutate(Ecodes = paste0(name, "E"))

pov_est <- get_acs(
  geography = "zcta",
  variables = c(
    pov_groups$name[c(1:3, 17)]
  ),
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
) %>% mutate(
  pct_pov = B17001_002E / B17001_001E
)

# median household income in the past 12 months (in 2016 inflation-adjusted dollars) 
income_est <- get_acs(
    geography = "zcta",
    variables = "B19013_001",
    state = "Texas",
    year = 2016,
    survey = "acs5",
    output = "wide"
  ) %>% mutate(
    med_income = B19013_001E
  ) %>%
  dplyr::select(c(
    "GEOID", "med_income"
  ))

# # weird % poverty values (... consider removing from data set)
# pop_est$Total_16[c(1271, 1656, 1674, 1768, 1916)]
# zip_poverty[c(1271, 1656, 1674, 1768, 1916),]

### 3. Educational attainment (B15) ###

ed_groups <- ACS_2016 %>%
  filter(str_detect(ACS_2016$name, "B15003")) %>%
  mutate(Ecodes = paste0(name, "E"))

zip_ed <- get_acs(
  geography = "zcta",
  variables = ed_groups$name[c(1, 17:25)],
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
)

ed_est <- zip_ed %>%
  mutate(
    ed_pop = B15003_001E, # total pop. 25 years or older
    pct_hs_grad = rowSums(across(ed_groups$Ecodes[17:25])) / B15003_001E, # % hs graduate or higher
  ) %>% dplyr::select(c(
    "GEOID", "ed_pop", "pct_hs_grad"
  )) 
  

### 4. Migration rate (B07) ###

move_groups <- ACS_2016 %>%
  filter(str_detect(ACS_2016$name, "B07001")) %>%
  mutate(Ecodes = paste0(name, "E"))

zip_move <- get_acs(
  geography = "zcta",
  variables = move_groups$name[c(1, 17)],
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
)

move_est <- zip_move %>%
  mutate(
    move_pop = B07001_001E,
    move_rate = (B07001_001E - B07001_017E) / B07001_001E
  ) %>%
  dplyr::select(c(
    "GEOID", "move_pop", "move_rate"
  ))


### 5. Health Insurance (B27) ###

ins_groups <- ACS_2016 %>%
  filter(str_detect(ACS_2016$name, "B27001")) %>%
  mutate(Ecodes = paste0(name, "E"))

zip_ins <- get_acs(
  geography = "zcta",
  variables = ins_groups$name[c(
    1,
    4, 7, 10, 13, 16, 19, 22, 25, 28,
    32, 35, 38, 41, 44, 47, 50, 53, 56
  )],
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
)

ins_est <- zip_ins %>%
  mutate(
    ins_pop = B27001_001E, # total pop. eligible for health insurance
    pct_ins = rowSums(across(ins_groups$Ecodes[c( # sum up people with health insurance
      4, 7, 10, 13, 16, 19, 22, 25, 28, # males
      32, 35, 38, 41, 44, 47, 50, 53, 56 # females
    )])) / B27001_001E,
  ) %>%
  dplyr::select(c(
    "GEOID", "ins_pop", "pct_ins"
  ))


### 6. Occupied Housing (B25) ###

housing_groups <- ACS_2016 %>%
  filter(str_detect(ACS_2016$name, "B25008")) %>%
  mutate(Ecodes = paste0(name, "E"))

zip_housing <- get_acs(
  geography = "zcta",
  variables = housing_groups$name[c(1, 3)],
  state = "Texas",
  year = 2016,
  survey = "acs5",
  output = "wide"
)

housing_est <- zip_housing %>%
  mutate(
    renter_housing = B25008_003E / B25008_001E
  ) %>%
  dplyr::select(c(
    "GEOID", "renter_housing"
  ))


### 7. Proportion of urban residents (decennial: H00200) ###

urban_vars <- SF_2010 %>%
  filter(str_detect(SF_2010$name, "H00200"))

zip_urban <- get_decennial(
  geography = "zcta",
  variables = urban_vars$name[1:2],
  state = "Texas",
  year = 2010,
  sumfile = "sf1",
  output = "wide"
)

urban_data <- zip_urban %>%
  mutate(
    urban_prop = H002002 / H002001
  ) %>%
  dplyr::select(c(
    "GEOID", "urban_prop"
  ))

# remove leading FIPS code from GEOID 
new_geoid <- substring(urban_data$GEOID, 3) 

# there are 4 extra zip codes that are not in ACS 5 ...
#   all three are mostly in other states (OK, NM, NM, NM)
setdiff(new_geoid, pop_est$GEOID)

# remove extra zips
matching_geoid <- new_geoid[-match(setdiff(new_geoid, pop_est$GEOID), new_geoid)]
urban_df <- urban_data[-match(setdiff(new_geoid, pop_est$GEOID), new_geoid), ]
urban_df$GEOID <- matching_geoid

###############################################################################
### III. Combine data into single table
###############################################################################

# combine via GEOID
zip_data_full <- pop_est %>%
  left_join(age_est, by = "GEOID") %>%
  left_join(pov_est, by = "GEOID") %>%
  left_join(income_est, by = "GEOID") %>%
  left_join(ed_est, by = "GEOID") %>%
  left_join(move_est, by = "GEOID") %>%
  left_join(ins_est, by = "GEOID") %>%
  left_join(housing_est, by = "GEOID") %>% 
  left_join(urban_df, by = "GEOID")

# keep relevant covariates
zip_data <- zip_data_full %>%
  mutate(
    log_pop = log(Total_16),
    female = pct_female,
    ped_prop = Ped_16 / Total_16,
    median_age = Median_age,
    white_prop = (White_ped_16 + White_adu_16) / Total_16,
    black_prop = (Black_ped_16 + Black_adu_16) / Total_16,
    hisp_prop = (Latinx_ped_16 + Latinx_adu_16) / Total_16, 
    hs_prop = pct_hs_grad, 
    pov_prop = pct_pov,
    log_income = log(med_income),
    move_rate = move_rate,
    ins_prop = pct_ins,
    renter_housing = renter_housing,
    urban_prop = urban_prop
  ) %>% 
  dplyr::select(c(
    "GEOID", "log_pop", "female", "ped_prop", "median_age", "white_prop", "black_prop",
    "hisp_prop", "hs_prop", "pov_prop", "log_income", "move_rate", "ins_prop",
    "renter_housing", "urban_prop"
  ))
  

###############################################################################
### IV. Add shape file, remove missing data
###############################################################################

zip_geoms <- get_acs(
  geography = "zcta",
  variables = "B01001_001",
  year = 2016,
  state = "Texas",
  geometry = TRUE,
  cb = FALSE
)

zip_geoms <- zip_geoms %>%
  dplyr::select(GEOID, NAME, variable, estimate) %>%
  tidyr::spread(variable, estimate) %>%  
  dplyr::select(c(
    "GEOID", "geometry"
  ))

# merge demographic and geometry data into a single object
census_df <- zip_geoms %>%
  left_join(zip_data, by = "GEOID")

# calculate area per zipcode
zip_area <- st_area(census_df)
units(zip_area) <- make_units(km^2)
zip_area <- as.vector(zip_area)

# add log(pop / km^2)
census_df <- census_df %>%
  mutate(
    log_density = log_pop - log(zip_area)
  ) %>%
  dplyr::select(c(
    "GEOID", "log_pop", "female", "ped_prop", "median_age", "white_prop",
    "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income", "move_rate", 
    "ins_prop", "renter_housing", "urban_prop",  "log_density", "geometry"
  ))


###############################################################################
### V. Save data
###############################################################################

saveRDS(census_df, file = here::here("data", "Census_2016_TxZCTA.RDS"))













