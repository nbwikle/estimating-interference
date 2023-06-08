
### make-facility-data.R
### authors: Cory Zigler, Nathan Wikle
###
### Create data structures from emissions data for coal-fired power plant 
###   facilities in the United States.

### read in monthly emissions data (merged with other stuff)
dat <- fread(here::here("data", "AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv"))

###########################
###   UNIT-LEVEL DATA   ###
###########################

# Data are on "Electricity Generating Units (EGUs)."  There are
#   There are usually >=1 units co-located to make up a single "facility."
dat$uID = paste(dat$Facility.ID..ORISPL., dat$Unit.ID, sep = "_")
dat$FacID = dat$Facility.ID..ORISPL.
dat$Year = as.numeric(dat$Year)
dat$Month = as.numeric(dat$Month)
dat[, year_month := paste(Year, Month, sep="_")]
setkeyv(dat, c("uID", "Year", "Month"))
setorderv(dat, c("uID", "Year", "Month"))
dat <- dat[,-1]  # remove unique row identifiers
dat <- unique(dat) # remove duplicate entries
dat <- dat[!duplicated(dat[,17:20]),]

# create region variable
dat$Region <- region(dat$FIPS)

# Might want to restrict to only certain types of Units.
#   E.g., restrict to Coal-fired units only.
dat_unit <- base::subset(dat, Fuel.Type..Primary..x == 'Coal') ## unit level, monthly
setkeyv(dat_unit, c("uID", "Year", "Month"))

# Define SO2 control strategies at the Unit Level.
dat_unit <- SO2controltechnologies(dat_unit)
dat_unit$AnySO2control <- FALSE
dat_unit$AnySO2control <- (apply(dat_unit[, c("DryLimeFGD", "DrySorbInj",
  "DualAlk", "MagOx", "SodiumBased", "WetLimeFGD", "WetLime", "OtherSO2",
  "FluidizedBed"), with=FALSE], 1, sum) >=1)
# with(dat_unit, mean(AnySO2control))

# Define NOx control strategies at the Unit Level.
dat_unit <- NOxcontroltechnologies(dat_unit)
dat_unit$NumNOxControls <- 0
dat_unit$NumNOxControls <- apply(dat_unit[, c("SCR", "SNCR", "LowNOxBurner",
  "OverFire", "Ammonia", "CombustMod", "WaterInj", "OtherNOx"), with = FALSE],
  1, sum, na.rm = TRUE)
# with(dat_unit, table(NumNOxControls))

## Group SCR and SnCR as one category
dat_unit$S_n_CR <- apply(dat_unit[, c("SCR", "SNCR"), with = FALSE],
                         1, any, na.rm = TRUE)
# with(dat_unit, mean(S_n_CR))

# Define Operating Capacity as the maximum number of hours per month times the 
#   Max Hourly Heat Input Rate.
maxopp <- with(dat_unit, tapply(Operating.Time, year_month, max, na.rm = TRUE))
# Note: maxopp is the max total number of hours operating time for each of 12 months

# Max.Hourly.HI.Rate..MMBtu.hr. is the heat input capacity (MMBtu/hour).
for (i in 1:length(maxopp)) {
  dat_unit[
    year_month == names(maxopp)[i],
    Capacity := maxopp[i] * Max.Hourly.HI.Rate..MMBtu.hr.
  ] ## MMBtu/month
}

# # Percent Capacity is Heat Input / Capacity
# dat_unit$PctCapacity = dat_unit$Heat.Input..MMBtu / dat_unit$Capacity
# with(dat_unit, mean(PctCapacity>1.5, na.rm = TRUE)) ## 0.1% capacity >1.5
# dat_unit$PctCapacity[dat_unit$PctCapacity > 1.5] = NA ## > 150% capacity
# 
# # Define Sulfur Content Variables at the Unit Level
# dat_unit[, Sulfur.Content := as.numeric(Sulfur.Content)]
# dat_unit = CreateSulfurCategories(dat_unit, sulfurvar = dat_unit[, Sulfur.Content])

# Output monthly unit-level data set
saveRDS(dat_unit, file = here::here("data", "MonthlyUnitData.RDS"))

# Make an annual unit data set
yearly_dat_unit = dat_unit[, list(Sulfur.Content = mean(as.numeric(Sulfur.Content), na.rm = TRUE),
  Has.Scrubber = (sum(Has.SO2.Scrub, na.rm=TRUE) > 6), # Has a scrubber at least half the year
  Initial.Year.of.Operation = Initial.Year.of.Operation[1],
  S_n_CR = (sum(S_n_CR, na.rm = TRUE) >=6), # SnCR for at least half the year
  NumNOxControls = mean(NumNOxControls,na.rm = TRUE), # avg number of nox controls in the year
  Operating.Time = sum(Operating.Time, na.rm = TRUE),
  SO2..tons. = sum(SO2..tons., na.rm = TRUE),
  NOx..tons. = sum(NOx..tons., na.rm = TRUE),
  CO2..short.tons. = sum(CO2..short.tons., na.rm = TRUE),
  Gross.Load..MW.h. = sum(Gross.Load..MW.h., na.rm = TRUE),
  Heat.Input..MMBtu. = sum(Heat.Input..MMBtu., na.rm = TRUE),
  pctCapacity = mean(Heat.Input..MMBtu./Capacity, na.rm = TRUE),
  Is.Phase2 = Is.Phase2[1],
  Facility.Latitude.x = Facility.Latitude.x[1],
  Facility.Longitude.x = Facility.Longitude.x[1],
  FIPS = FIPS[1],
  Region = Region[1]),
  by = c("uID", "Year")]

# Save Annual Unit Level Data Set 
saveRDS(yearly_dat_unit, file = here::here("data", "AnnualUnitData.RDS"))

###############################
###   FACILITY-LEVEL DATA   ###
###############################

# Make facility-level data
dat_unit[, Sulfur.Content := as.numeric(Sulfur.Content)]
dat_facility = dat_unit[, list(nunits = length(unique(Unit.ID)),
  nunits_withsulfur = sum(!is.na(Sulfur.Content)),
  pctunits_withsulfur = sum(!is.na(Sulfur.Content))/length(unique(Unit.ID)),
  meanSulfur = sum(Sulfur.Content * Heat.Input..MMBtu., na.rm = TRUE)/sum(Heat.Input..MMBtu.[!is.na(Heat.Input..MMBtu.) & !is.na(Sulfur.Content)]),
  ScrubbedFacility = sum(Has.SO2.Scrub, na.rm = TRUE) == length(unique(Unit.ID)), ## Scrubbed if all units have an SO2 control
  initialYear = Initial.Year.of.Operation[1],
  pctS_n_CR = sum(S_n_CR, na.rm = TRUE)/length(unique(Unit.ID)),
  totNumNOxControls = sum(NumNOxControls, na.rm = TRUE),
  totOpTime = sum(Operating.Time, na.rm = TRUE),
  totSO2emissions = sum(SO2..tons., na.rm = TRUE),
  totNOxemissions = sum(NOx..tons., na.rm = TRUE),
  totCO2emissions = sum(CO2..short.tons., na.rm = TRUE),
  totLoad = sum(Gross.Load..MW.h., na.rm = TRUE),
  totHeatInput = sum(Heat.Input..MMBtu., na.rm = TRUE),
  pctCapacity = sum(Heat.Input..MMBtu.)/sum(Capacity),
  Phase2 = Is.Phase2[1],
  Fac.Latitude = Facility.Latitude.x[1],
  Fac.Longitude = Facility.Longitude.x[1],
  Fac.FIPS = FIPS[1],
  Fac.Region = Region[1]),
  by = c("FacID", "Year", "Month")]

with(dat_facility, mean(pctCapacity > 1.5, na.rm = TRUE)) ## 0.1% capacity > 1.5
dat_facility$pctCapacity[dat_facility$pctCapacity > 1.5] = NA ## > 150% capacity

setkeyv(dat_facility, c("FacID", "Year", "Month"))
setorderv(dat_facility, c("FacID", "Year", "Month"))

# save monthly facility-Level data 
saveRDS(dat_facility, file = here::here("data", "MonthlyFacilityData.RDS"))

# make an annual facility data set
yearly_dat_facility = dat_facility[, list(nunits = nunits[1],
  nunits_withsulfur = max(nunits_withsulfur),
  pctunits_withsulfur = max(pctunits_withsulfur),
  meanSulfur = max(meanSulfur),
  ScrubbedFacility = any(ScrubbedFacility[Month<=6]),
  initialYear = initialYear[1],
  pctS_n_CR = max(pctS_n_CR),
  totNumNOxControls = max(totNumNOxControls),
  totOpTime = mean(totOpTime),
  totSO2emissions = mean(totSO2emissions),
  totNOxemissions = mean(totNOxemissions),
  totCO2emissions = mean(totCO2emissions),
  totLoad = mean(totLoad),
  totHeatInput = mean(totHeatInput),
  pctCapacity = mean(pctCapacity, na.rm = TRUE),
  Phase2 = Phase2[1],
  Fac.Latitude = Fac.Latitude[1],
  Fac.Longitude = Fac.Longitude[1],
  Fac.FIPS = Fac.FIPS[1],
  Fac.Region = Fac.Region[1]),
  by = c("FacID", "Year")]

# save annual facility-level data 
saveRDS(yearly_dat_facility, file = here::here("data", "AnnualFacilityData.RDS"))
