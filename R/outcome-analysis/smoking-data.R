### smoking-data.R
### Nathan Wikle

### Process the smoking data (source: https://doi.org/10.1186/1478-7954-12-5)

###############################################################################
### 1. Read data
###############################################################################

# smoking data
smoking.tb <- fread(here::here("data", "smokedatwithfips_1996-2012.csv"), drop = 1)

###############################################################################
### 2. Restrict to Texas counties, year == 2012
###############################################################################
smoking.tx <- smoking.tb[State == "TX" & year == 2012]

###############################################################################
### 3. Convert county values to zipcodes
###############################################################################

# county-to-zip crosswalk file created by the Missouri Census Data Center
#   (https://mcdc.missouri.edu/applications/geocorr2014.html). It uses 
#   population counts from the 2010 census.

# read crosswalk file
tx.crosswalk <- fread(here::here("data", "tx-zip-to-county.csv"))

# remove first row (column description)
col.description <- tx.crosswalk[1,]
tx.crosswalk <- tx.crosswalk[-1, ]
tx.crosswalk <- tx.crosswalk[order(zcta5, county)]

# create allocation matrix
zcta.alloc <- sparseMatrix(
  i = as.integer(as.factor(tx.crosswalk$zcta5)),
  j = as.integer(as.factor(tx.crosswalk$county)),
  dimnames = list(
    as.character(levels(as.factor(tx.crosswalk$zcta5))),
    as.character(levels(as.factor(tx.crosswalk$county)))
  ),
  x = as.numeric(tx.crosswalk$afact)
)

# calculate smokerate in each zcta
#   (weighted average based on population of zcta from each county)
zcta.smokerate <- matrix(zcta.alloc %*% smoking.tx[order(FIPS), smokerate])
colnames(zcta.smokerate) <- "smokerate"
rownames(zcta.smokerate) <- rownames(zcta.alloc)

###############################################################################
### 4. Save data
###############################################################################

# save data
write.table(
  zcta.smokerate,
  file = here::here("data", "tx-zcta-2012smokers.csv")
)

