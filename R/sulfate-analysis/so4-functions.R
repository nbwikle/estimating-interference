### so4-functions.R
### Nathan Wikle
###
### A collection of function definitions used throughout the sulfate analysis.

################################################################################
### 1. Functions for data cleaning.
################################################################################

### A. Cleaning emissions data

CreateSulfurCategories = function(dat, sulfurvar){
  # Author: Cory Zigler
  # Create sulfur categories based on EPA classifications.
  # Data comes from EIA_767: 
  #   http://www.eia.gov/survey/form/eia_767/instructions_form.pdf:
  # Sulfur Content and Ash Content, columns (d) and (e), (i) and (j), 
  #   report content to nearest 0.01 percent for sulfur and the nearest 0.1 % 
  #   for ash. 
  # From John Bachman on 8/28/2014: 
  # --- Looking up some numbers quickly on the net, I found this: "Subbituminous 
  # --- Wyoming coal is only 0.35 percent sulfur by weight, while Kentucky coal 
  # --- is 1.59 percent sulfur.” I certainly recall Indiana, West Va coals in 
  # --- the 2% range, and that western coal number sounds right.
  # From Barrett on 9/3/2014:
  # --- A Feb 1993 Report from EIA (Appendix B Table C1) creates categories: 
  # --- <=0.60 (Low Sulfur), 0.61 - 1.67 (Medium Sulfur), >=1.68 (High Sulfur).

  # EPA's sulfur categories
  SulfurCatCuts = c(0, 0.60, 1.67, max(sulfurvar, na.rm = TRUE))
  dat$SulfurCat = as.numeric(cut(sulfurvar, breaks = SulfurCatCuts, 
                                 right = TRUE, include.lowest = TRUE))
  
  print("Variable Created: SulfurCat, equal to 1,2,3 for Low, Medium, High Sulfur Coal, respectively.")
  
  dat[, LowScoal := SulfurCat==1]
  print("Variable Created: LowScoal, indicator that SulfurCat==1")
  dat[, MedScoal := SulfurCat==2]
  print("Variable Created: MedScoal, indicator that SulfurCat==2")
  dat[, HighScoal := SulfurCat==3]
  print("Variable Created: LowScoal, indicator that SulfurCat==3")
  
  return(dat)
}

NOxcontroltechnologies = function(dat){
  # Author: Cory Zigler.
  # Takes unit-level data set and creates several indicators denoting the
  #   presence/absence of various Nox control strategies.
  # Input:
  #   dat: data frame with NOx variables
  # Output:
  #   Data frame with new indicator variables.
  
  # create various indicator variables
  dat$SCR = FALSE
  dat$SCR[dat$NOx.Scrub1 == "Selective Catalytic Reduction" | 
            dat$NOx.Scrub2 == "Selective Catalytic Reduction" |  
            dat$NOx.Scrub3 == "Selective Catalytic Reduction" |
            dat$NOx.Scrub4 == "Selective Catalytic Reduction"] = TRUE
  print("Variable Created: SCR")
  
  dat$SNCR = FALSE
  dat$SNCR[dat$NOx.Scrub1 == "Selective Non-catalytic Reduction" | 
             dat$NOx.Scrub2 == "Selective Non-catalytic Reduction" |  
             dat$NOx.Scrub3 == "Selective Non-catalytic Reduction" |
             dat$NOx.Scrub4 == "Selective Non-catalytic Reduction"] = TRUE
  print("Variable Created: SNCR")
  
  dat$LowNOxBurner = FALSE
  dat$LowNOxBurner[substr(dat$NOx.Scrub1, 1, 7) == "Low NOx" | 
                     substr(dat$NOx.Scrub2, 1, 7) == "Low NOx" | 
                     substr(dat$NOx.Scrub3, 1, 7) == "Low NOx" | 
                     substr(dat$NOx.Scrub4, 1, 7) == "Low NOx" ] = TRUE
  print("Variable Created: LowNOxBurner")
  
  dat$OverFire = FALSE
  dat$OverFire[dat$NOx.Scrub1 == "Overfire Air" | 
                 dat$NOx.Scrub2 == "Overfire Air" | 
                 dat$NOx.Scrub3 == "Overfire Air" | 
                 dat$NOx.Scrub4 == "Overfire Air" ] = TRUE
  print("Variable Created: OverFire")
  
  dat$Ammonia = FALSE
  dat$Ammonia[dat$NOx.Scrub1 == "Ammonia Injection" | 
                dat$NOx.Scrub2 == "Ammonia Injection" | 
                dat$NOx.Scrub3 == "Ammonia Injection" | 
                dat$NOx.Scrub4 == "Ammonia Injection" ] = TRUE
  print("Variable Created: Ammonia")
  
  dat$CombustMod = FALSE
  dat$CombustMod[dat$NOx.Scrub1 == "Combustion Modification/Fuel Reburning" | 
                   dat$NOx.Scrub2 == "Combustion Modification/Fuel Reburning" | 
                   dat$NOx.Scrub3 == "Combustion Modification/Fuel Reburning" | 
                   dat$NOx.Scrub4 == "Combustion Modification/Fuel Reburning" ] = TRUE
  print("Variable Created: CombustMod")
  
  dat$WaterInj = FALSE
  dat$WaterInj[dat$NOx.Scrub1 == "Water Injection" | 
                 dat$NOx.Scrub2 == "Water Injection" | 
                 dat$NOx.Scrub3 == "Water Injection" | 
                 dat$NOx.Scrub4 == "Water Injection" ] = TRUE
  print("Variable Created: WaterInj")
  
  dat$OtherNOx = FALSE
  dat$OtherNOx[dat$NOx.Scrub1 == "Other" | 
                 dat$NOx.Scrub2 == "Other" | 
                 dat$NOx.Scrub3 == "Other" | 
                 dat$NOx.Scrub4 == "Other" ] = TRUE
  print("Variable Created: OtherNOx")
  
  return(dat)
}

SO2controltechnologies = function(dat){
  # Author: Cory Zigler.
  # Takes unit-level data set and creates several indicators denoting the
  #   presence/absence of various SO2 control strategies.
  # Input:
  #   dat: data frame with SO2 control variables.
  # Output:
  #   Add indicators denoting the presence/absence of SO2 control strategies.
  
  dat$DryLimeFGD = FALSE
  dat$DryLimeFGD[(dat$SO2.Scrub1 == "Dry Lime FGD") | 
                   (dat$SO2.Scrub2 %in% "Dry Lime FGD") |  
                   (dat$SO2.Scrub3 %in% "Dry Lime FGD") | 
                   (dat$SO2.Scrub4 %in% "Dry Lime FGD")] = TRUE
  print("Variable Created: DryLimeFGD")
  #with(dat, mean(DryLimeFGD))
  
  dat$DrySorbInj = FALSE
  dat$DrySorbInj[(dat$SO2.Scrub1 == "Dry Sorbent Injection") | 
                   (dat$SO2.Scrub2 %in% "Dry Sorbent Injection") |  
                   (dat$SO2.Scrub3 %in% "Dry Sorbent Injection") | 
                   (dat$SO2.Scrub4 %in% "Dry Sorbent Injection")] = TRUE
  print("Variable Created: DrySorbInj")
  #with(dat, mean(DrySorbInj))
  
  dat$DualAlk = FALSE
  dat$DualAlk[(dat$SO2.Scrub1 == "Dual Alkali") | 
                (dat$SO2.Scrub2 %in% "Dual Alkali") |  
                (dat$SO2.Scrub3 %in% "Dual Alkali") | 
                (dat$SO2.Scrub4 %in% "Dual Alkali")] = TRUE
  print("Variable Created: DualAlk")
  #with(dat, mean(DualAlk))
  
  dat$FluidizedBed = FALSE
  dat$FluidizedBed[(dat$SO2.Scrub1 == "Fluidized Bed Limestone Injection") | 
                     (dat$SO2.Scrub2 %in% "Fluidized Bed Limestone Injection") |  
                     (dat$SO2.Scrub3 %in% "Fluidized Bed Limestone Injection") | 
                     (dat$SO2.Scrub4 %in% "Fluidized Bed Limestone Injection")] = TRUE
  print("Variable Created: FluidizedBed")
  #with(dat, mean(FluidizedBed))
  
  
  dat$MagOx = FALSE
  dat$MagOx[(dat$SO2.Scrub1 == "Magnesium Oxide") | 
              (dat$SO2.Scrub2 %in% "Magnesium Oxide") |  
              (dat$SO2.Scrub3 %in% "Magnesium Oxide") | 
              (dat$SO2.Scrub4 %in% "Magnesium Oxide")] = TRUE
  print("Variable Created: MagOx")
  #with(dat, mean(MagOx))
  
  dat$SodiumBased = FALSE
  dat$SodiumBased[(dat$SO2.Scrub1 == "Sodium Based") | 
                    (dat$SO2.Scrub2 %in% "Sodium Based") |  
                    (dat$SO2.Scrub3 %in% "Sodium Based") | 
                    (dat$SO2.Scrub4 %in% "Sodium Based")] = TRUE
  print("Variable Created: SodiumBased")
  #with(dat, mean(SodiumBased))
  
  dat$WetLimeFGD = FALSE
  dat$WetLimeFGD[(dat$SO2.Scrub1 == "Wet Lime FGD") | 
                   (dat$SO2.Scrub2 %in% "Wet Lime FGD") |  
                   (dat$SO2.Scrub3 %in% "Wet Lime FGD") | 
                   (dat$SO2.Scrub4 %in% "Wet Lime FGD")] = TRUE
  print("Variable Created: WetLimeFGD")
  #with(dat, mean(WetLimeFGD))
  
  dat$WetLime = FALSE
  dat$WetLime[(dat$SO2.Scrub1 == "Wet Limestone") | 
                (dat$SO2.Scrub2 %in% "Wet Limestone") |  
                (dat$SO2.Scrub3 %in% "Wet Limestone") | 
                (dat$SO2.Scrub4 %in% "Wet Limestone")] = TRUE
  print("Variable Created: WetLime")
  #with(dat, mean(WetLime))
  
  dat$OtherSO2 = FALSE
  dat$OtherSO2[(dat$SO2.Scrub1 == "Other") | 
                 (dat$SO2.Scrub2 %in% "Other") |  
                 (dat$SO2.Scrub3 %in% "Other") | 
                 (dat$SO2.Scrub4 %in% "Other")] = TRUE
  print("Variable Created: OtherSO2")
  #with(dat, mean(OtherSO2))
  
  return(dat)
}

region <- function(countyfips){
  # Author: Cory Zigler
  # Converts FIPS to any of 7 regions used in EPA's 1996 PM Criteria Doc.
  # Input: 
  #   countyfips: FIPS codes.
  # Output:
  #   Region specified by FIPS.
  
  data(state.fips) # from maps
  
  # Define 7 regions by state names
  Northeast = c("ME", "NH", "VT", "NY", "PA", "DE", "NJ", "MD", "DC", 
                "VA", "MA", "CT", "RI")
  IndustrialMidwest = c("WV", "OH", "KY", "IN", "IL", "WI", "MI")
  Southeast = c("FL", "GA", "SC", "NC", "TN", "AL", "MS", "AR","LA")
  UpperMidwest = c("MN", "IA", "MO", "KS", "NE", "SD", "ND")
  Southwest = c("TX", "OK", "NM", "AZ")
  SouthernCalifornia = c("CA")
  Northwest = c("NV", "UT", "CO", "WY", "MT", "ID", "OR", "WA")
  
  # Define 7 regions by state FIPS codes (approximation)
  Northeast.fips = unique(state.fips$fips[state.fips$abb %in% Northeast])
  IndustrialMidwest.fips = unique(state.fips$fips[state.fips$abb %in% IndustrialMidwest])
  Southeast.fips = unique(state.fips$fips[state.fips$abb %in% Southeast])
  UpperMidwest.fips = unique(state.fips$fips[state.fips$abb %in% UpperMidwest])
  Southwest.fips = unique(state.fips$fips[state.fips$abb %in% Southwest])
  SouthernCalifornia.fips = unique(state.fips$fips[state.fips$abb %in% SouthernCalifornia])
  Northwest.fips = unique(state.fips$fips[state.fips$abb %in% Northwest])
  
  # Note: If FIPS codes have numerical values (i.e., some have 4 digits 
  #   without leading 0's), then run a following line.
  # countyfips <- as.character(formatC(countyfips, format="d", width=5, flag="0"))
  
  # convert county fips to regions
  State_FIPS <- as.numeric(substr(countyfips, 1, 2))
  new_region <- ifelse(State_FIPS %in% Northeast.fips, "Northeast",
    ifelse(State_FIPS %in% IndustrialMidwest.fips, "IndustrialMidwest",
    ifelse(State_FIPS %in% Southeast.fips, "Southeast",
    ifelse(State_FIPS %in% UpperMidwest.fips, "UpperMidwest",
    ifelse(State_FIPS %in% Southwest.fips, "Southwest",
    ifelse(State_FIPS %in% SouthernCalifornia.fips, "SouthernCalifornia", 
    "Northwest"))))))
  
  return(new_region)
}

cleanData <- function(dat){
  # This function is used to clean the "MonthlyUnitData.RDS"
  #   (coal-fired power plant facility) data.  The output includes the 
  #   same variables as the annual unit data.
  # Input:
  #   dat: data frame from "MonthlyUnitData.Rda"
  # Output:
  #   Cleaned data frame, only keeping important variables.
  
  # make pctCapacity variable:
  pct.capacity <- dat$Heat.Input..MMBtu./ dat$Capacity
  
  # # indices of columns to keep
  # inds1 <- c(83, 4, 9, 10, 46, 8, 106, 105, 25, 18, 19, 21, 23, 22)
  # inds2 <- c(51, 16, 17, 15, 86)
  
  # create new data frame with important columns
  new.dat <- cbind(dat[, c(83,4,9,10,46,8,106,105,25,18,19,21,23,22)], 
                   pct.capacity, 
                   dat[, c(51,16,17,15,86)])
  
  # rename columns
  colnames(new.dat) <- c("uID",
                         "Year",
                         "Month",
                         "Sulfur.Content",
                         "Has.Scrubber",
                         "Initial.Year.of.Operation",
                         "S_n_CR",
                         "NumNOxControls",
                         "Operating.Time",
                         "SO2..tons.",
                         "NOx..tons.",
                         "CO2..short.tons.",
                         "Gross.Load..MW.h.",
                         "Heat.Input..MMBtu.",
                         "pctCapacity",
                         "Is.Phase2",
                         "Facility.Latitude.x",
                         "Facility.Longitude.x",
                         "FIPS",
                         "Region")
  
  # return new data frame
  new.dat
}

unitNames <- function(annual = FALSE){
  # Create vector of column names for the emissions data, to better match 
  #   facility variable names.
  # Input: 
  #   annual: logical flag indicating if this is performed to be performed
  #     for annual or monthly data.
  # Output:
  #   Vector of column names.
  
  # column names
  c.names <- c("uID",
               "Year",
               "Month",
               "SulfurContent",
               "ScrubbedUnit",
               "initialYear",
               "S_n_CR",
               "NumNOxControls",
               "OpTime",
               "SO2emissions",
               "NOxemissions",
               "CO2emissions",
               "Load",
               "HeatInput",
               "pctCapacity",
               "Phase2",
               "Fac.Latitude",
               "Fac.Longitude",
               "Fac.FIPS",
               "Fac.Region")
  
  # remove "Month" name if data are annual
  if (annual){
    c.names <- c.names[-c(3)]
  }
  
  # return column names
  c.names
}

trimData <- function(dat, ex = c(-100, -81.5, 30.5, 41.7)){
  # Trim emissions facilities data to only include facilities within a given
  #   spatial extent.
  # Input: 
  #   dat: data frame with emissions data.
  #   ex: spatial extent (default = central USA).
  # Output
  #   Trimmed data frame, with unit and fac level data.

  # identify emissions units within spatial extent
  idx.unit <- which(dat$unit$Fac.Longitude > ex[1]
                    & dat$unit$Fac.Longitude < ex[2]
                    & dat$unit$Fac.Latitude > ex[3]
                    &  dat$unit$Fac.Latitude < ex[4])
  
  # identify emissions facilities within spatial extent
  idx.fac <- which(dat$fac$Fac.Longitude > ex[1]
                   & dat$fac$Fac.Longitude < ex[2]
                   & dat$fac$Fac.Latitude > ex[3]
                   &  dat$fac$Fac.Latitude < ex[4])
  
  # create new list of data frames (units and fac) found within ex
  final.dat <- list(unit = dat$unit[idx.unit,],
                    fac = dat$fac[idx.fac,])
  
  # return final.dat list
  final.dat
}

trimYear <- function(dat, year){
  # Trim emissions data to the specified year.
  # Input:
  #   dat: emissions data frame.
  #   year: year of interest.
  # Output:
  #   Data frame with emissions data for the specified year.
  
  # indices for units and facilities
  idx.unit <- which(dat$unit$Year == year)
  idx.fac <- which(dat$fac$Year == year)
  
  # new list of data frames, restricted to year
  new.dat <- list(unit = dat$unit[idx.unit,],
                  fac = dat$fac[idx.fac,])
  
  # return new list
  new.dat
}

createSimpleX <- function(r, emissions){
  # Create the X design matrix for the pollution data.
  # Input:
  #   r: a single raster object
  #   emissions: a data frame with power plant emissions data.
  # Output:
  #   A design vector, X, representing the total amount of emissions in 
  #     each grid of the raster.
  
  # number of cells in rectangle
  n.cells <- length(values(r))
  
  # create an SO2 variable using the emissions data
  em.vals <- rep(0, n.cells)
  em.inds <- raster::cellFromXY(r, cbind(emissions$Fac.Longitude, emissions$Fac.Latitude))
  
  # make sure unique inds are
  unique.inds <- unique(em.inds)
  
  for(curr.ind in unique.inds){
    spec.points <- em.inds == curr.ind
    em.vals[curr.ind] <- sum(emissions$totSO2emissions[spec.points])
  }
  
  #em.vals[em.inds] <- emissions$totSO2emissions
  
  # combine all values into a single vector
  X <- em.vals
  
  # return X
  return(X)
}

createMexicoX <- function(r, emissions){
  # Create the X design matrix for the pollution data.
  # Input:
  #   r: a single raster object
  #   emissions: a data frame with power plant emissions data.
  # Output:
  #   A design vector, X, representing the total amount of emissions in
  #     each grid of the raster.

  # number of cells in rectangle
  n.cells <- length(values(r))

  # create an SO2 variable using the emissions data
  em.vals <- rep(0, n.cells)
  em.inds <- raster::cellFromXY(r, cbind(emissions$longitude, emissions$latitude))

  # make sure unique inds are
  unique.inds <- unique(em.inds)

  for(curr.ind in unique.inds){
    spec.points <- em.inds == curr.ind
    em.vals[curr.ind] <- sum(emissions$total.SO2[spec.points]) / 12
  }

  #em.vals[em.inds] <- emissions$totSO2emissions

  # combine all values into a single vector
  X <- em.vals

  # return X
  return(X)
}

### B. Cleaning meteorology/landscape covariate data

yearRaster <- function(year.in, met = "temp", return.mean = TRUE){
  # Extract the yearly mean values for a specified meteorology raster.
  # Input:
  #   year.in: specifies which year of data to grab (e.g., 2011).
  #   met: type of meteorological data 
  #     (options include "temp", "precip", "humid", "v-wind", "uwind").
  # Output:
  #   Raster with specified meteorological data.
  
  # specify the input file
  if (met == "temp"){
    # temperature at 2m
    file.in <- "./data/air.2m.mon.mean.nc"
  } else if (met == "precip"){
    # total precip
    file.in <- "./data/apcp.mon.mean.nc"
  } else if (met == "humid"){
    # relative humidity at 2m
    file.in <- "./data/rhum.2m.mon.mean.nc"
  } else if (met == "v-wind"){
    # v (north-south) wind component
    file.in <- "./data/vwnd.10m.mon.mean.nc"
  } else if (met == "u-wind"){
    # u (east-west) wind component
    file.in <- "./data/uwnd.10m.mon.mean.nc"
  }
  
  # determine which raster bands to grab
  bands <- ((year.in - 1979) * 12) + 1:12
  
  # objects to save rasters and month values
  raster.month <- list()
  month.vals <- matrix(NA_real_, nrow = 96673, ncol = 12)
  
  for (i in 1:length(bands)){
    # grab raster for each month in the year of interest
    raster.month[[i]] <- raster(file.in, band = bands[i])
    month.vals[,i] <- values(raster.month[[i]])
  }
  
  if (met == "temp"){
    # convert temps to celsius
    month.vals <- month.vals - 273.15
  }
  
  if (return.mean){
  # take the annual mean for the given year
    annual.mean <- rowMeans(month.vals)
    annual.mean[is.na(annual.mean)] <- NA

    # save results in a raster
    annual.raster <- raster.month[[1]]
    values(annual.raster) <- annual.mean
  
    # return raster object
    return(annual.raster)

  } else {
    # return monthly data
    return(raster.month)
  }
  
}

trimToSO4 <- function(new.raster, so4.raster){
  # Trim environmental raster to match resolution and projection 
  #   of so4 raster.
  # Input:
  #   new.raster: raster with environmental covariate (to be changed).
  #   so4.raster: raster with so4 data (default res/proj).
  # Output:
  #   Environmental raster with res. and projection matching so4.raster.
  
  # change raster projection of SO4 to match new raster
  so <- projectRaster(so4.raster, 
                      crs = crs(projection(new.raster)), 
                      method = "bilinear")
  
  # crop new raster to match SO4 etent
  new.rs <- crop(new.raster, extent(so))
  
  # convert new raster to SO4 resolution and projection
  newc <- so
  values(newc) <- NA
  xy <- xyFromCell(so, 1:ncell(so))
  idx <- cellFromXY(new.rs, xy)
  new.vals <- values(new.rs)[idx]
  values(newc) <- new.vals
  new.rstr <- projectRaster(newc, 
                            crs = crs(projection(so4.raster)), 
                            method = "bilinear")
  
  # transfer values between non matching raster objects
  final <- resample(new.rstr, so4.raster)
  
  # return new raster object
  return(final)
}

################################################################################
### 2. Functions used in sulfate analysis.
################################################################################
 

### A. Create Finite Volume Method structures for given raster

FVMInds <- function(dims, bound.cond){
  # Create list denoting index for interior and boundary of raster with given
  #   spatial extent. This is used to create the FVM approximation of advection
  #   due to wind.
  # Input:
  #   dims: dimension of raster
  #   bound.cond: type of boundary condition (insulated or periodic).
  # Output:
  #   List with index of interior and boundary points.
  
  # rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]
  
  #######################
  ### interior points ###
  #######################
  interior <- matrix(NA_real_, nrow = (n.c - 2) * (n.r - 2), ncol = 5)
  colnames(interior) <- c("P", "N", "S", "E", "W")
  
  if (bound.cond == "insulated"){
    # create a diffusion matrix with insulated boundaries
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){
        
        # index in original scheme
        orig.ind <- (i - 1) * n.c + j
        
        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1
        
        cur.ind <- cur.ind + 1
      }
    }
    
    ### west boundary points
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(west.bound) <- c("P", "N", "S", "E")
    
    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### east boundary
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(east.bound) <- c("P", "N", "S", "W")
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      
      # index in original scheme
      orig.ind <- i * n.c
      
      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### north boundary
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(north.bound) <- c("P", "S", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j
      
      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + n.c
      north.bound[cur.ind, 3] <- orig.ind + 1
      north.bound[cur.ind, 4] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### south boundary
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(south.bound) <- c("P", "N", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)
      
      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind + 1
      south.bound[cur.ind, 4] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### corners
    nw <- c(1, NA_real_, 1 + n.c, 2, NA_real_)
    ne <- c(n.c, NA_real_, 2 * n.c, NA_real_, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, NA_real_, (n.r - 1) * n.c + 2, NA_real_)
    se <- c(n.r * n.c, (n.r - 1) * n.c, NA_real_, NA_real_, n.r * n.c - 1)
    
    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")
    
  } else if (bound.cond == "periodic"){
    
    # create a diffusion matrix with periodic boundaries
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){
        
        # index in original scheme
        orig.ind <- (i - 1) * n.c + j
        
        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1
        
        cur.ind <- cur.ind + 1
      }
    }
    
    ############################
    ### west boundary points ###
    ############################
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(west.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1
      west.bound[cur.ind, 5] <- orig.ind + (n.c - 1)
      
      cur.ind <- cur.ind + 1
    }
    
    #####################
    ### east boundary ###
    #####################
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(east.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      
      # index in original scheme
      orig.ind <- i * n.c
      
      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - (n.c - 1)
      east.bound[cur.ind, 5] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ######################
    ### north boundary ###
    ######################
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(north.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j
      
      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + (n.r - 1) * n.c
      north.bound[cur.ind, 3] <- orig.ind + n.c
      north.bound[cur.ind, 4] <- orig.ind + 1
      north.bound[cur.ind, 5] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ######################
    ### south boundary ###
    ######################
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(south.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)
      
      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind - (n.r - 1) * n.c
      south.bound[cur.ind, 4] <- orig.ind + 1
      south.bound[cur.ind, 5] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ###############
    ### corners ###
    ###############
    nw <- c(1, 1 + (n.r - 1) * n.c, 1 + n.c, 2, n.c)
    ne <- c(n.c, n.c * n.r, 2 * n.c, 1, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, 1, (n.r - 1) * n.c + 2, n.r * n.c)
    se <- c(n.r * n.c, (n.r - 1) * n.c, n.c, 1 + (n.r - 1) * n.c, n.r * n.c - 1)
    
    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")
  }
  
  # combine all indice matrices into a single list
  indices <- list(interior = interior,
                  north = north.bound,
                  south = south.bound,
                  east = east.bound,
                  west = west.bound,
                  corners = corners)
  
  # return indices
  return(indices)
}


diffFVM <- function(dims, inds, bound.cond){
  # Create finite volume matrix approximation to homogeneous diffusion.
  # Input:
  #   dims: dimension of raster
  #   inds: indices of interior/boundary points, created with FVMInds
  #   bound.cond: boundary conditions (insulated, periodic, Neumann)
  # Output:
  #   Sparse matrix used for discrete approximation of homogenous diffusion.

  # rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]

  # total number of cells
  n.cells <- n.r * n.c

  if (bound.cond == "insulated"){

    # i, j, and x values for interior points
    n.interior <- nrow(inds$interior)
    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(rep(4, n.interior), rep(-1, n.interior * 4)))

    # i, j, and x values for west boundary points
    n.west <- nrow(inds$west)
    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west),
                         c(rep(3, n.west), rep(-1, n.west * 3)))

    # i, j, and x values for east boundary points
    n.east <- nrow(inds$east)
    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east),
                         c(rep(3, n.east), rep(-1, n.east * 3)))

    # i, j, and x values for north boundary points
    n.north <- nrow(inds$north)
    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north),
                          c(rep(3, n.north), rep(-1, n.north * 3)))

    # i, j, and x values for south boundary points
    n.south <- nrow(inds$south)
    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south),
                          c(rep(3, n.south), rep(-1, n.south * 3)))

    # i, j, and x values for corner points

    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(rep(2, 4), rep(-1, 8)))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

  } else if (bound.cond == "periodic"){

    # i, j, and x values for interior points
    n.interior <- nrow(inds$interior)
    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(rep(4, n.interior), rep(-1, n.interior * 4)))


    # i, j, and x values for west boundary points
    n.west <- nrow(inds$west)
    sparse.west <- cbind(rep(inds$west[,1], 5),
                         as.vector(inds$west),
                         c(rep(4, n.west), rep(-1, n.west * 4)))

    # i, j, and x values for east boundary points
    n.east <- nrow(inds$east)
    sparse.east <- cbind(rep(inds$east[,1], 5),
                         as.vector(inds$east),
                         c(rep(4, n.east), rep(-1, n.east * 4)))

    # i, j, and x values for north boundary points
    n.north <- nrow(inds$north)
    sparse.north <- cbind(rep(inds$north[,1], 5),
                          as.vector(inds$north),
                          c(rep(4, n.north), rep(-1, n.north * 4)))

    # i, j, and x values for south boundary points
    n.south <- nrow(inds$south)
    sparse.south <- cbind(rep(inds$south[,1], 5),
                          as.vector(inds$south),
                          c(rep(4, n.south), rep(-1, n.south * 4)))

    # i, j, and x values for corner points
    sparse.corner <- cbind(rep(inds$corners[,1], 5),
                           as.vector(inds$corners),
                           c(rep(4, 4), rep(-1, 4 * 4)))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

  } else if (bound.cond == "Neumann"){

    # i, j, and x values for interior points
    n.interior <- nrow(inds$interior)
    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(rep(4, n.interior), rep(-1, n.interior * 4)))

    # i, j, and x values for west boundary points
    n.west <- nrow(inds$west)
    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west),
                         c(rep(4, n.west), rep(-1, n.west * 2), rep(-2, n.west)))

    # i, j, and x values for east boundary points
    n.east <- nrow(inds$east)
    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east),
                         c(rep(4, n.east), rep(-1, n.east * 2), rep(-2, n.east)))

    # i, j, and x values for north boundary points
    n.north <- nrow(inds$north)
    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north),
                          c(rep(4, n.north), rep(-2, n.north), rep(-1, n.north * 2)))

    # i, j, and x values for south boundary points
    n.south <- nrow(inds$south)
    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south),
                          c(rep(4, n.south), rep(-2, n.south), rep(-1, n.south * 2)))

    # i, j, and x values for corner points
    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(rep(4, 4), rep(-2, 8)))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)
  }

  # quickly fill the sparse matrix
  D.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
                        dims = c(n.cells, n.cells))

  # return the sparse diffusion matrix
  return(D.mat)
}


FVMInds <- function(dims, bound.cond){
  # Create list denoting index for interior and boundary of raster with given
  #   spatial extent. This is used to create the FVM approximation of advection
  #   due to wind.
  # Input:
  #   dims: dimension of raster
  #   bound.cond: type of boundary condition (insulated or periodic).
  # Output:
  #   List with index of interior and boundary points.

  # rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]

  #######################
  ### interior points ###
  #######################
  interior <- matrix(NA_real_, nrow = (n.c - 2) * (n.r - 2), ncol = 5)
  colnames(interior) <- c("P", "N", "S", "E", "W")

  if (bound.cond == "insulated"){
    # create a diffusion matrix with insulated boundaries

    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){

        # index in original scheme
        orig.ind <- (i - 1) * n.c + j

        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1

        cur.ind <- cur.ind + 1
      }
    }

    ### west boundary points
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(west.bound) <- c("P", "N", "S", "E")

    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1

      cur.ind <- cur.ind + 1
    }

    ### east boundary
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(east.bound) <- c("P", "N", "S", "W")

    cur.ind <- 1
    for (i in 2:(n.r-1)){

      # index in original scheme
      orig.ind <- i * n.c

      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ### north boundary
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(north.bound) <- c("P", "S", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j

      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + n.c
      north.bound[cur.ind, 3] <- orig.ind + 1
      north.bound[cur.ind, 4] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ### south boundary
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(south.bound) <- c("P", "N", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)

      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind + 1
      south.bound[cur.ind, 4] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ### corners
    nw <- c(1, NA_real_, 1 + n.c, 2, NA_real_)
    ne <- c(n.c, NA_real_, 2 * n.c, NA_real_, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, NA_real_, (n.r - 1) * n.c + 2, NA_real_)
    se <- c(n.r * n.c, (n.r - 1) * n.c, NA_real_, NA_real_, n.r * n.c - 1)

    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")

  } else if (bound.cond == "periodic"){

    # create a diffusion matrix with periodic boundaries

    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){

        # index in original scheme
        orig.ind <- (i - 1) * n.c + j

        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1

        cur.ind <- cur.ind + 1
      }
    }

    ############################
    ### west boundary points ###
    ############################
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(west.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1
      west.bound[cur.ind, 5] <- orig.ind + (n.c - 1)

      cur.ind <- cur.ind + 1
    }

    #####################
    ### east boundary ###
    #####################
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(east.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (i in 2:(n.r-1)){

      # index in original scheme
      orig.ind <- i * n.c

      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - (n.c - 1)
      east.bound[cur.ind, 5] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ######################
    ### north boundary ###
    ######################
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(north.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j

      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + (n.r - 1) * n.c
      north.bound[cur.ind, 3] <- orig.ind + n.c
      north.bound[cur.ind, 4] <- orig.ind + 1
      north.bound[cur.ind, 5] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ######################
    ### south boundary ###
    ######################
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(south.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)

      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind - (n.r - 1) * n.c
      south.bound[cur.ind, 4] <- orig.ind + 1
      south.bound[cur.ind, 5] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ###############
    ### corners ###
    ###############
    nw <- c(1, 1 + (n.r - 1) * n.c, 1 + n.c, 2, n.c)
    ne <- c(n.c, n.c * n.r, 2 * n.c, 1, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, 1, (n.r - 1) * n.c + 2, n.r * n.c)
    se <- c(n.r * n.c, (n.r - 1) * n.c, n.c, 1 + (n.r - 1) * n.c, n.r * n.c - 1)

    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")
  }

  # combine all indice matrices into a single list
  indices <- list(interior = interior,
                  north = north.bound,
                  south = south.bound,
                  east = east.bound,
                  west = west.bound,
                  corners = corners)

  # return indices
  return(indices)
}


interpWind <- function(r.old, r.new, direction, with.bound = TRUE){
  # Linear interpolation of wind velocity across surface boundaries.
  # Input:
  #   r.old: wind raster with larger spatial extent.
  #   r.new: wind raster with smaller spatial extent.
  #   direction: direction of wind (either "u" or "v").
  #   with.bound: boolean determining if velocities should be saved on
  #     the boundary.
  # Output:
  #   Matrix with interpolated matrix velocities.

  # grab xy coordinates from rasters
  xy.old <- coordinates(r.old)
  xy.new <- coordinates(r.new)

  # dimensions of rasters (n.rows by n.cols)
  dim.old <- dim(r.old)[1:2]
  dim.new <- dim(r.new)[1:2]

  # cell indices for trimmed case
  cell.inds <- cellFromXY(r.old, xy.new)

  # indices for the top boundary
  top.bound <- (cell.inds[1] - dim.old[2]) + (-1):(dim.new[2])

  # cell indices with boundaries
  all.bounds <- top.bound

  for (k in 1:(dim.new[1] + 1)){
    next.bound <- top.bound + (k * dim.old[2])
    all.bounds <- c(all.bounds, next.bound)
  }

  # wind values for all points (including boundaries)
  orig.vals <- values(r.old)[all.bounds]
  # matrix of face value velocities
  face.vals <- matrix(NA_real_, nrow = length(all.bounds), ncol = 2)

  # new dimensions
  n.rows <- dim.new[1] + 2
  n.cols <- dim.new[2] + 2

  if (direction == "u"){
    ### west to east

    # name columns appropriately
    colnames(face.vals) <- c("west", "east")

    for(i in 1:n.rows){
      for(j in 1:n.cols){

        # current indice
        cur.ind <- (i - 1) * n.cols + j

        if (j == 1){
          # no west face velocity
          west <- 0
          east <- sum(orig.vals[cur.ind + 0:1]) / 2
        } else if (j == n.cols){
          # no east face velocity
          west <-  sum(orig.vals[cur.ind + -1:0]) / 2
          east <- 0
        } else {
          # interior point
          west <- sum(orig.vals[cur.ind + -1:0]) / 2
          east <- sum(orig.vals[cur.ind + 0:1]) / 2
        }

        # add values to boundary matrix
        face.vals[cur.ind,] <- c(west, east)
      }
    }
  } else if (direction == "v"){
    # south to north
    # name columns appropriately
    colnames(face.vals) <- c("south", "north")

    for(i in 1:n.rows){
      for(j in 1:n.cols){

        # current indice
        cur.ind <- (i - 1) * n.cols + j

        if (i == 1){
          # no north face velocity
          north <- 0
          south <- (orig.vals[cur.ind] + orig.vals[cur.ind + n.cols]) / 2
        } else if (i == n.rows){
          # no south face velocity
          north <-  (orig.vals[cur.ind] + orig.vals[cur.ind - n.cols]) / 2
          south <- 0
        } else {
          # interior point
          north <-  (orig.vals[cur.ind] + orig.vals[cur.ind - n.cols]) / 2
          south <- (orig.vals[cur.ind] + orig.vals[cur.ind + n.cols]) / 2
        }

        # add values to boundary matrix
        face.vals[cur.ind,] <- c(south, north)
      }
    }
  }

  if (with.bound){
    # return face.vals as is
    return(face.vals)

  } else {
    # remove boundary points from face.vals
    bound.points <- c(1:n.cols,
                      (1:(n.rows-2) * n.cols + 1),
                      (2:(n.rows-1) * n.cols),
                      ((n.rows-1) * n.cols + 1:n.cols))

    return(face.vals[-bound.points,])
  }
}


windConvFVM <- function(dims, velocity.faces, inds, bc, bound = FALSE){
  # Create an advection matrix using the Finite Volume Method (FVM).
  # Input:
  #   dims: dimension of the raster
  #   velocity.faces: interpolated velocities on grid edges (from interpWind)
  #   inds: matrix indices of raster grid cells (from FVMinds).
  #   bc: boundary condition
  #   bound: not used?
  # Output:
  #   Sparse FVM matrix representation of advection due to wind.

  ## rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]

  # total number of cells
  n.cells <- n.r * n.c

  if (bc == "insulated"){

    ### i, j, and x values for interior points
    n.interior <- nrow(inds$interior)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$interior[,1], 2] -
      velocity.faces[inds$interior[,1], 1]
    F.diff.ns <- velocity.faces[inds$interior[,1], 4] -
      velocity.faces[inds$interior[,1], 3]
    F.diff.comb <- F.diff.ew + F.diff.ns

    a.w <- velocity.faces[inds$interior[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$interior[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$interior[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$interior[,1], 4]; a.n[a.n < 0] <- 0

    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(F.diff.comb + a.n + a.s + a.e + a.w, -a.n, -a.s, -a.e, -a.w))


    ### i, j, and x values for west boundary points
    n.west <- nrow(inds$west)

    F.diff.ns <- velocity.faces[inds$west[,1], 4] -
      velocity.faces[inds$west[,1], 3]
    F.diff.comb <- F.diff.ns + velocity.faces[inds$west[,1], 2]

    a.e <- -velocity.faces[inds$west[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$west[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- velocity.faces[inds$west[,1], 4]; a.n[a.n > 0] <- 0

    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west[,-5]),
                         c(F.diff.comb + a.n + a.s + a.e, -a.n, -a.s, -a.e))

    ### i, j, and x values for east boundary points
    n.east <- nrow(inds$east)

    # coefficient calculations
    F.diff.ns <- velocity.faces[inds$east[,1], 4] -
      velocity.faces[inds$east[,1], 3]
    F.diff.comb <- F.diff.ns - velocity.faces[inds$east[,1], 1]

    a.w <- velocity.faces[inds$east[,1], 1]; a.w[a.w < 0] <- 0
    a.s <- velocity.faces[inds$east[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$east[,1], 4]; a.n[a.n < 0] <- 0

    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east[,-4]),
                         c(F.diff.comb + a.n + a.s + a.w, -a.n, -a.s, -a.w))

    ### i, j, and x values for north boundary points
    n.north <- nrow(inds$north)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$north[,1], 2] -
      velocity.faces[inds$north[,1], 1]
    F.diff.comb <- F.diff.ew - velocity.faces[inds$north[,1], 3]

    a.w <- velocity.faces[inds$north[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$north[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$north[,1], 3]; a.s[a.s < 0] <- 0

    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north[,-2]),
                          c(F.diff.comb + a.s + a.e + a.w, -a.s, -a.e, -a.w))

    ### i, j, and x values for south boundary points
    n.south <- nrow(inds$south)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$south[,1], 2] -
      velocity.faces[inds$south[,1], 1]
    F.diff.comb <- F.diff.ew + velocity.faces[inds$south[,1], 4]

    a.w <- velocity.faces[inds$south[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$south[,1], 2]; a.e[a.e < 0] <- 0
    a.n <- -velocity.faces[inds$south[,1], 4]; a.n[a.n < 0] <- 0

    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south[,-3]),
                          c(F.diff.comb + a.n + a.e + a.w, -a.n, -a.e, -a.w))

    ### i, j, and x values for corner points

    # northwest values
    nw.ae <- max(c(0, -velocity.faces[inds$corners[1,1], 2]))
    nw.as <- max(c(0, velocity.faces[inds$corners[1,1], 3]))
    nw.ap <- nw.ae + nw.as + velocity.faces[inds$corners[1,1], 2] -
      velocity.faces[inds$corners[1,1], 3]

    # northeast values
    ne.aw <- max(c(0, velocity.faces[inds$corners[2,1], 1]))
    ne.as <- max(c(0, velocity.faces[inds$corners[2,1], 3]))
    ne.ap <- ne.aw + ne.as - velocity.faces[inds$corners[2,1], 1] -
      velocity.faces[inds$corners[2,1], 3]

    # southwest values
    sw.ae <- max(c(0, -velocity.faces[inds$corners[3,1], 2]))
    sw.an <- max(c(0, -velocity.faces[inds$corners[3,1], 4]))
    sw.ap <- sw.ae + sw.an + velocity.faces[inds$corners[3,1], 2] +
      velocity.faces[inds$corners[3,1], 4]

    # southeast values
    se.aw <- max(c(0, velocity.faces[inds$corners[4,1], 1]))
    se.an <- max(c(0, -velocity.faces[inds$corners[4,1], 4]))
    se.ap <- se.aw + se.an - velocity.faces[inds$corners[4,1], 1] +
      velocity.faces[inds$corners[4,1], 4]

    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(nw.ap, ne.ap, sw.ap, se.ap, -sw.an, -se.an, -nw.as,
                             -ne.as, -nw.ae, -sw.ae, -ne.aw, -se.aw))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

    # quickly fill the sparse matrix
    A.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
                          dims = c(n.cells, n.cells))

  } else if (bc == "endless world"){

    ### i, j, and x values for interior points
    n.interior <- nrow(inds$interior)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$interior[,1], 2] -
      velocity.faces[inds$interior[,1], 1]
    F.diff.ns <- velocity.faces[inds$interior[,1], 4] -
      velocity.faces[inds$interior[,1], 3]
    F.diff.comb <- F.diff.ew + F.diff.ns

    a.w <- velocity.faces[inds$interior[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$interior[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$interior[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$interior[,1], 4]; a.n[a.n < 0] <- 0

    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(F.diff.comb + a.n + a.s + a.e + a.w, -a.n, -a.s, -a.e, -a.w))


    ### i, j, and x values for west boundary points
    n.west <- nrow(inds$west)

    # coefficient values
    F.diff.ns <- velocity.faces[inds$west[,1], 4] -
      velocity.faces[inds$west[,1], 3]
    F.diff.comb <- F.diff.ns + velocity.faces[inds$west[,1], 2]

    a.e <- -velocity.faces[inds$west[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$west[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- velocity.faces[inds$west[,1], 4]; a.n[a.n > 0] <- 0

    west.faces <- velocity.faces[inds$west[,1], 1]
    diag.vals <- F.diff.comb + a.n + a.s + a.e
    diag.vals[west.faces < 0] <- (diag.vals - west.faces)[west.faces < 0]

    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west[,-5]),
                         c(diag.vals, -a.n, -a.s, -a.e))

    ### i, j, and x values for east boundary points
    n.east <- nrow(inds$east)

    # coefficient calculations
    F.diff.ns <- velocity.faces[inds$east[,1], 4] -
      velocity.faces[inds$east[,1], 3]
    F.diff.comb <- F.diff.ns - velocity.faces[inds$east[,1], 1]

    a.w <- velocity.faces[inds$east[,1], 1]; a.w[a.w < 0] <- 0
    a.s <- velocity.faces[inds$east[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$east[,1], 4]; a.n[a.n < 0] <- 0

    # if east face value is positive, change diagonal value to include it leaving system
    east.faces <- velocity.faces[inds$east[,1],2]
    diag.vals <- F.diff.comb + a.n + a.s + a.w
    diag.vals[east.faces > 0] <- (diag.vals + east.faces)[east.faces > 0]

    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east[,-4]),
                         c(diag.vals, -a.n, -a.s, -a.w))

    ### i, j, and x values for north boundary points
    n.north <- nrow(inds$north)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$north[,1], 2] -
      velocity.faces[inds$north[,1], 1]
    F.diff.comb <- F.diff.ew - velocity.faces[inds$north[,1], 3]

    a.w <- velocity.faces[inds$north[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$north[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$north[,1], 3]; a.s[a.s < 0] <- 0

    north.faces <- velocity.faces[inds$north[,1], 4]
    diag.vals <- F.diff.comb + a.s + a.e + a.w
    diag.vals[north.faces > 0] <- (diag.vals + north.faces)[north.faces > 0]

    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north[,-2]),
                          c(diag.vals, -a.s, -a.e, -a.w))

    ### i, j, and x values for south boundary points
    n.south <- nrow(inds$south)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$south[,1], 2] -
      velocity.faces[inds$south[,1], 1]
    F.diff.comb <- F.diff.ew + velocity.faces[inds$south[,1], 4]

    a.w <- velocity.faces[inds$south[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$south[,1], 2]; a.e[a.e < 0] <- 0
    a.n <- -velocity.faces[inds$south[,1], 4]; a.n[a.n < 0] <- 0

    south.faces <- velocity.faces[inds$south[,1], 3]
    diag.vals <- F.diff.comb + a.n + a.e + a.w
    diag.vals[south.faces < 0] <- (diag.vals - south.faces)[south.faces < 0]

    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south[,-3]),
                          c(F.diff.comb + a.n + a.e + a.w, -a.n, -a.e, -a.w))

    ### i, j, and x values for corner points

    ## northwest values
    nw.ae <- max(c(0, -velocity.faces[inds$corners[1,1], 2]))
    nw.as <- max(c(0, velocity.faces[inds$corners[1,1], 3]))
    nw.ap <- nw.ae + nw.as + velocity.faces[inds$corners[1,1], 2] -
      velocity.faces[inds$corners[1,1], 3]

    # north face
    nw.nf <- velocity.faces[inds$corners[1,1], 4]
    # west face
    nw.wf <- velocity.faces[inds$corners[1,1], 1]

    if (nw.nf > 0){
      # add to diagonal
      nw.ap <- nw.ap + nw.nf
    }

    if (nw.wf < 0){
      # subtract from diagonal
      nw.ap <- nw.ap - nw.wf
    }

    ## northeast values
    ne.aw <- max(c(0, velocity.faces[inds$corners[2,1], 1]))
    ne.as <- max(c(0, velocity.faces[inds$corners[2,1], 3]))
    ne.ap <- ne.aw + ne.as - velocity.faces[inds$corners[2,1], 1] -
      velocity.faces[inds$corners[2,1], 3]

    # north face
    ne.nf <- velocity.faces[inds$corners[2,1], 4]
    # east face
    ne.ef <- velocity.faces[inds$corners[2,1], 2]

    if (ne.nf > 0){
      ne.ap <- ne.ap + ne.nf
    }

    if (ne.ef > 0){
      ne.ap <- ne.ap + ne.ef
    }


    ## southwest values
    sw.ae <- max(c(0, -velocity.faces[inds$corners[3,1], 2]))
    sw.an <- max(c(0, -velocity.faces[inds$corners[3,1], 4]))
    sw.ap <- sw.ae + sw.an + velocity.faces[inds$corners[3,1], 2] +
      velocity.faces[inds$corners[3,1], 4]

    # south face
    sw.sf <- velocity.faces[inds$corners[3,1],3]
    # west face
    sw.wf <- velocity.faces[inds$corners[3,1],1]

    if (sw.sf < 0){
      sw.ap <- sw.ap - sw.sf
    }

    if (sw.wf < 0){
      sw.ap <- sw.ap - sw.wf
    }

    ## southeast values
    se.aw <- max(c(0, velocity.faces[inds$corners[4,1], 1]))
    se.an <- max(c(0, -velocity.faces[inds$corners[4,1], 4]))
    se.ap <- se.aw + se.an - velocity.faces[inds$corners[4,1], 1] +
      velocity.faces[inds$corners[4,1], 4]

    # south face
    se.sf <- velocity.faces[inds$corners[4,1], 3]
    # east face
    se.ef <- velocity.faces[inds$corners[4,1], 2]

    if (se.sf < 0){
      se.ap <- se.ap - se.sf
    }

    if (se.ef > 0){
      se.ap <- se.ap + se.ef
    }

    # combine all values together
    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(nw.ap, ne.ap, sw.ap, se.ap, -sw.an, -se.an, -nw.as,
                             -ne.as, -nw.ae, -sw.ae, -ne.aw, -se.aw))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

    # quickly fill the sparse matrix
    A.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
                          dims = c(n.cells, n.cells))
  }

  return(A.mat)
}

createC <- function(big, small, dims, inds, bound.cond = "endless world"){
  # Create matrix fo discretized advection process, due to wind.
  # Input:
  #   big: raster with wind data, large spatial extent
  #   small: raster with wind data, small spatial extent
  #   dims: dimension of the raster
  #   inds: list of indices, created with FVMInds
  #   bound.cond: boundary condition (periodic or endless world)
  # Output:
  #   Sparse matrix representing advection due to wind, as specified by
  #     wind velocities found in 'small'.

  big.uwind <- big$uwind
  big.vwind <- big$vwind

  small.uwind <- small$uwind
  small.vwind <- small$vwind

  u.faces <- interpWind(big.uwind, small.uwind, "u", with.bound = FALSE)
  v.faces <- interpWind(big.vwind, small.vwind, "v", with.bound = FALSE)
  velocity <- cbind(u.faces, v.faces)

  # wind vmatrix
  C <- windConvFVM(dims, velocity, inds, bc = bound.cond, bound = FALSE)
  return(C)
}


### B. time averaged, coupled SO2-SO4, with wind

ll.data.Wind <- function(theta, y, X, mats, int = FALSE){
  # Calculate the log-likelihood given sampled parameters for coupled
  #   SO2-SO4 time-averaged model (WITH WIND!).
  # Input:
  #   theta: sampled parameters
  #   y: observed data
  #   X: vector of SOZ emissions
  #   mats: matrix operators representing diffusion, advection, etc.
  # Output:
  #   Log-likelihood of the coupled SO2-SO4 model, with wind.

  # grab parameters
  gamma <- theta[1]
  xi <- theta[2]
  beta <- theta[3]
  s2 <- theta[4]
  delta <- theta[5]
  alpha <- theta[6]

  if (int){
    b0 <- theta[7]
  }

  # create matrices
  D <- mats$D
  C <- mats$C
  n.cells <- nrow(D)
  A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
  A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)


  if (int){
    y.star <- y - b0
  } else {
    y.star <- y
  }

  # constanst for log-likelihood calc.
  A1.y <- A.1 %*% y.star
  sources <- solve(A.2, X * beta)
  sources.cp <- t(sources) %*% sources
  log.det.A1 <- as.numeric(determinant(A.1)$modulus)

  # log-likelihood
  ll <- (-n.cells * log(s2) / 2) + log.det.A1 -
    1 / 2 / s2 * (t(A1.y) %*% A1.y  - 2 * t(A1.y) %*% sources + sources.cp)

  # return ll
  return(as.numeric(ll))
}

ymissingSample <- function(y, missing, X, mats, theta, int = FALSE){

  # grab parameters
  gamma <- theta[1]
  xi <- theta[2]
  beta <- theta[3]
  s2 <- theta[4]
  delta <- theta[5]
  alpha <- theta[6]

  if (int){
    b0 <- theta[7]
  }

  # create matrices
  D <- mats$D
  C <- mats$C
  n.cells <- nrow(D)
  A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
  A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)

  # joint mean and precision matrix
  sources <- solve(A.2, X * beta)

  if (int){
    mu <- b0 + solve(A.1, sources)
  } else {
    mu <- solve(A.1, sources)
  }
  Q <- (t(A.1) %*% A.1) / s2

  # conditional mean and precision matrix
  Q.cond <- Q[missing, missing]
  mu.cond <- mu[missing] - (solve(Q.cond, Q[missing, !missing]) %*% (y[!missing] - mu[!missing]))

  # sample y.missing
  y[missing] <- rnorm.Q(Q.cond, mu.cond, zero.constraint = FALSE, canon = FALSE)

  return(y)
}


coupledMCMC.Wind <- function(N, y, X, mats, params, priors, update, iter, n.y = 1000, burnin = 0, intercept = FALSE){
  # MCMC sampler for the coupled SO2-SO4 model (WITH WIND!).
  # Input:
  #   N: number of MCMC samples.
  #   y: observed data.
  #   X: vector of SO2 outputs.
  #   mats: matrix operators approximating physical processes.
  #   params: vector of initial parameters.
  #   priors: hyperparameters for prior distribution.
  #   update: number of iterations at which to update MCMC proposal dist
  #     (uses Shaby and Wells (2009) algorithm).
  #   iter: when to print out iterations.
  # Output:
  #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC
  # proposals.

  # number of observations
  n.cells <- length(y)

  # missing observations
  missing <- is.na(y)
  n.missing <- sum(missing)

  # initial parameter values
  gamma.k <- params$gamma
  xi.k <- params$xi
  beta.k <- params$beta
  s2.k <- params$s2
  delta.k <- params$delta
  alpha.k <- params$alpha

  if (intercept){
    # read in intercept
    b0.k <- params$b0
    # all parameters in a vector
    theta.k <- c(gamma.k, xi.k, beta.k, s2.k, delta.k, alpha.k, b0.k)
  } else {
    # all parameters in a vector
    theta.k <- c(gamma.k, xi.k, beta.k, s2.k, delta.k, alpha.k)
  }

  # sample y.0 (missing) | theta.0
  y.k <- ymissingSample(y, missing, X, mats, theta.k, int = intercept)

  # matrix to hold parameters
  samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))

  if (intercept) {
    colnames(samples) <- c("gamma", "xi", "beta", "s2", "delta", "alpha", "intercept")
  } else {
    colnames(samples) <- c("gamma", "xi", "beta", "s2", "delta", "alpha")
  }

  # matrix to hold imputed missing y data
  y.samples <- matrix(NA_real_, nrow = n.y, ncol = n.missing)
  y.s.idx <- sample((burnin + 1):N, size = n.y, replace = FALSE)
  y.s.idx <- sort(y.s.idx)

  # vector to hold correction factor samples for DIC calculation
  dic.vec <- rep(NA_real_, length = N)

  # calculate log-likelihood for current values
  ll.curr <- ll.data.Wind(theta.k, y.k, X, mats, int = intercept)

  ### create adaptive structures

  # matrix to hold samples
  update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))

  # default constants
  c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1

  # structure for gamma
  ad.gamma <- list()
  ad.gamma$s0 <- 1
  ad.gamma$var <- 10
  ad.gamma$n.jumps <- 0
  ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)

  # structure for xi
  ad.xi <- list()
  ad.xi$s0 <- 1
  ad.xi$var <- 1
  ad.xi$n.jumps <- 0
  ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)

  # structure for s2
  ad.s2 <- list()
  ad.s2$s0 <- 1
  ad.s2$var <- 10
  ad.s2$n.jumps <- 0
  ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)

  # structure for alpha
  ad.alpha <- list()
  ad.alpha$s0 <- 1
  ad.alpha$var <- 1
  ad.alpha$n.jumps <- 0
  ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)

  ### generate mcmc samples
  for (k in 1:N){

    # a. sample gamma
    gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)

    if (gamma.star > 0){

      # calculate new log-likelihood
      theta.star <- theta.k
      theta.star[1] <- gamma.star
      ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)

      # determine log-MH ratio
      mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
      mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
      mh.ratio <- mh.num - mh.denom

      if (log(runif(1)) < mh.ratio){
        # accept new gamma
        gamma.k <- gamma.star
        theta.k <- theta.star
        ll.curr <- ll.star
        ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
      }
    }

    # b. sample xi
    xi.star <- rnorm(n = 1, mean = xi.k, sd = ad.xi$sd)

    if (xi.star > 0){

      # calculate new log-likelihood
      theta.star <- theta.k
      theta.star[2] <- xi.star
      ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)

      # determine log-MH ratio
      mh.num <- ll.star + dexp(xi.star, rate = priors$xi, log = TRUE)
      mh.denom <- ll.curr + dexp(xi.k, rate = priors$xi, log = TRUE)
      mh.ratio <- mh.num - mh.denom

      if (log(runif(1)) < mh.ratio){
        # accept new xi
        xi.k <- xi.star
        theta.k <- theta.star
        ll.curr <- ll.star
        ad.xi$n.jumps <- ad.xi$n.jumps + 1
      }
    }

    # sample alpha
    alpha.star <- rnorm(n = 1, mean = alpha.k, sd = ad.alpha$sd)

    if (alpha.star > 0){
      # calculate new log-likelihood
      theta.star <- theta.k
      theta.star[6] <- alpha.star
      ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)

      # determine log-MH ratio
      mh.num <- ll.star + dnorm(alpha.star, mean = 0, sd = priors$alpha.sd, log = TRUE)
      mh.denom <- ll.curr + dnorm(alpha.k, mean = 0, sd = priors$alpha.sd, log = TRUE)
      mh.ratio <- mh.num - mh.denom

      if (log(runif(1)) < mh.ratio){
        # accept new alpha
        alpha.k <- alpha.star
        theta.k <- theta.star
        ll.curr <- ll.star
        ad.alpha$n.jumps <- ad.alpha$n.jumps + 1
      }

    }

    # c. sample s2
    s2.star <- rnorm(n = 1, mean = s2.k, sd = ad.s2$sd)

    if (s2.star > 0){

      # calculate new log-likelihood
      theta.star <- theta.k
      theta.star[4] <- s2.star
      ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)

      # determine log-MH ratio
      mh.num <- ll.star + dexp(s2.star, rate = priors$s2, log = TRUE)
      mh.denom <- ll.curr + dexp(s2.k, rate = priors$s2, log = TRUE)
      mh.ratio <- mh.num - mh.denom

      if (log(runif(1)) < mh.ratio){
        # accept new s2
        s2.k <- s2.star
        theta.k <- theta.star
        ll.curr <- ll.star
        ad.s2$n.jumps <- ad.s2$n.jumps + 1
      }
    }

    # d. sample beta (and b0, if applicable)

    if (intercept) {

      # truncated normal prior on beta
      D <- mats$D
      C <- mats$C
      n.cells <- nrow(D)
      A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
      A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
      A2inv.X <- solve(A.2, X)
      A1.y <- A.1 %*% (y.k - b0.k) # subtract intercept from y.k

      beta.var <- 1 / (as.numeric((t(A2inv.X) %*% A2inv.X / s2.k)) + (1 / priors$beta))
      beta.mean <- beta.var * as.numeric(t(A2inv.X) %*% A1.y) / s2.k

      # beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
      beta.k <- rtruncnorm(n = 1, a = 0, b = Inf, mean = beta.mean, sd = sqrt(beta.var))

    } else {

      # truncated normal prior on beta
      D <- mats$D
      C <- mats$C
      n.cells <- nrow(D)
      A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
      A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
      A2inv.X <- solve(A.2, X)
      A1.y <- A.1 %*% y.k

      beta.var <- 1 / (as.numeric((t(A2inv.X) %*% A2inv.X / s2.k)) + (1 / priors$beta))
      beta.mean <- beta.var * as.numeric(t(A2inv.X) %*% A1.y) / s2.k

      # beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
      beta.k <- rtruncnorm(n = 1, a = 0, b = Inf, mean = beta.mean, sd = sqrt(beta.var))
    }

    # save sample
    theta.k[3] <- beta.k

    # d. (v2) sample b0 (if applicable)
    if (intercept){

      # truncated normal prior on beta
      D <- mats$D
      C <- mats$C
      n.cells <- nrow(D)
      A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
      A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
      A2inv.X <- solve(A.2, X)
      A1.y <- A.1 %*% y.k
      A1.u <- A.1 %*% rep(1, n.cells)

      # variance
      var.b0 <- as.numeric(1 / (((t(A1.u) %*% A1.u) / s2.k) + (1 / priors$b0)))

      # mean
      mu.comp1 <- t(A1.y) %*% A1.u
      mu.comp2 <- beta.k * (t(A2inv.X) %*% A1.u)
      mu.b0 <- var.b0 * as.numeric((mu.comp1 - mu.comp2) / s2.k)

      # sample b0
      b0.k <- rtruncnorm(n = 1, a = 0, b = Inf, mean = mu.b0, sd = sqrt(var.b0))

      # save b0
      theta.k[7] <- b0.k
    }

    # calculate new likelihood
    ll.curr <- ll.data.Wind(theta.k, y.k, X, mats, int = intercept)

    # save ll for DIC calc
    dic.vec[k] <- -2 * ll.curr

    # e. save theta.k
    samples[k,] <- theta.k

    # f. sample y.k
    y.k <- ymissingSample(y.k, missing, X, mats, theta.k, int = intercept)

    # g. update adaptive structures

     # determine if y.k should be saved
    y.in.samp <- (y.s.idx == k)
    if (sum(y.in.samp) > 0){
      # save y.k[missing]
      y.samples[which(y.in.samp), ] <- y.k[missing]
    }

    if (k %% update > 0){
      # add theta.k to sample matrix
      update.samples[k %% update, ] <- theta.k
    } else {
      # add theta.k to sample matrix
      update.samples[update, ] <- theta.k

      # update constants
      gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1

      ### update gamma adaptive structure
      s.hat <- var(update.samples[,1])
      r.hat <- ad.gamma$n.jumps / update
      log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
      ad.gamma$var <- exp(log.var)
      ad.gamma$s0 <- ad.gamma$s0 +
        gamma1 * (s.hat - ad.gamma$s0)
      ad.gamma$n.jumps <- 0
      ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)

      ### update xi adaptive structure
      s.hat <- var(update.samples[,2])
      r.hat <- ad.xi$n.jumps / update
      log.var <- log(ad.xi$var) + gamma2 * (r.hat - r.opt)
      ad.xi$var <- exp(log.var)

      if (s.hat == 0 & gamma1 == 1){
        ad.xi$s0 <- ad.xi$s0 / 2
      } else {
        ad.xi$s0 <- ad.xi$s0 +
          gamma1 * (s.hat - ad.xi$s0)
      }
      ad.xi$n.jumps <- 0
      ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)

      ### update s2 adaptive structure
      s.hat <- var(update.samples[,4])
      r.hat <- ad.s2$n.jumps / update
      log.var <- log(ad.s2$var) + gamma2 * (r.hat - r.opt)
      ad.s2$var <- exp(log.var)
      ad.s2$s0 <- ad.s2$s0 +
        gamma1 * (s.hat - ad.s2$s0)
      ad.s2$n.jumps <- 0
      ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)

      ### update alpha adaptive structure
      s.hat <- var(update.samples[,6])
      r.hat <- ad.alpha$n.jumps / update
      log.var <- log(ad.alpha$var) + gamma2 * (r.hat - r.opt)
      ad.alpha$var <- exp(log.var)
      ad.alpha$s0 <- ad.alpha$s0 +
        gamma1 * (s.hat - ad.alpha$s0)
      ad.alpha$n.jumps <- 0
      ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)

      # update number of adaptive steps
      time <- time + 1
    }

    # f. print out the current iteration
    if (iter > 0){
      if (k %% iter == 0){
        print(paste("iteration:", k))
      }
    }
  }

  # save adaptive structures
  ad.str <- list(gamma = ad.gamma,
                 xi = ad.xi,
                 alpha = ad.alpha,
                 s2 = ad.s2)

  # return samples
  results <- list(samples = samples, y.missing = y.samples, dic = dic.vec, ad.str = ad.str)
  return(results)
}


simCoupled.Wind <- function(theta, mats, X){
  # Simulate coupled time-averaged model given samples from posterior
  #   (WITH WIND!).
  # Input:
  #   theta: model parameters
  #   mats: matrix structures used in model
  #   X: vector of SO2 emissions
  # Output:
  #   Simulated SO4 surface (as vector).

  # grab parameters
  gamma <- theta[1]
  xi <- theta[2]
  beta <- theta[3]
  s2 <- theta[4]
  delta <- theta[5]
  alpha <- theta[6]

  D <- mats$D
  C <- mats$C

  n.cells <- nrow(D)
  A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
  A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)


  sources <- solve(A.2, X * beta)
  mu <- solve(A.1, sources)

  Q <- t(A.1) %*% A.1 / s2

  sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)

  return(as.vector(sample))
}



















# 
# 
# ### A. time averaged, coupled SO2-SO4, no wind
# 
# ll.data <- function(theta, y, X, mats){
#   # Calculate the log-likelihood given sampled parameters for coupled 
#   #   SO2-SO4 time-averaged model (ASSUMES NO WIND!).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the coupled SO2-SO4 model, without wind.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   
#   # create matrices
#   D <- mats$D
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + Diagonal(n.cells, 1)
#   
#   # constanst for log-likelihood calc.
#   A1.y <- A.1 %*% y
#   sources <- solve(A.2, X * beta)
#   sources.cp <- t(sources) %*% sources
#   log.det.A1 <- as.numeric(determinant(A.1)$modulus)
#   
#   # log-likelihood
#   ll <- (-n.cells * log(s2) / 2) + log.det.A1 -
#     1 / 2 / s2 * (t(A1.y) %*% A1.y  - 2 * t(A1.y) %*% sources + sources.cp)
#   
#   # return ll
#   return(as.numeric(ll))
# }
# 
# coupledMCMC <- function(N, y, X, mats, params, priors, update, iter){
#   # MCMC sampler for the coupled SO2-SO4 time-averaged model (ASSUMES NO WIND!).
#   # Input:
#   #   N: number of MCMC samples.
#   #   y: observed data.
#   #   X: vector of SO2 outputs.
#   #   mats: matrix operators approximating physical processes.
#   #   params: vector of initial parameters.
#   #   priors: hyperparameters for prior distribution.
#   #   update: number of iterations at which to update MCMC proposal dist
#   #     (uses Shaby and Wells (2009) algorithm).
#   #   iter: when to print out iterations.
#   # Output:
#   #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC 
#   # proposals.
#     
#   n.cells <- length(y)
# 
#   # initial parameter values
#   gamma.k <- params$gamma
#   xi.k <- params$xi
#   beta.k <- params$beta
#   s2.k <- params$s2
#   delta.k <- params$delta
#   
#   # all parameters in a vector
#   theta.k <- c(gamma.k, xi.k, beta.k, s2.k, delta.k)
#   
#   # matrix to hold parameters
#   samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))
#   colnames(samples) <- c("gamma", "xi", "beta", "s2", "delta")
#   
#   # vector to hold correction factor samples for DIC calculation
#   dic.vec <- rep(NA_real_, length = N)
#   
#   # calculate log-likelihood for current values
#   ll.curr <- ll.data(theta.k, y, X, mats)
#   
#   ### create adaptive structures
#   
#   # matrix to hold samples
#   update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))
#   
#   # default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1
#   
#   # structure for gamma
#   ad.gamma <- list()
#   ad.gamma$s0 <- 1
#   ad.gamma$var <- 10
#   ad.gamma$n.jumps <- 0
#   ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#   
#   # structure for xi
#   ad.xi <- list()
#   ad.xi$s0 <- 1
#   ad.xi$var <- 1
#   ad.xi$n.jumps <- 0
#   ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#   
#   # structure for s2
#   ad.s2 <- list()
#   ad.s2$s0 <- 1
#   ad.s2$var <- 10
#   ad.s2$n.jumps <- 0
#   ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#   
#   ### generate mcmc samples
#   for (k in 1:N){
#     
#     # a. sample gamma
#     gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)
#     
#     if (gamma.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[1] <- gamma.star
#       ll.star <- ll.data(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new gamma
#         gamma.k <- gamma.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
#       }
#     }
#     
#     # b. sample xi
#     xi.star <- rnorm(n = 1, mean = xi.k, sd = ad.xi$sd)
#     
#     if (xi.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[2] <- xi.star
#       ll.star <- ll.data(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(xi.star, rate = priors$xi, log = TRUE)
#       mh.denom <- ll.curr + dexp(xi.k, rate = priors$xi, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new xi
#         xi.k <- xi.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.xi$n.jumps <- ad.xi$n.jumps + 1
#       }
#     }
#     
#     # c. sample s2
#     s2.star <- rnorm(n = 1, mean = s2.k, sd = ad.s2$sd)
#     
#     if (s2.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[4] <- s2.star
#       ll.star <- ll.data(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(s2.star, rate = priors$s2, log = TRUE)
#       mh.denom <- ll.curr + dexp(s2.k, rate = priors$s2, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new s2
#         s2.k <- s2.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.s2$n.jumps <- ad.s2$n.jumps + 1
#       }
#     }
#     
#     # d. sample beta
#     
#     D <- mats$D
#     n.cells <- nrow(D)
#     A.1 <- gamma.k * D + Diagonal(n.cells, delta.k)
#     A.2 <- (gamma.k / xi.k) * D + Diagonal(n.cells, 1)
#     A2inv.X <- solve(A.2, X)
#     A1.y <- A.1 %*% y
#     
#     beta.var <- 1 / (as.numeric((t(A2inv.X) %*% A2inv.X / s2.k)) + (1 / priors$beta))
#     beta.mean <- beta.var * as.numeric(t(A2inv.X) %*% A1.y) / s2.k
#     
#     beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
#     
#     # calculate new likelihood
#     theta.k[3] <- beta.k
#     ll.curr <- ll.data(theta.k, y, X, mats)
#     
#     # save ll for DIC calc
#     dic.vec[k] <- -2 * ll.curr
#     
#     # d. save theta.k
#     samples[k,] <- theta.k
#     
#     # e. update adaptive structures
#     
#     if (k %% update > 0){
#       # add theta.k to sample matrix
#       update.samples[k %% update, ] <- theta.k
#     } else {
#       # add theta.k to sample matrix
#       update.samples[update, ] <- theta.k
#       
#       # update constants
#       gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1
#       
#       ### update gamma adaptive structure
#       s.hat <- var(update.samples[,1])
#       r.hat <- ad.gamma$n.jumps / update
#       log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
#       ad.gamma$var <- exp(log.var)
#       ad.gamma$s0 <- ad.gamma$s0 +
#         gamma1 * (s.hat - ad.gamma$s0)
#       ad.gamma$n.jumps <- 0
#       ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#       
#       ### update xi adaptive structure
#       s.hat <- var(update.samples[,2])
#       r.hat <- ad.xi$n.jumps / update
#       log.var <- log(ad.xi$var) + gamma2 * (r.hat - r.opt)
#       ad.xi$var <- exp(log.var)
#       ad.xi$s0 <- ad.xi$s0 +
#         gamma1 * (s.hat - ad.xi$s0)
#       ad.xi$n.jumps <- 0
#       ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#       
#       ### update s2 adaptive structure
#       s.hat <- var(update.samples[,4])
#       r.hat <- ad.s2$n.jumps / update
#       log.var <- log(ad.s2$var) + gamma2 * (r.hat - r.opt)
#       ad.s2$var <- exp(log.var)
#       ad.s2$s0 <- ad.s2$s0 +
#         gamma1 * (s.hat - ad.s2$s0)
#       ad.s2$n.jumps <- 0
#       ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#       
#       # update number of adaptive steps
#       time <- time + 1
#     }
#     
#     
#     # f. print out the current iteration
#     if (iter > 0){
#       if (k %% iter == 0){
#         print(paste("iteration:", k))
#       }
#     }
#   }
#   
#   # save adaptive structures
#   ad.str <- list(gamma = ad.gamma,
#                  xi = ad.xi,
#                  s2 = ad.s2)
#   
#   # return samples
#   results <- list(samples = samples, dic = dic.vec, ad.str = ad.str)
#   
#   # return samples
#   return(results)
# }
# 
# simCoupled <- function(theta, mats, X){
#   # Simulate coupled time-averaged model given samples from posterior 
#   #   (NO WIND!).
#   # Input:
#   #   theta: model parameters
#   #   mats: matrix structures used in model
#   #   X: vector of SO2 emissions
#   # Output:
#   #   Simulated SO4 surface (as vector).
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   
#   # diffusion matrix
#   D <- mats$D
#   
#   # operator matrices
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + Diagonal(n.cells, 1)
#   
#   # SO2 surface
#   sources <- solve(A.2, X * beta)
#   # mean of SO4 model
#   mu <- solve(A.1, sources)
#   
#   # precision matrix
#   Q <- (A.1 %*% A.1) / s2
#   
#   # sample SO4 surface
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# ### B. time averaged, coupled SO2-SO4, with wind
# 
# ll.data.Wind <- function(theta, y, X, mats, int = FALSE){
#   # Calculate the log-likelihood given sampled parameters for coupled 
#   #   SO2-SO4 time-averaged model (WITH WIND!).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the coupled SO2-SO4 model, with wind.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
# 
#   if (int){
#     b0 <- theta[7]
#   }
#   
#   # create matrices
#   D <- mats$D
#   C <- mats$C
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
# 
# 
#   if (int){
#     y.star <- y - b0
#   } else {
#     y.star <- y
#   }
#   
#   # constanst for log-likelihood calc.
#   A1.y <- A.1 %*% y.star
#   sources <- solve(A.2, X * beta)
#   sources.cp <- t(sources) %*% sources
#   log.det.A1 <- as.numeric(determinant(A.1)$modulus)
#   
#   # log-likelihood
#   ll <- (-n.cells * log(s2) / 2) + log.det.A1 -
#     1 / 2 / s2 * (t(A1.y) %*% A1.y  - 2 * t(A1.y) %*% sources + sources.cp)
#   
#   # return ll
#   return(as.numeric(ll))
# }
# 
# ymissingSample <- function(y, missing, X, mats, theta, int = FALSE){
# 
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
# 
#   if (int){
#     b0 <- theta[7]
#   }
#   
#   # create matrices
#   D <- mats$D
#   C <- mats$C
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
# 
#   # joint mean and precision matrix
#   sources <- solve(A.2, X * beta)
# 
#   if (int){
#     mu <- b0 + solve(A.1, sources)
#   } else {
#     mu <- solve(A.1, sources)
#   }
#   Q <- (t(A.1) %*% A.1) / s2
# 
#   # conditional mean and precision matrix
#   Q.cond <- Q[missing, missing]
#   mu.cond <- mu[missing] - (solve(Q.cond, Q[missing, !missing]) %*% (y[!missing] - mu[!missing]))
# 
#   # sample y.missing
#   y[missing] <- rnorm.Q(Q.cond, mu.cond, zero.constraint = FALSE, canon = FALSE)
# 
#   return(y)
# }
# 
# coupledMCMC.Wind <- function(N, y, X, mats, params, priors, update, iter, n.y = 1000, burnin = 0, intercept = FALSE){
#   # MCMC sampler for the coupled SO2-SO4 model (WITH WIND!).
#   # Input:
#   #   N: number of MCMC samples.
#   #   y: observed data.
#   #   X: vector of SO2 outputs.
#   #   mats: matrix operators approximating physical processes.
#   #   params: vector of initial parameters.
#   #   priors: hyperparameters for prior distribution.
#   #   update: number of iterations at which to update MCMC proposal dist
#   #     (uses Shaby and Wells (2009) algorithm).
#   #   iter: when to print out iterations.
#   # Output:
#   #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC
#   # proposals.
# 
#   # number of observations
#   n.cells <- length(y)
#   
#   # missing observations
#   missing <- is.na(y)
#   n.missing <- sum(missing)
# 
#   # initial parameter values
#   gamma.k <- params$gamma
#   xi.k <- params$xi
#   beta.k <- params$beta
#   s2.k <- params$s2
#   delta.k <- params$delta
#   alpha.k <- params$alpha
# 
#   if (intercept){
#     # read in intercept
#     b0.k <- params$b0
#     # all parameters in a vector
#     theta.k <- c(gamma.k, xi.k, beta.k, s2.k, delta.k, alpha.k, b0.k)
#   } else {
#     # all parameters in a vector
#     theta.k <- c(gamma.k, xi.k, beta.k, s2.k, delta.k, alpha.k)
#   }
# 
#   # sample y.0 (missing) | theta.0
#   y.k <- ymissingSample(y, missing, X, mats, theta.k, int = intercept)
# 
#   # matrix to hold parameters
#   samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))
# 
#   if (intercept) {
#     colnames(samples) <- c("gamma", "xi", "beta", "s2", "delta", "alpha", "intercept")
#   } else {
#     colnames(samples) <- c("gamma", "xi", "beta", "s2", "delta", "alpha")
#   }
#   
#   # matrix to hold imputed missing y data
#   y.samples <- matrix(NA_real_, nrow = n.y, ncol = n.missing)
#   y.s.idx <- sample((burnin + 1):N, size = n.y, replace = FALSE)
#   y.s.idx <- sort(y.s.idx)
# 
#   # vector to hold correction factor samples for DIC calculation
#   dic.vec <- rep(NA_real_, length = N)
#   
#   # calculate log-likelihood for current values
#   ll.curr <- ll.data.Wind(theta.k, y.k, X, mats, int = intercept)
#   
#   ### create adaptive structures
#   
#   # matrix to hold samples
#   update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))
#   
#   # default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1
#   
#   # structure for gamma
#   ad.gamma <- list()
#   ad.gamma$s0 <- 1
#   ad.gamma$var <- 10
#   ad.gamma$n.jumps <- 0
#   ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#   
#   # structure for xi
#   ad.xi <- list()
#   ad.xi$s0 <- 1
#   ad.xi$var <- 1
#   ad.xi$n.jumps <- 0
#   ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#   
#   # structure for s2
#   ad.s2 <- list()
#   ad.s2$s0 <- 1
#   ad.s2$var <- 10
#   ad.s2$n.jumps <- 0
#   ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#   
#   # structure for alpha
#   ad.alpha <- list()
#   ad.alpha$s0 <- 1
#   ad.alpha$var <- 1
#   ad.alpha$n.jumps <- 0
#   ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)
#   
#   ### generate mcmc samples
#   for (k in 1:N){
#     
#     # a. sample gamma
#     gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)
#     
#     if (gamma.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[1] <- gamma.star
#       ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new gamma
#         gamma.k <- gamma.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
#       }
#     }
#     
#     # b. sample xi
#     xi.star <- rnorm(n = 1, mean = xi.k, sd = ad.xi$sd)
#     
#     if (xi.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[2] <- xi.star
#       ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(xi.star, rate = priors$xi, log = TRUE)
#       mh.denom <- ll.curr + dexp(xi.k, rate = priors$xi, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new xi
#         xi.k <- xi.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.xi$n.jumps <- ad.xi$n.jumps + 1
#       }
#     }
#     
#     # sample alpha
#     alpha.star <- rnorm(n = 1, mean = alpha.k, sd = ad.alpha$sd)
#     
#     if (alpha.star > 0){
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[6] <- alpha.star
#       ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(alpha.star, mean = 0, sd = priors$alpha.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(alpha.k, mean = 0, sd = priors$alpha.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new alpha
#         alpha.k <- alpha.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.alpha$n.jumps <- ad.alpha$n.jumps + 1
#       }
#       
#     }
#     
#     # c. sample s2
#     s2.star <- rnorm(n = 1, mean = s2.k, sd = ad.s2$sd)
#     
#     if (s2.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[4] <- s2.star
#       ll.star <- ll.data.Wind(theta.star, y.k, X, mats, int = intercept)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(s2.star, rate = priors$s2, log = TRUE)
#       mh.denom <- ll.curr + dexp(s2.k, rate = priors$s2, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new s2
#         s2.k <- s2.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.s2$n.jumps <- ad.s2$n.jumps + 1
#       }
#     }
#     
#     # d. sample beta (and b0, if applicable)
# 
#     if (intercept) {
# 
#       # truncated normal prior on beta
#       D <- mats$D
#       C <- mats$C
#       n.cells <- nrow(D)
#       A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
#       A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
#       A2inv.X <- solve(A.2, X)
#       A1.y <- A.1 %*% (y.k - b0.k) # subtract intercept from y.k
#       
#       beta.var <- 1 / (as.numeric((t(A2inv.X) %*% A2inv.X / s2.k)) + (1 / priors$beta))
#       beta.mean <- beta.var * as.numeric(t(A2inv.X) %*% A1.y) / s2.k
# 
#       # beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
#       beta.k <- rtruncnorm(n = 1, a = 0, b = Inf, mean = beta.mean, sd = sqrt(beta.var))
# 
#     } else {
# 
#       # truncated normal prior on beta
#       D <- mats$D
#       C <- mats$C
#       n.cells <- nrow(D)
#       A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
#       A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
#       A2inv.X <- solve(A.2, X)
#       A1.y <- A.1 %*% y.k
# 
#       beta.var <- 1 / (as.numeric((t(A2inv.X) %*% A2inv.X / s2.k)) + (1 / priors$beta))
#       beta.mean <- beta.var * as.numeric(t(A2inv.X) %*% A1.y) / s2.k
# 
#       # beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
#       beta.k <- rtruncnorm(n = 1, a = 0, b = Inf, mean = beta.mean, sd = sqrt(beta.var))
#     }
#     
#     # save sample
#     theta.k[3] <- beta.k
# 
#     # d. (v2) sample b0 (if applicable)
#     if (intercept){
# 
#       # truncated normal prior on beta
#       D <- mats$D
#       C <- mats$C
#       n.cells <- nrow(D)
#       A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
#       A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
#       A2inv.X <- solve(A.2, X)
#       A1.y <- A.1 %*% y.k
#       A1.u <- A.1 %*% rep(1, n.cells)
# 
#       # variance
#       var.b0 <- as.numeric(1 / (((t(A1.u) %*% A1.u) / s2.k) + (1 / priors$b0)))
# 
#       # mean
#       mu.comp1 <- t(A1.y) %*% A1.u
#       mu.comp2 <- beta.k * (t(A2inv.X) %*% A1.u)
#       mu.b0 <- var.b0 * as.numeric((mu.comp1 - mu.comp2) / s2.k)
# 
#       # sample b0
#       b0.k <- rtruncnorm(n = 1, a = 0, b = Inf, mean = mu.b0, sd = sqrt(var.b0))
# 
#       # save b0
#       theta.k[7] <- b0.k
#     }
# 
#     # calculate new likelihood
#     ll.curr <- ll.data.Wind(theta.k, y.k, X, mats, int = intercept)
#     
#     # save ll for DIC calc
#     dic.vec[k] <- -2 * ll.curr
#     
#     # e. save theta.k
#     samples[k,] <- theta.k
# 
#     # f. sample y.k
#     y.k <- ymissingSample(y.k, missing, X, mats, theta.k, int = intercept)
#     
#     # g. update adaptive structures
# 
#      # determine if y.k should be saved
#     y.in.samp <- (y.s.idx == k)
#     if (sum(y.in.samp) > 0){
#       # save y.k[missing]
#       y.samples[which(y.in.samp), ] <- y.k[missing]
#     }
#       
#     if (k %% update > 0){
#       # add theta.k to sample matrix
#       update.samples[k %% update, ] <- theta.k
#     } else {
#       # add theta.k to sample matrix
#       update.samples[update, ] <- theta.k
# 
#       # update constants
#       gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1
#       
#       ### update gamma adaptive structure
#       s.hat <- var(update.samples[,1])
#       r.hat <- ad.gamma$n.jumps / update
#       log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
#       ad.gamma$var <- exp(log.var)
#       ad.gamma$s0 <- ad.gamma$s0 +
#         gamma1 * (s.hat - ad.gamma$s0)
#       ad.gamma$n.jumps <- 0
#       ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#       
#       ### update xi adaptive structure
#       s.hat <- var(update.samples[,2])
#       r.hat <- ad.xi$n.jumps / update
#       log.var <- log(ad.xi$var) + gamma2 * (r.hat - r.opt)
#       ad.xi$var <- exp(log.var)
# 
#       if (s.hat == 0 & gamma1 == 1){
#         ad.xi$s0 <- ad.xi$s0 / 2
#       } else {
#         ad.xi$s0 <- ad.xi$s0 +
#           gamma1 * (s.hat - ad.xi$s0)
#       }
#       ad.xi$n.jumps <- 0
#       ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
# 
#       ### update s2 adaptive structure
#       s.hat <- var(update.samples[,4])
#       r.hat <- ad.s2$n.jumps / update
#       log.var <- log(ad.s2$var) + gamma2 * (r.hat - r.opt)
#       ad.s2$var <- exp(log.var)
#       ad.s2$s0 <- ad.s2$s0 +
#         gamma1 * (s.hat - ad.s2$s0)
#       ad.s2$n.jumps <- 0
#       ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#       
#       ### update alpha adaptive structure
#       s.hat <- var(update.samples[,6])
#       r.hat <- ad.alpha$n.jumps / update
#       log.var <- log(ad.alpha$var) + gamma2 * (r.hat - r.opt)
#       ad.alpha$var <- exp(log.var)
#       ad.alpha$s0 <- ad.alpha$s0 +
#         gamma1 * (s.hat - ad.alpha$s0)
#       ad.alpha$n.jumps <- 0
#       ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)
#       
#       # update number of adaptive steps
#       time <- time + 1
#     }
#     
#     # f. print out the current iteration
#     if (iter > 0){
#       if (k %% iter == 0){
#         print(paste("iteration:", k))
#       }
#     }
#   }
#   
#   # save adaptive structures
#   ad.str <- list(gamma = ad.gamma,
#                  xi = ad.xi,
#                  alpha = ad.alpha,
#                  s2 = ad.s2)
#   
#   # return samples
#   results <- list(samples = samples, y.missing = y.samples, dic = dic.vec, ad.str = ad.str)
#   return(results)
# }
# 
# simCoupled.Wind <- function(theta, mats, X){
#   # Simulate coupled time-averaged model given samples from posterior 
#   #   (WITH WIND!).
#   # Input:
#   #   theta: model parameters
#   #   mats: matrix structures used in model
#   #   X: vector of SO2 emissions
#   # Output:
#   #   Simulated SO4 surface (as vector).
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
#   
#   D <- mats$D
#   C <- mats$C
#   
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#   
#   
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
#   
#   Q <- t(A.1) %*% A.1 / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# ### C. snapshot model, coupled SO2-SO4, no wind
# 
# ll.data.snap <- function(theta, y, X, mats){
#   # Calculate the log-likelihood given sampled parameters for coupled 
#   #   SO2-SO4 snapshot model (no wind).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the snapshot, coupled SO2-SO4 model, no wind.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   
#   # create matrices
#   D <- mats$D
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + Diagonal(n.cells, 1)
#   
#   # constanst for log-likelihood calc.
#   A1.y <- A.1 %*% y
#   sources <- solve(A.2, X * beta)
#   sources.cp <- t(sources) %*% solve(A.1, sources)
#   log.det.A1 <- as.numeric(determinant(A.1)$modulus)
#   
#   # log-likelihood
#   ll <- (n.cells * log(2) / 2) - (n.cells * log(s2) / 2) + (log.det.A1 / 2) -
#     1 / s2 * (t(y) %*% A1.y  - 2 * t(y) %*% sources + sources.cp)
#   
#   # return ll
#   return(as.numeric(ll))
# }
# 
# coupledMCMC.snap <- function(N, y, X, mats, params, priors, update, iter){
#   # MCMC sampler for the coupled SO2-SO4 snapshot model (no wind).
#   # Input:
#   #   N: number of MCMC samples.
#   #   y: observed data.
#   #   X: vector of SO2 outputs.
#   #   mats: matrix operators approximating physical processes.
#   #   params: vector of initial parameters.
#   #   priors: hyperparameters for prior distribution.
#   #   update: number of iterations at which to update MCMC proposal dist
#   #     (uses Shaby and Wells (2009) algorithm).
#   #   iter: when to print out iterations.
#   # Output:
#   #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC 
#   # proposals.
#   
#   n.cells <- length(y)
#   
#   # initial parameter values
#   gamma.k <- params$gamma
#   xi.k <- params$xi
#   beta.k <- params$beta
#   s2.k <- params$s2
#   delta.k <- params$delta
#   
#   # all parameters in a vector
#   theta.k <- c(gamma.k, xi.k, beta.k, s2.k, delta.k)
#   
#   # matrix to hold parameters
#   samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))
#   colnames(samples) <- c("gamma", "xi", "beta", "s2", "delta")
#   
#   # vector to hold correction factor samples for DIC calculation
#   dic.vec <- rep(NA_real_, length = N)
#   
#   # calculate log-likelihood for current values
#   ll.curr <- ll.data.snap(theta.k, y, X, mats)
#   
#   ### create adaptive structures
#   
#   # matrix to hold samples
#   update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))
#   
#   # default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1
#   
#   # structure for gamma
#   ad.gamma <- list()
#   ad.gamma$s0 <- 1
#   ad.gamma$var <- 10
#   ad.gamma$n.jumps <- 0
#   ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#   
#   # structure for xi
#   ad.xi <- list()
#   ad.xi$s0 <- 1
#   ad.xi$var <- 1
#   ad.xi$n.jumps <- 0
#   ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#   
#   # structure for s2
#   ad.s2 <- list()
#   ad.s2$s0 <- 1
#   ad.s2$var <- 1
#   ad.s2$n.jumps <- 0
#   ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#   
#   ### generate mcmc samples
#   for (k in 1:N){
#     
#     # a. sample gamma
#     gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)
#     
#     if (gamma.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[1] <- gamma.star
#       ll.star <- ll.data.snap(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new gamma
#         gamma.k <- gamma.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
#       }
#     }
#     
#     # b. sample xi
#     xi.star <- rnorm(n = 1, mean = xi.k, sd = ad.xi$sd)
#     
#     if (xi.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[2] <- xi.star
#       ll.star <- ll.data.snap(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(xi.star, rate = priors$xi, log = TRUE)
#       mh.denom <- ll.curr + dexp(xi.k, rate = priors$xi, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new xi
#         xi.k <- xi.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.xi$n.jumps <- ad.xi$n.jumps + 1
#       }
#     }
#     
#     # c. sample s2
#     s2.star <- rnorm(n = 1, mean = s2.k, sd = ad.s2$sd)
#     
#     if (s2.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[4] <- s2.star
#       ll.star <- ll.data.snap(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(s2.star, rate = priors$s2, log = TRUE)
#       mh.denom <- ll.curr + dexp(s2.k, rate = priors$s2, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new s2
#         s2.k <- s2.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.s2$n.jumps <- ad.s2$n.jumps + 1
#       }
#     }
#     
#     # d. sample beta
#     
#     D <- mats$D
#     n.cells <- nrow(D)
#     A.1 <- gamma.k * D + Diagonal(n.cells, delta.k)
#     A.2 <- (gamma.k / xi.k) * D + Diagonal(n.cells, 1)
#     A2inv.X <- solve(A.2, X)
#     
#     beta.var <- 1 / (2 * as.numeric((t(A2inv.X) %*% solve(A.1, A2inv.X) / s2.k)) +
#                        (1 / priors$beta))
#     beta.mean <- beta.var * 2 * as.numeric(t(A2inv.X) %*% y) / s2.k
#     
#     beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
#     
#     # calculate new likelihood
#     theta.k[3] <- beta.k
#     ll.curr <- ll.data.snap(theta.k, y, X, mats)
#     
#     # save ll for DIC calc
#     dic.vec[k] <- -2 * ll.curr
#     
#     # d. save theta.k
#     samples[k,] <- theta.k
#     
#     # e. update adaptive structures
#     
#     if (k %% update > 0){
#       # add theta.k to sample matrix
#       update.samples[k %% update, ] <- theta.k
#     } else {
#       # add theta.k to sample matrix
#       update.samples[update, ] <- theta.k
#       
#       # update constants
#       gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1
#       
#       ### update gamma adaptive structure
#       s.hat <- var(update.samples[,1]) + 0.00001
#       r.hat <- ad.gamma$n.jumps / update
#       log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
#       ad.gamma$var <- exp(log.var)
#       ad.gamma$s0 <- ad.gamma$s0 +
#         gamma1 * (s.hat - ad.gamma$s0)
#       ad.gamma$n.jumps <- 0
#       ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#       
#       ### update xi adaptive structure
#       s.hat <- var(update.samples[,2]) + 0.00001
#       r.hat <- ad.xi$n.jumps / update
#       log.var <- log(ad.xi$var) + gamma2 * (r.hat - r.opt)
#       ad.xi$var <- exp(log.var)
#       ad.xi$s0 <- ad.xi$s0 +
#         gamma1 * (s.hat - ad.xi$s0)
#       ad.xi$n.jumps <- 0
#       ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#       
#       ### update s2 adaptive structure
#       s.hat <- var(update.samples[,4]) + 0.00001
#       r.hat <- ad.s2$n.jumps / update
#       log.var <- log(ad.s2$var) + gamma2 * (r.hat - r.opt)
#       ad.s2$var <- exp(log.var)
#       ad.s2$s0 <- ad.s2$s0 +
#         gamma1 * (s.hat - ad.s2$s0)
#       ad.s2$n.jumps <- 0
#       ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#       
#       # update number of adaptive steps
#       time <- time + 1
#     }
#     
#     # f. print out the current iteration
#     if (iter > 0){
#       if (k %% iter == 0){
#         print(paste("iteration:", k))
#       }
#     }
#   }
#   
#   # save adaptive structures
#   ad.str <- list(gamma = ad.gamma,
#                  xi = ad.xi,
#                  s2 = ad.s2)
#   
#   # return samples
#   results <- list(samples = samples, dic = dic.vec, ad.str = ad.str)
#   
#   # return samples
#   return(results)
# }
# 
# simCoupled.snap <- function(theta, mats, X){
#   # Simulate coupled snapshot model given samples from posterior (no wind).
#   # Input:
#   #   theta: model parameters
#   #   mats: matrix structures used in model
#   #   X: vector of SO2 emissions
#   # Output:
#   #   Simulated SO4 surface (as vector).
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   
#   D <- mats$D
#   
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + Diagonal(n.cells, 1)
#   
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
#   
#   Q <- 2 * A.1 / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# ### E. time averaged, uncoupled SO4 model, with wind
# 
# ll.data.uncoupled.W <- function(theta, y, X, mats){
#   # Calculate the log-likelihood given sampled parameters for uncoupled 
#   #   SO2-SO4 time-averaged model (WITH WIND!).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the time-averaged, uncoupled SO2-SO4 model, with wind.
#   
#   # grab parameters
#   gamma <- theta[1]
#   beta <- theta[2]
#   s2 <- theta[3]
#   delta <- theta[4]
#   alpha <- theta[5]
#   
#   # create matrices
#   D <- mats$D
#   C <- mats$C
#   n.cells <- nrow(D)
#   A <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   
#   # constanst for log-likelihood calc.
#   A.y <- A %*% y
#   X.b <- X * beta
#   
#   like.1 <- A.y - X.b
#   cp.like <- t(like.1) %*% like.1
#   log.det.A <- as.numeric(determinant(A)$modulus)
#   
#   # log-likelihood
#   ll <- (-n.cells * log(s2) / 2) + log.det.A -
#     1 / 2 / s2 * (cp.like)
#   
#   # return ll
#   return(as.numeric(ll))
# }
# 
# uncoupledMCMC.W <- function(N, y, X, mats, params, priors, update, iter){
#   # MCMC sampler for the uncoupled SO2-SO4 time-averaged model (with wind).
#   # Input:
#   #   N: number of MCMC samples.
#   #   y: observed data.
#   #   X: vector of SO2 outputs.
#   #   mats: matrix operators approximating physical processes.
#   #   params: vector of initial parameters.
#   #   priors: hyperparameters for prior distribution.
#   #   update: number of iterations at which to update MCMC proposal dist
#   #     (uses Shaby and Wells (2009) algorithm).
#   #   iter: when to print out iterations.
#   # Output:
#   #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC 
#   # proposals.
#   
#   n.cells <- length(y)
#   
#   # initial parameter values
#   gamma.k <- params$gamma
#   beta.k <- params$beta
#   s2.k <- params$s2
#   delta.k <- params$delta
#   alpha.k <- params$alpha
#   
#   # all parameters in a vector
#   theta.k <- c(gamma.k, beta.k, s2.k, delta.k, alpha.k)
#   
#   # matrix to hold parameters
#   samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))
#   colnames(samples) <- c("gamma", "beta", "s2", "delta", "alpha")
#   
#   # vector to hold correction factor samples for DIC calculation
#   dic.vec <- rep(NA_real_, length = N)
#   
#   # calculate log-likelihood for current values
#   ll.curr <- ll.data.uncoupled.W(theta.k, y, X, mats)
#   
#   ### create adaptive structures
#   
#   # matrix to hold samples
#   update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))
#   
#   # default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1
#   
#   # structure for gamma
#   ad.gamma <- list()
#   ad.gamma$s0 <- 1
#   ad.gamma$var <- 10
#   ad.gamma$n.jumps <- 0
#   ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#   
#   # structure for alpha
#   ad.alpha <- list()
#   ad.alpha$s0 <- 1
#   ad.alpha$var <- 1
#   ad.alpha$n.jumps <- 0
#   ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)
#   
#   # structure for s2
#   ad.s2 <- list()
#   ad.s2$s0 <- 1
#   ad.s2$var <- 10
#   ad.s2$n.jumps <- 0
#   ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#   
#   ### generate mcmc samples
#   for (k in 1:N){
#     
#     # a. sample gamma
#     gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)
#     
#     if (gamma.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[1] <- gamma.star
#       ll.star <- ll.data.uncoupled.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new gamma
#         gamma.k <- gamma.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
#       }
#     }
#     
#     # b. sample alpha
#     alpha.star <- rnorm(n = 1, mean = alpha.k, sd = ad.alpha$sd)
#     
#     if (alpha.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[5] <- alpha.star
#       ll.star <- ll.data.uncoupled.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(alpha.star, mean = 0, sd = priors$alpha.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(alpha.k, mean = 0, sd = priors$alpha.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new alpha
#         alpha.k <- alpha.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.alpha$n.jumps <- ad.alpha$n.jumps + 1
#       }
#     }
#     
#     # c. sample s2
#     s2.star <- rnorm(n = 1, mean = s2.k, sd = ad.s2$sd)
#     
#     if (s2.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[3] <- s2.star
#       ll.star <- ll.data.uncoupled.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(s2.star, rate = priors$s2, log = TRUE)
#       mh.denom <- ll.curr + dexp(s2.k, rate = priors$s2, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new s2
#         s2.k <- s2.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.s2$n.jumps <- ad.s2$n.jumps + 1
#       }
#     }
#     
#     # d. sample beta
#     
#     D <- mats$D
#     C <- mats$C
#     n.cells <- nrow(D)
#     A <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
#     A.y <- A %*% y
#     
#     beta.var <- 1 / (as.numeric((t(X) %*% X / s2.k)) + (1 / priors$beta))
#     beta.mean <- beta.var * as.numeric(t(X) %*% A.y) / s2.k
#     
#     beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
#     
#     # calculate new likelihood
#     theta.k[2] <- beta.k
#     ll.curr <- ll.data.uncoupled.W(theta.k, y, X, mats)
#     
#     # save ll for DIC calc
#     dic.vec[k] <- -2 * ll.curr
#     
#     # d. save theta.k
#     samples[k,] <- theta.k
#     
#     # e. update adaptive structures
#     
#     if (k %% update > 0){
#       # add theta.k to sample matrix
#       update.samples[k %% update, ] <- theta.k
#     } else {
#       # add theta.k to sample matrix
#       update.samples[update, ] <- theta.k
#       
#       # update constants
#       gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1
#       
#       ### update gamma adaptive structure
#       s.hat <- var(update.samples[,1]) + 0.000001
#       r.hat <- ad.gamma$n.jumps / update
#       log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
#       ad.gamma$var <- exp(log.var)
#       ad.gamma$s0 <- ad.gamma$s0 +
#         gamma1 * (s.hat - ad.gamma$s0)
#       ad.gamma$n.jumps <- 0
#       ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#       
#       ### update alpha adaptive structure
#       s.hat <- var(update.samples[,5]) + 0.000001
#       r.hat <- ad.alpha$n.jumps / update
#       log.var <- log(ad.alpha$var) + gamma2 * (r.hat - r.opt)
#       ad.alpha$var <- exp(log.var)
#       ad.alpha$s0 <- ad.alpha$s0 +
#         gamma1 * (s.hat - ad.alpha$s0)
#       ad.alpha$n.jumps <- 0
#       ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)
#       
#       ### update s2 adaptive structure
#       s.hat <- var(update.samples[,3]) + 0.000001
#       r.hat <- ad.s2$n.jumps / update
#       log.var <- log(ad.s2$var) + gamma2 * (r.hat - r.opt)
#       ad.s2$var <- exp(log.var)
#       ad.s2$s0 <- ad.s2$s0 +
#         gamma1 * (s.hat - ad.s2$s0)
#       ad.s2$n.jumps <- 0
#       ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#       
#       # update number of adaptive steps
#       time <- time + 1
#     }
#     
#     # f. print out the current iteration
#     if (iter > 0){
#       if (k %% iter == 0){
#         print(paste("iteration:", k))
#       }
#     }
#   }
#   
#   # save adaptive structures
#   ad.str <- list(gamma = ad.gamma,
#                  alpha = ad.alpha,
#                  s2 = ad.s2)
#   
#   # return samples
#   results <- list(samples = samples, dic = dic.vec, ad.str = ad.str)
#   
#   # return samples
#   return(results)
# }
# 
# simUncoupled.W <- function(theta, mats, X){
#   # Simulate uncoupled time-averaged model given samples from posterior 
#   #   (with wind).
#   # Input:
#   #   theta: model parameters
#   #   mats: matrix structures used in model
#   #   X: vector of SO2 emissions
#   # Output:
#   #   Simulated SO4 surface (as vector).
#   
#   # grab parameters
#   gamma <- theta[1]
#   beta <- theta[2]
#   s2 <- theta[3]
#   delta <- theta[4]
#   alpha <- theta[5]
#   
#   D <- mats$D
#   C <- mats$C
#   
#   n.cells <- nrow(D)
#   A <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   
#   mu <- solve(A, X * beta)
#   
#   Q <- (t(A) %*% A) / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# ### F. snapshot, uncoupled SO4 model, no wind 
# 
# ll.data.uncoupled.snap <- function(theta, y, X, mats){
#   # Calculate the log-likelihood given sampled parameters for uncoupled 
#   #   SO2-SO4 snapshot model (WITH WIND!).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the snaphsot, uncoupled SO2-SO4 model, with wind.
#   
#   # grab parameters
#   gamma <- theta[1]
#   beta <- theta[2]
#   s2 <- theta[3]
#   delta <- theta[4]
#   
#   # create matrices
#   D <- mats$D
#   n.cells <- nrow(D)
#   A <- gamma * D +  Diagonal(n.cells, delta)
#   
#   # constanst for log-likelihood calc.
#   A.y <- A %*% y
#   X.b <- X * beta
#   log.det.A <- as.numeric(determinant(A)$modulus)
#   
#   # log-likelihood
#   ll <- (n.cells * log(2) / 2) - (n.cells * log(s2) / 2) + (log.det.A / 2) -
#     1 / s2 * (t(y) %*% A.y - 2 * t(y) %*% X.b + t(X.b) %*% solve(A, X.b))
#   
#   # return ll
#   return(as.numeric(ll))
# }
# 
# uncoupledMCMC.snap <- function(N, y, X, mats, params, priors, update, iter){
#   # MCMC sampler for the uncoupled SO2-SO4 snapshot model (no wind).
#   # Input:
#   #   N: number of MCMC samples.
#   #   y: observed data.
#   #   X: vector of SO2 outputs.
#   #   mats: matrix operators approximating physical processes.
#   #   params: vector of initial parameters.
#   #   priors: hyperparameters for prior distribution.
#   #   update: number of iterations at which to update MCMC proposal dist
#   #     (uses Shaby and Wells (2009) algorithm).
#   #   iter: when to print out iterations.
#   # Output:
#   #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC 
#   # proposals.
#   
#   n.cells <- length(y)
#   
#   # initial parameter values
#   gamma.k <- params$gamma
#   beta.k <- params$beta
#   s2.k <- params$s2
#   delta.k <- params$delta
#   
#   # all parameters in a vector
#   theta.k <- c(gamma.k, beta.k, s2.k, delta.k)
#   
#   # matrix to hold parameters
#   samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))
#   colnames(samples) <- c("gamma", "beta", "s2", "delta")
#   
#   # vector to hold correction factor samples for DIC calculation
#   dic.vec <- rep(NA_real_, length = N)
#   
#   # calculate log-likelihood for current values
#   ll.curr <- ll.data.uncoupled.snap(theta.k, y, X, mats)
#   
#   ### create adaptive structures
#   
#   # matrix to hold samples
#   update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))
#   
#   # default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1
#   
#   # structure for gamma
#   ad.gamma <- list()
#   ad.gamma$s0 <- 1
#   ad.gamma$var <- 10
#   ad.gamma$n.jumps <- 0
#   ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#   
#   # structure for s2
#   ad.s2 <- list()
#   ad.s2$s0 <- 1
#   ad.s2$var <- 10
#   ad.s2$n.jumps <- 0
#   ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#   
#   ### generate mcmc samples
#   for (k in 1:N){
#     
#     # a. sample gamma
#     gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)
#     
#     if (gamma.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[1] <- gamma.star
#       ll.star <- ll.data.uncoupled.snap(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new gamma
#         gamma.k <- gamma.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
#       }
#     }
#     
#     # c. sample s2
#     s2.star <- rnorm(n = 1, mean = s2.k, sd = ad.s2$sd)
#     
#     if (s2.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[3] <- s2.star
#       ll.star <- ll.data.uncoupled.snap(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(s2.star, rate = priors$s2, log = TRUE)
#       mh.denom <- ll.curr + dexp(s2.k, rate = priors$s2, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new s2
#         s2.k <- s2.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.s2$n.jumps <- ad.s2$n.jumps + 1
#       }
#     }
#     
#     # d. sample beta
#     
#     D <- mats$D
#     n.cells <- nrow(D)
#     A <- gamma.k * D + Diagonal(n.cells, delta.k)
#     A.y <- A %*% y
#     
#     beta.var <- 1 / (2 * as.numeric((t(X) %*% solve(A, X)) / s2.k) + (1 / priors$beta))
#     beta.mean <- beta.var * 2 * as.numeric(t(X) %*% y) / s2.k
#     
#     beta.k <- rnorm(n = 1, mean = beta.mean, sd = sqrt(beta.var))
#     
#     # calculate new likelihood
#     theta.k[2] <- beta.k
#     ll.curr <- ll.data.uncoupled.snap(theta.k, y, X, mats)
#     
#     # save ll for DIC calc
#     dic.vec[k] <- -2 * ll.curr
#     
#     # d. save theta.k
#     samples[k,] <- theta.k
#     
#     # e. update adaptive structures
#     
#     if (k %% update > 0){
#       # add theta.k to sample matrix
#       update.samples[k %% update, ] <- theta.k
#     } else {
#       # add theta.k to sample matrix
#       update.samples[update, ] <- theta.k
#       
#       # update constants
#       gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1
#       
#       ### update gamma adaptive structure
#       s.hat <- var(update.samples[,1]) + 0.000001
#       r.hat <- ad.gamma$n.jumps / update
#       log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
#       ad.gamma$var <- exp(log.var)
#       ad.gamma$s0 <- ad.gamma$s0 +
#         gamma1 * (s.hat - ad.gamma$s0)
#       ad.gamma$n.jumps <- 0
#       ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#       
#       ### update s2 adaptive structure
#       s.hat <- var(update.samples[,3]) + 0.000001
#       r.hat <- ad.s2$n.jumps / update
#       log.var <- log(ad.s2$var) + gamma2 * (r.hat - r.opt)
#       ad.s2$var <- exp(log.var)
#       ad.s2$s0 <- ad.s2$s0 +
#         gamma1 * (s.hat - ad.s2$s0)
#       ad.s2$n.jumps <- 0
#       ad.s2$sd <- sqrt(ad.s2$var * ad.s2$s0)
#       
#       # update number of adaptive steps
#       time <- time + 1
#     }
#     
#     # f. print out the current iteration
#     if (iter > 0){
#       if (k %% iter == 0){
#         print(paste("iteration:", k))
#       }
#     }
#   }
#   
#   # save adaptive structures
#   ad.str <- list(gamma = ad.gamma,
#                  s2 = ad.s2)
#   
#   # return samples
#   results <- list(samples = samples, dic = dic.vec, ad.str = ad.str)
#   
#   # return samples
#   return(results)
# }
# 
# simUncoupled.snap <- function(theta, mats, X){
#   # Simulate uncoupled snapshot model given samples from posterior (no wind).
#   # Input:
#   #   theta: model parameters
#   #   mats: matrix structures used in model
#   #   X: vector of SO2 emissions
#   # Output:
#   #   Simulated SO4 surface (as vector).
#   
#   # grab parameters
#   gamma <- theta[1]
#   beta <- theta[2]
#   s2 <- theta[3]
#   delta <- theta[4]
#   
#   D <- mats$D
#   
#   n.cells <- nrow(D)
#   A <- gamma * D + Diagonal(n.cells, delta)
#   
#   mu <- solve(A, X * beta)
#   
#   Q <- 2 * A / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# ### G. Finite Volume Method Functions
# 
# interpWind <- function(r.old, r.new, direction, with.bound = TRUE){
#   # Linear interpolation of wind velocity across surface boundaries.
#   # Input:
#   #   r.old: wind raster with larger spatial extent.
#   #   r.new: wind raster with smaller spatial extent.
#   #   direction: direction of wind (either "u" or "v").
#   #   with.bound: boolean determining if velocities should be saved on
#   #     the boundary.
#   # Output:
#   #   Matrix with interpolated matrix velocities.
# 
#   # grab xy coordinates from rasters
#   xy.old <- coordinates(r.old)
#   xy.new <- coordinates(r.new)
# 
#   # dimensions of rasters (n.rows by n.cols)
#   dim.old <- dim(r.old)[1:2]
#   dim.new <- dim(r.new)[1:2]
# 
#   # cell indices for trimmed case
#   cell.inds <- cellFromXY(r.old, xy.new)
# 
#   # indices for the top boundary
#   top.bound <- (cell.inds[1] - dim.old[2]) + (-1):(dim.new[2])
# 
#   # cell indices with boundaries
#   all.bounds <- top.bound
# 
#   for (k in 1:(dim.new[1] + 1)){
#     next.bound <- top.bound + (k * dim.old[2])
#     all.bounds <- c(all.bounds, next.bound)
#   }
# 
#   # wind values for all points (including boundaries)
#   orig.vals <- values(r.old)[all.bounds]
#   # matrix of face value velocities
#   face.vals <- matrix(NA_real_, nrow = length(all.bounds), ncol = 2)
# 
#   # new dimensions
#   n.rows <- dim.new[1] + 2
#   n.cols <- dim.new[2] + 2
# 
#   if (direction == "u"){
#     ### west to east
# 
#     # name columns appropriately
#     colnames(face.vals) <- c("west", "east")
# 
#     for(i in 1:n.rows){
#       for(j in 1:n.cols){
# 
#         # current indice
#         cur.ind <- (i - 1) * n.cols + j
# 
#         if (j == 1){
#           # no west face velocity
#           west <- 0
#           east <- sum(orig.vals[cur.ind + 0:1]) / 2
#         } else if (j == n.cols){
#           # no east face velocity
#           west <-  sum(orig.vals[cur.ind + -1:0]) / 2
#           east <- 0
#         } else {
#           # interior point
#           west <- sum(orig.vals[cur.ind + -1:0]) / 2
#           east <- sum(orig.vals[cur.ind + 0:1]) / 2
#         }
# 
#         # add values to boundary matrix
#         face.vals[cur.ind,] <- c(west, east)
#       }
#     }
#   } else if (direction == "v"){
#     # south to north
#     # name columns appropriately
#     colnames(face.vals) <- c("south", "north")
# 
#     for(i in 1:n.rows){
#       for(j in 1:n.cols){
# 
#         # current indice
#         cur.ind <- (i - 1) * n.cols + j
# 
#         if (i == 1){
#           # no north face velocity
#           north <- 0
#           south <- (orig.vals[cur.ind] + orig.vals[cur.ind + n.cols]) / 2
#         } else if (i == n.rows){
#           # no south face velocity
#           north <-  (orig.vals[cur.ind] + orig.vals[cur.ind - n.cols]) / 2
#           south <- 0
#         } else {
#           # interior point
#           north <-  (orig.vals[cur.ind] + orig.vals[cur.ind - n.cols]) / 2
#           south <- (orig.vals[cur.ind] + orig.vals[cur.ind + n.cols]) / 2
#         }
# 
#         # add values to boundary matrix
#         face.vals[cur.ind,] <- c(south, north)
#       }
#     }
#   }
# 
#   if (with.bound){
#     # return face.vals as is
#     return(face.vals)
# 
#   } else {
#     # remove boundary points from face.vals
#     bound.points <- c(1:n.cols,
#                       (1:(n.rows-2) * n.cols + 1),
#                       (2:(n.rows-1) * n.cols),
#                       ((n.rows-1) * n.cols + 1:n.cols))
# 
#     return(face.vals[-bound.points,])
#   }
# }
# 
# windConvFVM <- function(dims, velocity.faces, inds, bc, bound = FALSE){
#   # Create an advection matrix using the Finite Volume Method (FVM).
#   # Input:
#   #   dims: dimension of the raster
#   #   velocity.faces: interpolated velocities on grid edges (from interpWind)
#   #   inds: matrix indices of raster grid cells (from FVMinds).
#   #   bc: boundary condition
#   #   bound: not used?
#   # Output:
#   #   Sparse FVM matrix representation of advection due to wind.
# 
#   ## rows and column totals
#   n.r <- dims[1]
#   n.c <- dims[2]
# 
#   # total number of cells
#   n.cells <- n.r * n.c
# 
#   if (bc == "insulated"){
# 
#     ### i, j, and x values for interior points
#     n.interior <- nrow(inds$interior)
# 
#     # coefficient calculations
#     F.diff.ew <- velocity.faces[inds$interior[,1], 2] -
#       velocity.faces[inds$interior[,1], 1]
#     F.diff.ns <- velocity.faces[inds$interior[,1], 4] -
#       velocity.faces[inds$interior[,1], 3]
#     F.diff.comb <- F.diff.ew + F.diff.ns
# 
#     a.w <- velocity.faces[inds$interior[,1], 1]; a.w[a.w < 0] <- 0
#     a.e <- -velocity.faces[inds$interior[,1], 2]; a.e[a.e < 0] <- 0
#     a.s <- velocity.faces[inds$interior[,1], 3]; a.s[a.s < 0] <- 0
#     a.n <- -velocity.faces[inds$interior[,1], 4]; a.n[a.n < 0] <- 0
# 
#     sparse.interior <- cbind(rep(inds$interior[,1], 5),
#                              as.vector(inds$interior),
#                              c(F.diff.comb + a.n + a.s + a.e + a.w, -a.n, -a.s, -a.e, -a.w))
# 
# 
#     ### i, j, and x values for west boundary points
#     n.west <- nrow(inds$west)
# 
#     F.diff.ns <- velocity.faces[inds$west[,1], 4] -
#       velocity.faces[inds$west[,1], 3]
#     F.diff.comb <- F.diff.ns + velocity.faces[inds$west[,1], 2]
# 
#     a.e <- -velocity.faces[inds$west[,1], 2]; a.e[a.e < 0] <- 0
#     a.s <- velocity.faces[inds$west[,1], 3]; a.s[a.s < 0] <- 0
#     a.n <- velocity.faces[inds$west[,1], 4]; a.n[a.n > 0] <- 0
# 
#     sparse.west <- cbind(rep(inds$west[,1], 4),
#                          as.vector(inds$west[,-5]),
#                          c(F.diff.comb + a.n + a.s + a.e, -a.n, -a.s, -a.e))
# 
#     ### i, j, and x values for east boundary points
#     n.east <- nrow(inds$east)
# 
#     # coefficient calculations
#     F.diff.ns <- velocity.faces[inds$east[,1], 4] -
#       velocity.faces[inds$east[,1], 3]
#     F.diff.comb <- F.diff.ns - velocity.faces[inds$east[,1], 1]
# 
#     a.w <- velocity.faces[inds$east[,1], 1]; a.w[a.w < 0] <- 0
#     a.s <- velocity.faces[inds$east[,1], 3]; a.s[a.s < 0] <- 0
#     a.n <- -velocity.faces[inds$east[,1], 4]; a.n[a.n < 0] <- 0
# 
#     sparse.east <- cbind(rep(inds$east[,1], 4),
#                          as.vector(inds$east[,-4]),
#                          c(F.diff.comb + a.n + a.s + a.w, -a.n, -a.s, -a.w))
# 
#     ### i, j, and x values for north boundary points
#     n.north <- nrow(inds$north)
# 
#     # coefficient calculations
#     F.diff.ew <- velocity.faces[inds$north[,1], 2] -
#       velocity.faces[inds$north[,1], 1]
#     F.diff.comb <- F.diff.ew - velocity.faces[inds$north[,1], 3]
# 
#     a.w <- velocity.faces[inds$north[,1], 1]; a.w[a.w < 0] <- 0
#     a.e <- -velocity.faces[inds$north[,1], 2]; a.e[a.e < 0] <- 0
#     a.s <- velocity.faces[inds$north[,1], 3]; a.s[a.s < 0] <- 0
# 
#     sparse.north <- cbind(rep(inds$north[,1], 4),
#                           as.vector(inds$north[,-2]),
#                           c(F.diff.comb + a.s + a.e + a.w, -a.s, -a.e, -a.w))
# 
#     ### i, j, and x values for south boundary points
#     n.south <- nrow(inds$south)
# 
#     # coefficient calculations
#     F.diff.ew <- velocity.faces[inds$south[,1], 2] -
#       velocity.faces[inds$south[,1], 1]
#     F.diff.comb <- F.diff.ew + velocity.faces[inds$south[,1], 4]
# 
#     a.w <- velocity.faces[inds$south[,1], 1]; a.w[a.w < 0] <- 0
#     a.e <- -velocity.faces[inds$south[,1], 2]; a.e[a.e < 0] <- 0
#     a.n <- -velocity.faces[inds$south[,1], 4]; a.n[a.n < 0] <- 0
# 
#     sparse.south <- cbind(rep(inds$south[,1], 4),
#                           as.vector(inds$south[,-3]),
#                           c(F.diff.comb + a.n + a.e + a.w, -a.n, -a.e, -a.w))
# 
#     ### i, j, and x values for corner points
# 
#     # northwest values
#     nw.ae <- max(c(0, -velocity.faces[inds$corners[1,1], 2]))
#     nw.as <- max(c(0, velocity.faces[inds$corners[1,1], 3]))
#     nw.ap <- nw.ae + nw.as + velocity.faces[inds$corners[1,1], 2] -
#       velocity.faces[inds$corners[1,1], 3]
# 
#     # northeast values
#     ne.aw <- max(c(0, velocity.faces[inds$corners[2,1], 1]))
#     ne.as <- max(c(0, velocity.faces[inds$corners[2,1], 3]))
#     ne.ap <- ne.aw + ne.as - velocity.faces[inds$corners[2,1], 1] -
#       velocity.faces[inds$corners[2,1], 3]
# 
#     # southwest values
#     sw.ae <- max(c(0, -velocity.faces[inds$corners[3,1], 2]))
#     sw.an <- max(c(0, -velocity.faces[inds$corners[3,1], 4]))
#     sw.ap <- sw.ae + sw.an + velocity.faces[inds$corners[3,1], 2] +
#       velocity.faces[inds$corners[3,1], 4]
# 
#     # southeast values
#     se.aw <- max(c(0, velocity.faces[inds$corners[4,1], 1]))
#     se.an <- max(c(0, -velocity.faces[inds$corners[4,1], 4]))
#     se.ap <- se.aw + se.an - velocity.faces[inds$corners[4,1], 1] +
#       velocity.faces[inds$corners[4,1], 4]
# 
#     sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
#                              inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
#                            c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
#                              inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
#                            c(nw.ap, ne.ap, sw.ap, se.ap, -sw.an, -se.an, -nw.as,
#                              -ne.as, -nw.ae, -sw.ae, -ne.aw, -se.aw))
# 
#     # combine i, j, and x matrices into single matrix
#     sparse.ijx <- rbind(sparse.interior,
#                         sparse.west,
#                         sparse.east,
#                         sparse.north,
#                         sparse.south,
#                         sparse.corner)
# 
#     # quickly fill the sparse matrix
#     A.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
#                           dims = c(n.cells, n.cells))
# 
#   } else if (bc == "endless world"){
# 
#     ### i, j, and x values for interior points
#     n.interior <- nrow(inds$interior)
# 
#     # coefficient calculations
#     F.diff.ew <- velocity.faces[inds$interior[,1], 2] -
#       velocity.faces[inds$interior[,1], 1]
#     F.diff.ns <- velocity.faces[inds$interior[,1], 4] -
#       velocity.faces[inds$interior[,1], 3]
#     F.diff.comb <- F.diff.ew + F.diff.ns
# 
#     a.w <- velocity.faces[inds$interior[,1], 1]; a.w[a.w < 0] <- 0
#     a.e <- -velocity.faces[inds$interior[,1], 2]; a.e[a.e < 0] <- 0
#     a.s <- velocity.faces[inds$interior[,1], 3]; a.s[a.s < 0] <- 0
#     a.n <- -velocity.faces[inds$interior[,1], 4]; a.n[a.n < 0] <- 0
# 
#     sparse.interior <- cbind(rep(inds$interior[,1], 5),
#                              as.vector(inds$interior),
#                              c(F.diff.comb + a.n + a.s + a.e + a.w, -a.n, -a.s, -a.e, -a.w))
# 
# 
#     ### i, j, and x values for west boundary points
#     n.west <- nrow(inds$west)
# 
#     # coefficient values
#     F.diff.ns <- velocity.faces[inds$west[,1], 4] -
#       velocity.faces[inds$west[,1], 3]
#     F.diff.comb <- F.diff.ns + velocity.faces[inds$west[,1], 2]
# 
#     a.e <- -velocity.faces[inds$west[,1], 2]; a.e[a.e < 0] <- 0
#     a.s <- velocity.faces[inds$west[,1], 3]; a.s[a.s < 0] <- 0
#     a.n <- velocity.faces[inds$west[,1], 4]; a.n[a.n > 0] <- 0
# 
#     west.faces <- velocity.faces[inds$west[,1], 1]
#     diag.vals <- F.diff.comb + a.n + a.s + a.e
#     diag.vals[west.faces < 0] <- (diag.vals - west.faces)[west.faces < 0]
# 
#     sparse.west <- cbind(rep(inds$west[,1], 4),
#                          as.vector(inds$west[,-5]),
#                          c(diag.vals, -a.n, -a.s, -a.e))
# 
#     ### i, j, and x values for east boundary points
#     n.east <- nrow(inds$east)
# 
#     # coefficient calculations
#     F.diff.ns <- velocity.faces[inds$east[,1], 4] -
#       velocity.faces[inds$east[,1], 3]
#     F.diff.comb <- F.diff.ns - velocity.faces[inds$east[,1], 1]
# 
#     a.w <- velocity.faces[inds$east[,1], 1]; a.w[a.w < 0] <- 0
#     a.s <- velocity.faces[inds$east[,1], 3]; a.s[a.s < 0] <- 0
#     a.n <- -velocity.faces[inds$east[,1], 4]; a.n[a.n < 0] <- 0
# 
#     # if east face value is positive, change diagonal value to include it leaving system
#     east.faces <- velocity.faces[inds$east[,1],2]
#     diag.vals <- F.diff.comb + a.n + a.s + a.w
#     diag.vals[east.faces > 0] <- (diag.vals + east.faces)[east.faces > 0]
# 
#     sparse.east <- cbind(rep(inds$east[,1], 4),
#                          as.vector(inds$east[,-4]),
#                          c(diag.vals, -a.n, -a.s, -a.w))
# 
#     ### i, j, and x values for north boundary points
#     n.north <- nrow(inds$north)
# 
#     # coefficient calculations
#     F.diff.ew <- velocity.faces[inds$north[,1], 2] -
#       velocity.faces[inds$north[,1], 1]
#     F.diff.comb <- F.diff.ew - velocity.faces[inds$north[,1], 3]
# 
#     a.w <- velocity.faces[inds$north[,1], 1]; a.w[a.w < 0] <- 0
#     a.e <- -velocity.faces[inds$north[,1], 2]; a.e[a.e < 0] <- 0
#     a.s <- velocity.faces[inds$north[,1], 3]; a.s[a.s < 0] <- 0
# 
#     north.faces <- velocity.faces[inds$north[,1], 4]
#     diag.vals <- F.diff.comb + a.s + a.e + a.w
#     diag.vals[north.faces > 0] <- (diag.vals + north.faces)[north.faces > 0]
# 
#     sparse.north <- cbind(rep(inds$north[,1], 4),
#                           as.vector(inds$north[,-2]),
#                           c(diag.vals, -a.s, -a.e, -a.w))
# 
#     ### i, j, and x values for south boundary points
#     n.south <- nrow(inds$south)
# 
#     # coefficient calculations
#     F.diff.ew <- velocity.faces[inds$south[,1], 2] -
#       velocity.faces[inds$south[,1], 1]
#     F.diff.comb <- F.diff.ew + velocity.faces[inds$south[,1], 4]
# 
#     a.w <- velocity.faces[inds$south[,1], 1]; a.w[a.w < 0] <- 0
#     a.e <- -velocity.faces[inds$south[,1], 2]; a.e[a.e < 0] <- 0
#     a.n <- -velocity.faces[inds$south[,1], 4]; a.n[a.n < 0] <- 0
# 
#     south.faces <- velocity.faces[inds$south[,1], 3]
#     diag.vals <- F.diff.comb + a.n + a.e + a.w
#     diag.vals[south.faces < 0] <- (diag.vals - south.faces)[south.faces < 0]
# 
#     sparse.south <- cbind(rep(inds$south[,1], 4),
#                           as.vector(inds$south[,-3]),
#                           c(F.diff.comb + a.n + a.e + a.w, -a.n, -a.e, -a.w))
# 
#     ### i, j, and x values for corner points
# 
#     ## northwest values
#     nw.ae <- max(c(0, -velocity.faces[inds$corners[1,1], 2]))
#     nw.as <- max(c(0, velocity.faces[inds$corners[1,1], 3]))
#     nw.ap <- nw.ae + nw.as + velocity.faces[inds$corners[1,1], 2] -
#       velocity.faces[inds$corners[1,1], 3]
# 
#     # north face
#     nw.nf <- velocity.faces[inds$corners[1,1], 4]
#     # west face
#     nw.wf <- velocity.faces[inds$corners[1,1], 1]
# 
#     if (nw.nf > 0){
#       # add to diagonal
#       nw.ap <- nw.ap + nw.nf
#     }
# 
#     if (nw.wf < 0){
#       # subtract from diagonal
#       nw.ap <- nw.ap - nw.wf
#     }
# 
#     ## northeast values
#     ne.aw <- max(c(0, velocity.faces[inds$corners[2,1], 1]))
#     ne.as <- max(c(0, velocity.faces[inds$corners[2,1], 3]))
#     ne.ap <- ne.aw + ne.as - velocity.faces[inds$corners[2,1], 1] -
#       velocity.faces[inds$corners[2,1], 3]
# 
#     # north face
#     ne.nf <- velocity.faces[inds$corners[2,1], 4]
#     # east face
#     ne.ef <- velocity.faces[inds$corners[2,1], 2]
# 
#     if (ne.nf > 0){
#       ne.ap <- ne.ap + ne.nf
#     }
# 
#     if (ne.ef > 0){
#       ne.ap <- ne.ap + ne.ef
#     }
# 
# 
#     ## southwest values
#     sw.ae <- max(c(0, -velocity.faces[inds$corners[3,1], 2]))
#     sw.an <- max(c(0, -velocity.faces[inds$corners[3,1], 4]))
#     sw.ap <- sw.ae + sw.an + velocity.faces[inds$corners[3,1], 2] +
#       velocity.faces[inds$corners[3,1], 4]
# 
#     # south face
#     sw.sf <- velocity.faces[inds$corners[3,1],3]
#     # west face
#     sw.wf <- velocity.faces[inds$corners[3,1],1]
# 
#     if (sw.sf < 0){
#       sw.ap <- sw.ap - sw.sf
#     }
# 
#     if (sw.wf < 0){
#       sw.ap <- sw.ap - sw.wf
#     }
# 
#     ## southeast values
#     se.aw <- max(c(0, velocity.faces[inds$corners[4,1], 1]))
#     se.an <- max(c(0, -velocity.faces[inds$corners[4,1], 4]))
#     se.ap <- se.aw + se.an - velocity.faces[inds$corners[4,1], 1] +
#       velocity.faces[inds$corners[4,1], 4]
# 
#     # south face
#     se.sf <- velocity.faces[inds$corners[4,1], 3]
#     # east face
#     se.ef <- velocity.faces[inds$corners[4,1], 2]
# 
#     if (se.sf < 0){
#       se.ap <- se.ap - se.sf
#     }
# 
#     if (se.ef > 0){
#       se.ap <- se.ap + se.ef
#     }
# 
#     # combine all values together
#     sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
#                              inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
#                            c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
#                              inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
#                            c(nw.ap, ne.ap, sw.ap, se.ap, -sw.an, -se.an, -nw.as,
#                              -ne.as, -nw.ae, -sw.ae, -ne.aw, -se.aw))
# 
#     # combine i, j, and x matrices into single matrix
#     sparse.ijx <- rbind(sparse.interior,
#                         sparse.west,
#                         sparse.east,
#                         sparse.north,
#                         sparse.south,
#                         sparse.corner)
# 
#     # quickly fill the sparse matrix
#     A.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
#                           dims = c(n.cells, n.cells))
#   }
# 
#   return(A.mat)
# }
# 
# diffFVM <- function(dims, inds, bound.cond){
#   # Create finite volume matrix approximation to homogeneous diffusion.
#   # Input:
#   #   dims: dimension of raster
#   #   inds: indices of interior/boundary points, created with FVMInds
#   #   bound.cond: boundary conditions (insulated, periodic, Neumann)
#   # Output:
#   #   Sparse matrix used for discrete approximation of homogenous diffusion.
#   
#   # rows and column totals
#   n.r <- dims[1]
#   n.c <- dims[2]
#   
#   # total number of cells
#   n.cells <- n.r * n.c
#   
#   if (bound.cond == "insulated"){
#     
#     # i, j, and x values for interior points
#     n.interior <- nrow(inds$interior)
#     sparse.interior <- cbind(rep(inds$interior[,1], 5), 
#                              as.vector(inds$interior),
#                              c(rep(4, n.interior), rep(-1, n.interior * 4)))
#     
#     # i, j, and x values for west boundary points
#     n.west <- nrow(inds$west)
#     sparse.west <- cbind(rep(inds$west[,1], 4),
#                          as.vector(inds$west),
#                          c(rep(3, n.west), rep(-1, n.west * 3)))
#     
#     # i, j, and x values for east boundary points
#     n.east <- nrow(inds$east)
#     sparse.east <- cbind(rep(inds$east[,1], 4),
#                          as.vector(inds$east),
#                          c(rep(3, n.east), rep(-1, n.east * 3)))
#     
#     # i, j, and x values for north boundary points
#     n.north <- nrow(inds$north)
#     sparse.north <- cbind(rep(inds$north[,1], 4),
#                           as.vector(inds$north),
#                           c(rep(3, n.north), rep(-1, n.north * 3)))
#     
#     # i, j, and x values for south boundary points
#     n.south <- nrow(inds$south)
#     sparse.south <- cbind(rep(inds$south[,1], 4),
#                           as.vector(inds$south),
#                           c(rep(3, n.south), rep(-1, n.south * 3)))
#     
#     # i, j, and x values for corner points
#     
#     sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
#                              inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
#                            c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
#                              inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
#                            c(rep(2, 4), rep(-1, 8)))
#     
#     # combine i, j, and x matrices into single matrix
#     sparse.ijx <- rbind(sparse.interior,
#                         sparse.west,
#                         sparse.east,
#                         sparse.north,
#                         sparse.south,
#                         sparse.corner)
#     
#   } else if (bound.cond == "periodic"){
#     
#     # i, j, and x values for interior points
#     n.interior <- nrow(inds$interior)
#     sparse.interior <- cbind(rep(inds$interior[,1], 5), 
#                              as.vector(inds$interior),
#                              c(rep(4, n.interior), rep(-1, n.interior * 4)))
#     
#     
#     # i, j, and x values for west boundary points
#     n.west <- nrow(inds$west)
#     sparse.west <- cbind(rep(inds$west[,1], 5),
#                          as.vector(inds$west),
#                          c(rep(4, n.west), rep(-1, n.west * 4)))
#     
#     # i, j, and x values for east boundary points
#     n.east <- nrow(inds$east)
#     sparse.east <- cbind(rep(inds$east[,1], 5),
#                          as.vector(inds$east),
#                          c(rep(4, n.east), rep(-1, n.east * 4)))
#     
#     # i, j, and x values for north boundary points
#     n.north <- nrow(inds$north)
#     sparse.north <- cbind(rep(inds$north[,1], 5),
#                           as.vector(inds$north),
#                           c(rep(4, n.north), rep(-1, n.north * 4)))
#     
#     # i, j, and x values for south boundary points
#     n.south <- nrow(inds$south)
#     sparse.south <- cbind(rep(inds$south[,1], 5),
#                           as.vector(inds$south),
#                           c(rep(4, n.south), rep(-1, n.south * 4)))
#     
#     # i, j, and x values for corner points
#     sparse.corner <- cbind(rep(inds$corners[,1], 5),
#                            as.vector(inds$corners),
#                            c(rep(4, 4), rep(-1, 4 * 4)))
#     
#     # combine i, j, and x matrices into single matrix
#     sparse.ijx <- rbind(sparse.interior,
#                         sparse.west,
#                         sparse.east,
#                         sparse.north,
#                         sparse.south,
#                         sparse.corner)
#     
#   } else if (bound.cond == "Neumann"){
# 
#     # i, j, and x values for interior points
#     n.interior <- nrow(inds$interior)
#     sparse.interior <- cbind(rep(inds$interior[,1], 5), 
#                              as.vector(inds$interior),
#                              c(rep(4, n.interior), rep(-1, n.interior * 4)))
#     
#     # i, j, and x values for west boundary points
#     n.west <- nrow(inds$west)
#     sparse.west <- cbind(rep(inds$west[,1], 4),
#                          as.vector(inds$west),
#                          c(rep(4, n.west), rep(-1, n.west * 2), rep(-2, n.west)))
#     
#     # i, j, and x values for east boundary points
#     n.east <- nrow(inds$east)
#     sparse.east <- cbind(rep(inds$east[,1], 4),
#                          as.vector(inds$east),
#                          c(rep(4, n.east), rep(-1, n.east * 2), rep(-2, n.east)))
#     
#     # i, j, and x values for north boundary points
#     n.north <- nrow(inds$north)
#     sparse.north <- cbind(rep(inds$north[,1], 4),
#                           as.vector(inds$north),
#                           c(rep(4, n.north), rep(-2, n.north), rep(-1, n.north * 2)))
#     
#     # i, j, and x values for south boundary points
#     n.south <- nrow(inds$south)
#     sparse.south <- cbind(rep(inds$south[,1], 4),
#                           as.vector(inds$south),
#                           c(rep(4, n.south), rep(-2, n.south), rep(-1, n.south * 2)))
#     
#     # i, j, and x values for corner points
#     sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
#                              inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
#                            c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
#                              inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
#                            c(rep(4, 4), rep(-2, 8)))
#     
#     # combine i, j, and x matrices into single matrix
#     sparse.ijx <- rbind(sparse.interior,
#                         sparse.west,
#                         sparse.east,
#                         sparse.north,
#                         sparse.south,
#                         sparse.corner)
#   }
#   
#   # quickly fill the sparse matrix
#   D.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
#                         dims = c(n.cells, n.cells))
#   
#   # return the sparse diffusion matrix
#   return(D.mat)
# }
# 
# FVMInds <- function(dims, bound.cond){
#   # Create list denoting index for interior and boundary of raster with given
#   #   spatial extent. This is used to create the FVM approximation of advection
#   #   due to wind.
#   # Input:
#   #   dims: dimension of raster
#   #   bound.cond: type of boundary condition (insulated or periodic).
#   # Output:
#   #   List with index of interior and boundary points.
#   
#   # rows and column totals
#   n.r <- dims[1]
#   n.c <- dims[2]
#   
#   #######################
#   ### interior points ###
#   #######################
#   interior <- matrix(NA_real_, nrow = (n.c - 2) * (n.r - 2), ncol = 5)
#   colnames(interior) <- c("P", "N", "S", "E", "W")
#   
#   if (bound.cond == "insulated"){
#     # create a diffusion matrix with insulated boundaries
#     
#     cur.ind <- 1
#     for (i in 2:(n.r-1)){
#       for (j in 2:(n.c - 1)){
#         
#         # index in original scheme
#         orig.ind <- (i - 1) * n.c + j
#         
#         # indices (with neighbors) for interior
#         interior[cur.ind, 1] <- orig.ind
#         interior[cur.ind, 2] <- orig.ind - n.c
#         interior[cur.ind, 3] <- orig.ind + n.c
#         interior[cur.ind, 4] <- orig.ind + 1
#         interior[cur.ind, 5] <- orig.ind - 1
#         
#         cur.ind <- cur.ind + 1
#       }
#     }
#     
#     ### west boundary points
#     west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
#     colnames(west.bound) <- c("P", "N", "S", "E")
#     
#     cur.ind <- 1
#     for (i in 2:(n.r - 1)){
#       # index in original scheme
#       orig.ind <- (i - 1) * n.c + 1
#       # indices (with neighbors) for west boundary
#       west.bound[cur.ind, 1] <- orig.ind
#       west.bound[cur.ind, 2] <- orig.ind - n.c
#       west.bound[cur.ind, 3] <- orig.ind + n.c
#       west.bound[cur.ind, 4] <- orig.ind + 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ### east boundary
#     east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
#     colnames(east.bound) <- c("P", "N", "S", "W")
#     
#     cur.ind <- 1
#     for (i in 2:(n.r-1)){
#       
#       # index in original scheme
#       orig.ind <- i * n.c
#       
#       # indices (with neighbors) for interior
#       east.bound[cur.ind, 1] <- orig.ind
#       east.bound[cur.ind, 2] <- orig.ind - n.c
#       east.bound[cur.ind, 3] <- orig.ind + n.c
#       east.bound[cur.ind, 4] <- orig.ind - 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ### north boundary
#     north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
#     colnames(north.bound) <- c("P", "S", "E", "W")
#     
#     cur.ind <- 1
#     for (j in 2:(n.c-1)){
#       
#       # index in original scheme
#       orig.ind <- j
#       
#       # indices (with neighbors) for interior
#       north.bound[cur.ind, 1] <- orig.ind
#       north.bound[cur.ind, 2] <- orig.ind + n.c
#       north.bound[cur.ind, 3] <- orig.ind + 1
#       north.bound[cur.ind, 4] <- orig.ind - 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ### south boundary
#     south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
#     colnames(south.bound) <- c("P", "N", "E", "W")
#     
#     cur.ind <- 1
#     for (j in 2:(n.c-1)){
#       
#       # index in original scheme
#       orig.ind <- j + n.c * (n.r - 1)
#       
#       # indices (with neighbors) for interior
#       south.bound[cur.ind, 1] <- orig.ind
#       south.bound[cur.ind, 2] <- orig.ind - n.c
#       south.bound[cur.ind, 3] <- orig.ind + 1
#       south.bound[cur.ind, 4] <- orig.ind - 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ### corners
#     nw <- c(1, NA_real_, 1 + n.c, 2, NA_real_)
#     ne <- c(n.c, NA_real_, 2 * n.c, NA_real_, n.c - 1)
#     sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, NA_real_, (n.r - 1) * n.c + 2, NA_real_)
#     se <- c(n.r * n.c, (n.r - 1) * n.c, NA_real_, NA_real_, n.r * n.c - 1)
#     
#     corners <- rbind(nw, ne, sw, se)
#     colnames(corners) <- c("P", "N", "S", "E", "W")
#     
#   } else if (bound.cond == "periodic"){
#     
#     # create a diffusion matrix with periodic boundaries
#     
#     cur.ind <- 1
#     for (i in 2:(n.r-1)){
#       for (j in 2:(n.c - 1)){
#         
#         # index in original scheme
#         orig.ind <- (i - 1) * n.c + j
#         
#         # indices (with neighbors) for interior
#         interior[cur.ind, 1] <- orig.ind
#         interior[cur.ind, 2] <- orig.ind - n.c
#         interior[cur.ind, 3] <- orig.ind + n.c
#         interior[cur.ind, 4] <- orig.ind + 1
#         interior[cur.ind, 5] <- orig.ind - 1
#         
#         cur.ind <- cur.ind + 1
#       }
#     }
#     
#     ############################
#     ### west boundary points ###
#     ############################
#     west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
#     colnames(west.bound) <- c("P", "N", "S", "E", "W")
#     
#     cur.ind <- 1
#     for (i in 2:(n.r - 1)){
#       # index in original scheme
#       orig.ind <- (i - 1) * n.c + 1
#       # indices (with neighbors) for west boundary
#       west.bound[cur.ind, 1] <- orig.ind
#       west.bound[cur.ind, 2] <- orig.ind - n.c
#       west.bound[cur.ind, 3] <- orig.ind + n.c
#       west.bound[cur.ind, 4] <- orig.ind + 1
#       west.bound[cur.ind, 5] <- orig.ind + (n.c - 1)
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     #####################
#     ### east boundary ###
#     #####################
#     east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
#     colnames(east.bound) <- c("P", "N", "S", "E", "W")
#     
#     cur.ind <- 1
#     for (i in 2:(n.r-1)){
#       
#       # index in original scheme
#       orig.ind <- i * n.c
#       
#       # indices (with neighbors) for interior
#       east.bound[cur.ind, 1] <- orig.ind
#       east.bound[cur.ind, 2] <- orig.ind - n.c
#       east.bound[cur.ind, 3] <- orig.ind + n.c
#       east.bound[cur.ind, 4] <- orig.ind - (n.c - 1)
#       east.bound[cur.ind, 5] <- orig.ind - 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ######################
#     ### north boundary ###
#     ######################
#     north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
#     colnames(north.bound) <- c("P", "N", "S", "E", "W")
#     
#     cur.ind <- 1
#     for (j in 2:(n.c-1)){
#       
#       # index in original scheme
#       orig.ind <- j
#       
#       # indices (with neighbors) for interior
#       north.bound[cur.ind, 1] <- orig.ind
#       north.bound[cur.ind, 2] <- orig.ind + (n.r - 1) * n.c
#       north.bound[cur.ind, 3] <- orig.ind + n.c
#       north.bound[cur.ind, 4] <- orig.ind + 1
#       north.bound[cur.ind, 5] <- orig.ind - 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ######################
#     ### south boundary ###
#     ######################
#     south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
#     colnames(south.bound) <- c("P", "N", "S", "E", "W")
#     
#     cur.ind <- 1
#     for (j in 2:(n.c-1)){
#       
#       # index in original scheme
#       orig.ind <- j + n.c * (n.r - 1)
#       
#       # indices (with neighbors) for interior
#       south.bound[cur.ind, 1] <- orig.ind
#       south.bound[cur.ind, 2] <- orig.ind - n.c
#       south.bound[cur.ind, 3] <- orig.ind - (n.r - 1) * n.c
#       south.bound[cur.ind, 4] <- orig.ind + 1
#       south.bound[cur.ind, 5] <- orig.ind - 1
#       
#       cur.ind <- cur.ind + 1
#     }
#     
#     ###############
#     ### corners ###
#     ###############
#     nw <- c(1, 1 + (n.r - 1) * n.c, 1 + n.c, 2, n.c)
#     ne <- c(n.c, n.c * n.r, 2 * n.c, 1, n.c - 1)
#     sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, 1, (n.r - 1) * n.c + 2, n.r * n.c)
#     se <- c(n.r * n.c, (n.r - 1) * n.c, n.c, 1 + (n.r - 1) * n.c, n.r * n.c - 1)
#     
#     corners <- rbind(nw, ne, sw, se)
#     colnames(corners) <- c("P", "N", "S", "E", "W")
#   }
#   
#   # combine all indice matrices into a single list
#   indices <- list(interior = interior,
#                   north = north.bound,
#                   south = south.bound,
#                   east = east.bound,
#                   west = west.bound,
#                   corners = corners)
#   
#   # return indices
#   return(indices)
# }
# 
# createC <- function(big, small, dims, inds, bound.cond = "endless world"){
#   # Create matrix fo discretized advection process, due to wind.
#   # Input:
#   #   big: raster with wind data, large spatial extent
#   #   small: raster with wind data, small spatial extent
#   #   dims: dimension of the raster
#   #   inds: list of indices, created with FVMInds
#   #   bound.cond: boundary condition (periodic or endless world)
#   # Output:
#   #   Sparse matrix representing advection due to wind, as specified by 
#   #     wind velocities found in 'small'.
#   
#   big.uwind <- big$uwind
#   big.vwind <- big$vwind
#   
#   small.uwind <- small$uwind
#   small.vwind <- small$vwind
#   
#   u.faces <- interpWind(big.uwind, small.uwind, "u", with.bound = FALSE)
#   v.faces <- interpWind(big.vwind, small.vwind, "v", with.bound = FALSE)
#   velocity <- cbind(u.faces, v.faces)
#   
#   # wind vmatrix
#   C <- windConvFVM(dims, velocity, inds, bc = bound.cond, bound = FALSE)
#   return(C)
# }
# 
# ################################################################################
# ### 3. Functions used to compile results.
# ################################################################################
# 
# dic.calc <- function(y, X, burnin, samples, dic.vec, mats, model, intercept = FALSE){
#   # Calculate DIC using posterior output for a given type of model.
#   # Input:
#   #   y: response
#   #   X: covariates
#   #   burnin: burnin level
#   #   samples: samples (from mcmc output)
#   #   div.vec: vector of likelihood values (from mcmc output)
#   #   mats: list of FVM matrices
#   #   model: which likelihood to use ("avg2", "avg1", "snap2", "snap1")
#   # Output:
#   #   DIC for given posterior samples.
#   
#   D.bar <- mean(dic.vec[-burnin])
#   theta.hat <- colMeans(samples[-burnin,])
#   
#   if (model == "avg2"){
#     D.hat <- -2 * ll.data.Wind(theta.hat, y, X, mats, int = intercept)
#   } else if (model == "avg1"){
#     D.hat <- -2 * ll.data.uncoupled.W(theta.hat, y, X, mats)
#   } else if (model == "snap2"){
#     D.hat <- -2 * ll.data.snap(theta.hat, y, X, mats)
#   } else {
#     D.hat <- -2 * ll.data.uncoupled.snap(theta.hat, y, X, mats)
#   }
#   
#   DIC <- (2 * D.bar) - D.hat
#   return(DIC)
# }
# 
# fittedMeanSO4Surface <- function(theta, mats, X){
#   # Estimate the expected annual SO4 concentration, using parameter estimates
#   #   fit to the coupled, time-averaged model.
#   # Input:
#   #   theta: parameter values
#   #   mats: list of FVM matrices
#   #   X: vector of SO2 emissions
#   # Output:
#   #   Vector of estimated SO4 concentrations.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
#   
#   D <- mats$D
#   C <- mats$C
#   
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#   
#   
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
#   return(as.vector(mu))
# }
# 
# simCoupled.Wind <- function(theta, mats, X){
#   # Simulate from the time-averaged, coupled model, using given parameter and 
#   #   X values.
#   # Input:
#   #   theta: model parameters
#   #   mats: list of FVM matrices
#   #   X: vector of SO2 emissions
#   # Output:
#   #   A single realization of annual SO4 concentrations simulated from the 
#   #     time-averaged, coupled model.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
#   
#   D <- mats$D
#   C <- mats$C
#   
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#   
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
#   
#   Q <- t(A.1) %*% A.1 / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# simUncoupled.W <- function(theta, mats, X){
#   # Simulate from the time-averaged, uncoupled model, using given parameter and 
#   #   X values.
#   # Input:
#   #   theta: model parameters
#   #   mats: list of FVM matrices
#   #   X: vector of SO2 emissions
#   # Output:
#   #   A single realization of annual SO4 concentrations simulated from the 
#   #     time-averaged, uncoupled model.
#   
#   # grab parameters
#   gamma <- theta[1]
#   beta <- theta[2]
#   s2 <- theta[3]
#   delta <- theta[4]
#   alpha <- theta[5]
#   
#   D <- mats$D
#   C <- mats$C
#   
#   n.cells <- nrow(D)
#   A <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   
#   mu <- solve(A, X * beta)
#   
#   Q <- (t(A) %*% A) / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# simCoupled.snap <- function(theta, mats, X){
#   # Simulate from the snapshot, coupled model, using given parameter and 
#   #   X values.
#   # Input:
#   #   theta: model parameters
#   #   mats: list of FVM matrices
#   #   X: vector of SO2 emissions
#   # Output:
#   #   A single realization of annual SO4 concentrations simulated from the 
#   #     snapshot, coupled model.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   
#   D <- mats$D
#   
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + Diagonal(n.cells, 1)
#   
#   
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
#   
#   Q <- 2 * A.1 / s2
#   
#   sample <- rnorm.Q(Q, mu, zero.constraint = FALSE, canon = FALSE)
#   
#   return(as.vector(sample))
# }
# 
# 
# expCalc <- function(N, X.orig, X.new, samples, burnin = 0, mats.exp, pop.r){
#   # Calculate population-normalized exposure to SO4 pollution using draws from
#   #   posterior likelihood.
#   # Input:
#   #   N: number of posterior draws/simulations
#   #   X.orig: vector of original SO2 emissions levels
#   #   X.new: vector of hypothetical SO2 emissions after scrubber
#   #   samples: list of posterior samples, from mcmc output
#   #   burnin: size of burnin
#   #   mats.exp: FVM matrices used in discretization
#   #   pop.r: raster of population density
#   # Output:
#   #   Estimated population-normalized exposure to SO4 pollution.
#   
#   if (burnin > 0){
#     samps <- samples[-c(1:burnin),]
#   } else {
#     samps <- samples
#   }
# 
#   # sample indices
#   s.inds <- sample.int(nrow(samps), size = N, replace = FALSE)
#   
#   log.pop.vec <- log(values(pop.r))
#   norm.pop <- exp(log.pop.vec - logSumExp(log.pop.vec))
#   
#   exposure <- matrix(NA_real_, nrow = N, ncol = ncol(X.new) + 1)
#   
#   for (k in 1:N){
#     
#     # sample theta
#     theta.k <- samples[s.inds[k],]
#     
#     # exposure in original system
#     orig.sim <- simCoupled.Wind(theta = theta.k, mats = mats.exp, X = X.orig)
#     exposure[,1] <- sum(orig.sim * norm.pop)
#     
#     for (j in 1:ncol(X.new)){
#       new.sim <- simCoupled.Wind(theta = theta.k, mats = mats.exp, X = X.new[,j])
#       exposure[k,1 + j] <- sum(new.sim * norm.pop)
#     }
#   }
#   
#   return(exposure)
# }
# 
# breakpoint.creation <- function(ras.obj, n, min.val = NA, max.val = NA){
#     # This function creates a similar color scale for plotting across rasters.
#     # Input:
#     #   ras.obj: a raster object to be plotted
#     #   n: number of color breaks
#     #   min.val: minimum value of color scale (default = NA)
#     #   max.val: maximum value of color scale (default = NA)
# 
#     if (is.na(min.val)){
#         min.val <- min(values(ras.obj))
#     }
# 
#     if (is.na(max.val)){
#         max.val <- max(values(ras.obj))
#     }
# 
#     seq(from = min.val, to = max.val, length.out = n + 1)
# }
# 
# # creates a list of colors from #F2F2F2FF to #F1DEDCFF
# colfunc1 <- colorRampPalette(c("#F2F2F2FF", "#F1DEDCFF"))
# 
# # creates a list of blue colors, from #F2F2F2FF to #08306B
# colfunc2 <- function(n){
#     new.blues <- c("#F2F2F2FF", "#F7FBFF", "#DEEBF7", "#C6DBEF", 
#         "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
# 
#     colorRampPalette(new.blues)(n)
# }
# 
# indep.mean <- function(samples, X, mats, burnin = 0){
#   # Creates surface plots from posterior samples, assuming independent
#   #   space-time white noise for the driving process of the SDE
#   
#   # determine posterior means:
#   if (burnin > 0){
#     theta.hat <- colMeans(samples[-c(1:burnin),])
#   } else {
#     theta.hat <- colMeans(samples)
#   }
# 
#   # grab parameters
#   gamma <- theta.hat[1]
#   eta <- theta.hat[2]
#   beta <- theta.hat[3]
#   s2 <- theta.hat[4]
#   delta <- theta.hat[5]
#   alpha <- theta.hat[6]
#   
#   D <- mats$D
#   C <- mats$C
#   
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / eta) * D + (alpha / eta) * C + Diagonal(n.cells, 1)
#   
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
# 
#   return(mu)
# }
# 
# 
# ################################################################################
# ### 4. Functions for fitting phenomenological model.
# ################################################################################
# 
# phenom.mean <- function(theta, D, x, u, v, vrsn = 1){
#   # Mean function for our phenomenological model, which is used for comparison.
#   #   The surface is composed of Gaussian functions, centered and scaled by
#   #   power plant locations and output. 
#   # Input:
#   #   theta: vector of parameters
#   #   D: distance matrix
#   #   x: vector of power plant outputs mapped to each raster element
#   #   u: vector of u-wind components
#   #   v: vector of v-wind components
#   #   vrsn: type of phenom. mean structure 
#   #     (1 = constant variance, 2 = variance scaled by power plant output)
#   # Output:
#   #   Mean surface of phenomenological model; vector mapped to raster.
# 
#   # parameters
#   alpha <- theta[1]
#   s2.m <- theta[2]  
#   beta.0 <- theta[3]
#   beta.1 <- theta[4]
#   beta.2 <- theta[5]
# 
#   # location of power plants
#   locs <- which(x > 0)
#   
#   if (vrsn == 1){
#     # Gaussian surface, with constant variance for each Gaussian function
#     f <- colSums(alpha * sqrt(x[locs] / 1000) * exp(-D[locs,]^2 / 2 / s2.m ))
#   } else {
#     # Gaussian surface, with variance scaled by power plant output
#     f <- colSums(alpha * sqrt(x[locs] / 1000) * exp(-D[locs,]^2 / 2 / s2.m / x[locs]))
#   }
# 
#   # mean surface
#   mu <- rep(beta.0, length(f)) + f + beta.1 * v + beta.2 * u
# 
#   return(mu)
# }
# 
# phenom.prec <- function(theta, n, W.p = NA, M.inv = NA, type = "iid"){
#   # Create the precision matrix of the phenomenological model.
#   # Input:
#   #   theta: vector of parameters
#   #   n: size of n x n precision matrix
#   #   W.p: weight matrix
#   #   M.inv: adjacency matrix (for CAR and SAR construction)
#   #   type: type of covariance (options: "iid", "CAR", "SAR")
#   # Output:
#   #   Precision matrix as defined by "type".
# 
#   if (type == "iid"){
#     # create IID precision matrix
#     s2 <- theta[1]
#     Omega <- Diagonal(n, 1 / s2)
#   } else if (type == "CAR"){
#     # create CAR precision matrix
#     s2 <- theta[1]
#     rho <- theta[2]
# 
#     I.mat <- Diagonal(n)
#     Omega <- M.inv %*% (I.mat - rho * W.p) / s2
# 
#   } else if (type == "SAR"){
#     # create SAR precision matrix
#     s2 <- theta[1]
#     rho <- theta[2]
# 
#     I.mat <- Diagonal(n)
#     Omega <- (I.mat - rho * t(W.p)) %*%  (I.mat - rho * W.p) / s2
#   }
# 
#   return(Omega)
# }
# 
# inBounds <- function(theta, type){
#   # Determine if covariance parameters have proper support.
#   # Input:
#   #   theta: vector of parameters
#   #   type: type of covariance (options = "iid", "CAR", or "SAR")
#   # Output:
#   #   Boolean specifying if proposed parameters have proper support.
# 
#   # parameters
#   alpha <- theta[1]
#   s2.m <- theta[2]  
#   beta.0 <- theta[3]
#   beta.1 <- theta[4]
#   beta.2 <- theta[5]
# 
#   # type of covariance
#   if (type == "iid"){
#     s2.cov <- theta[6]
#     rho <- 0.5
#   } else {
#     s2.cov <- theta[6]
#     rho <- theta[7]
#   }
# 
#   # check if parameters have proper support
#   if (alpha < 0){
#     return(FALSE)
#   } else if (s2.m < 0){
#     return(FALSE)
#   } else if (s2.cov < 0){
#     return(FALSE)
#   } else if (rho <= 0 | rho >= 1){
#     return(FALSE)
#   } else {
#     return(TRUE)
#   }
# }
# 
# phenomLikelihood <- function(theta, D, y, x, u, v, W.plus = NA, M.inv = NA, type = "iid", vrsn = 1){
#   # Calculate the likelihood of the phenomenological model.
#   # Input: 
#   #   theta: vector of parameters
#   #   D: distance matrix
#   #   y: response vector
#   #   x: vector of power plant outputs mapped to each raster element
#   #   u: vector of u-wind components
#   #   v: vector of v-wind components
#   #   W.plus: weights matrix (covar)
#   #   M.inv: adjacency matrix (covar)
#   #   type: type of covariance ("iid", "CAR", "SAR")
#   #   vrsn: type of phenom. mean structure 
#   #     (1 = constant variance, 2 = variance scaled by power plant output)
#   # Output: 
#   #   Log-likelihood of phenomenological model (vrsn = 1 or 2).
#   
#   # number of responses
#   n.y <- length(y)
#   
#   # split parameters
#   theta.mu <- theta[1:5]
#   theta.cov <- theta[-c(1:5)]
# 
#   # mean
#   mu <- phenom.mean(theta.mu, D = D, x = x, u = u, v = v, vrsn = vrsn)  
# 
#   # precision matrix
#   Omega <- phenom.prec(theta.cov, n = n.y, W.p = W.plus, M.inv = M.inv, type = type)
# 
#   # log-likelihood:
#   if(inBounds(theta, type)){
#     resid <- as.vector(y - mu)
#     QF <- t(resid) %*% (Omega %*% resid)
#     ll <- as.numeric(-(n.y * log(2 * pi) / 2) + (as.numeric(determinant(Omega)$modulus) / 2) - (QF / 2))
#   } else {
#     ll <- -Inf
#   }
#   return(ll)
# }
# 
# phenomFit <- function(usa, theta.0, cov.type = "iid", mean.vrsn = 1, ...) {
#   # Maximum likelihood optimization of the phenomenological model.
#   # Input:
#   #   usa: data raster
#   #   theta.0: vector with initial parameter values
#   #   cov.type: type of covariance
#   #   mean.vrsn:  type of phenom. mean structure
#   #     (1 = constant variance, 2 = variance scaled by power plant output)
#   # Output:
#   #   MLE parameter estimates.
# 
#   ## set-up
# 
#   # X and y values
#   y.2011 <- values(usa$so4$so4.2011)
#   X.2011 <- usa$X$X.2011
# 
#   # wind
#   u <- values(usa$wind.small$wind.2011$uwind)
#   v <- values(usa$wind.small$wind.2011$vwind)
# 
#   # distance
#   r.coords <- coordinates(usa$so4$so4.2011)
#   D.mat <- distm(r.coords, fun = distHaversine) / 1000
# 
#   if (cov.type == "iid") {
#     # assume independent variance
#     W.plus = NA
#     M.inv = NA
#   } else {
#     # weights matrix for autoregressive covariances
#     adj <- raster::adjacent(usa$so4$so4.2011, cells = 1:length(y.2011), pairs = TRUE, directions = 4)
#     W <- Matrix::sparseMatrix(i = adj[, 1], j = adj[, 2], x = 1)
#     M.inv <- Matrix::Diagonal(n = nrow(W), x = rowSums(W))
#     W.plus <- solve(M.inv, W)
#   }
# 
#   # call to optim
#   mle <- optim(par = theta.0, function(x) {
#     -phenomLikelihood(x, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, mean.vrsn)
#   }, ...)
# 
#   # result
#   return(mle)
# }
# 
# phenom.dic <- function(usa, mcmc.res, burnin = 0, cov.type, mean.vrsn = 2){
#   # Calculate the DIC using MCMC samples from the phenomenological model.
#   # Input: 
#   #   usa: data raster
#   #   mcmc.res: list with MCMC samples
#   #   burnin: size of burning
#   #   cov.type: type of covariance 
#   #   mean.vrsn: type of phenom. mean structure 
#   #     (1 = constant variance, 2 = variance scaled by power plant output)
#   # Output:
#   #   DIC for the phenomenological model fit. 
# 
#   ## set-up
# 
#   # X and y values
#   y.2011 <- values(usa$so4$so4.2011)
#   X.2011 <- usa$X$X.2011
# 
#   # wind
#   u <- values(usa$wind.small$wind.2011$uwind)
#   v <- values(usa$wind.small$wind.2011$vwind)
# 
#   # distance
#   r.coords <- coordinates(usa$so4$so4.2011)
#   D.mat <- distm(r.coords, fun = distHaversine) / 1000  
# 
#   if (cov.type == "iid"){
#     # assume independent variance
#     W.plus = NA; M.inv = NA
#   } else {
#     # weights matrix for autoregressive covariances
#     adj <- raster::adjacent(usa$so4$so4.2011, cells = 1:length(y.2011), pairs = TRUE, directions = 4) 
#     W <- Matrix::sparseMatrix(i = adj[,1], j = adj[,2], x = 1)
#     M.inv <- Matrix::Diagonal(n = nrow(W), x = rowSums(W))
#     W.plus <- solve(M.inv, W)
#   }
# 
#   ### theta.hat
#   theta.hat <- colMeans(mcmc.res$samples[-c(0:burnin),])
# 
#   # # split parameters
#   # theta.hat.mu <- theta.hat[1:5]
#   # theta.hat.cov <- theta.hat[-c(1:5)]
# 
#   # # mean
#   # mu.hat <- phenom.mean(theta.hat.mu, D = D.mat, x = X.2011, u = u, v = v, vrsn = mean.vrsn)  
#   # # precision matrix
#   # Omega.hat <- phenom.prec(theta.hat.cov, n = n.y, W.p = W.plus, M.inv = Minv, type = cov.type)
#   
#   ll.hat <- phenomLikelihood(theta.hat, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, 2)
# 
#   D.bar <- -2 * mean(mcmc.res$loglik[-c(0:burnin)])
#   D.hat <- -2 * ll.hat
# 
#   # return DIC
#   2 * D.bar - D.hat
# }
# 
# beta.prior.ll <- function(beta, s2.beta = 100){
#   # Beta prior density evaluation.
#   # Input: 
#   #   beta: vector of parameter values
#   #   s2.beta: variance hyperparameter (default = 100)
#   # Output: 
#   #   Beta prior density evaluation.
#   sum(dnorm(beta, mean = 0, sd = sqrt(s2.beta), log = TRUE))
# }
# 
# alpha.prior.ll <- function(alpha, s2.alpha = 100){
#   # Alpha prior density evaluation.
#   # Input: 
#   #   alpha: parameter value
#   #   s2.alpha: variance hyperparameter (default = 100)
#   # Output: 
#   #   Alpha prior density evaluation.
#   dnorm(alpha, mean = 0, sd = sqrt(s2.alpha), log = TRUE) + log(2)
# }
# 
# s2.mu.prior.ll <- function(s2.mu, s2 = 250){
#   # s2.mu prior density evaluation.
#   # Input: 
#   #   s2.mu: parameter value
#   #   s2.beta: variance hyperparameter (default = 250)
#   # Output: 
#   #   s2.mu prior density evaluation.
#   dnorm(s2.mu, mean = 0, sd = sqrt(s2), log = TRUE) + log(2)
# }
# 
# phenom.mcmc <- function(usa, n.mcmc, theta.0, cov.type = "iid", update = 100, print.iter = 10){
#   # Generate MCMC samples from the phenomenological model.
#   # Input: 
#   #   usa: data raster
#   #   n.mcmc: number of MCMC samples
#   #   theta.0: initial parameter values
#   #   cov.type: type of covariance ("iid", "CAR", "SAR")
#   #   update: when to update adapative tuning proposal distribution
#   #   print.iter: value specifying how often to print updates
#   # Output: 
#   #   Generates a list, with mcmc samples ("samples"), log-likelihood evaluations ("loglik"), 
#   #     and adaptive tuning information ("tuning").
# 
#   ### Set-Up
# 
#   # X and y values
#   y.2011 <- values(usa$so4$so4.2011)
#   X.2011 <- usa$X$X.2011
#   n.y <- length(y.2011)
# 
#   # wind
#   u <- values(usa$wind.small$wind.2011$uwind)
#   v <- values(usa$wind.small$wind.2011$vwind)
# 
#   # distance
#   r.coords <- coordinates(usa$so4$so4.2011)
#   D.mat <- distm(r.coords, fun = distHaversine) / 1000  
# 
#   # covariance structure
#   if (cov.type == "iid"){
#     # assume independent variance
#     W.plus = NA; M.inv = NA
#   } else {
#     # weights matrix for autoregressive covariances
#     adj <- raster::adjacent(usa$so4$so4.2011, cells = 1:length(y.2011), pairs = TRUE, directions = 4) 
#     W <- Matrix::sparseMatrix(i = adj[,1], j = adj[,2], x = 1)
#     M.inv <- Matrix::Diagonal(n = nrow(W), x = rowSums(W))
#     W.plus <- solve(M.inv, W)
#   }
# 
#   # mcmc structures
#   theta.mat <- matrix(NA_real_, nrow = n.mcmc, ncol = length(theta.0))
#   ll.vals <- rep(NA_real_, length = n.mcmc)
#   theta.k <- theta.0 # initialize theta
#   mu.inds <- 1:5; alpha.ind <- 1; s2.mu.ind <- 2; beta.inds <- 3:5; s2.var.ind <- 6
#   if (cov.type != "iid"){
#     rho.ind <- 7
#   }
# 
#   ### adaptive tuning
# 
#   ## default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234 
# 
#   ## initial tuning parameters
#   adapt <- list()
#   adapt$time <- 1
# 
#   # beta
#   adapt$beta <- list()
#   adapt$beta$S0 <- diag(c(0.1, 0.00001, 0.00001))
#   adapt$beta$var <- 0.05
#   adapt$beta$n.jumps <- 0
#   adapt$beta$S <- adapt$beta$var * adapt$beta$S0
#   
#   # alpha
#   adapt$alpha <- list()
#   adapt$alpha$s0 <- 1
#   adapt$alpha$var <- .01
#   adapt$alpha$n.jumps <- 0
#   adapt$alpha$sd <- sqrt(adapt$alpha$var * adapt$alpha$s0)
#   
#   # s2.mu
#   adapt$s2.mu <- list()
#   adapt$s2.mu$s0 <- 1
#   adapt$s2.mu$var <- 0.5
#   adapt$s2.mu$n.jumps <- 0
#   adapt$s2.mu$sd <- sqrt(adapt$s2.mu$var * adapt$s2.mu$s0)
# 
#   if (cov.type != "iid"){
#     # rho
#     adapt$rho <- list()
#     adapt$rho$s0 <- 1
#     adapt$rho$var <- .01
#     adapt$rho$n.jumps <- 0
#     adapt$rho$sd <- sqrt(adapt$rho$var * adapt$rho$s0)
#   }
# 
#   ### MCMC sampler
#   for (k in 1:n.mcmc){
# 
#     ## a. sample s2.nu
# 
#     # mean
#     mu.k <- phenom.mean(theta.k[mu.inds], D = D.mat, x = X.2011, u = u, v = v, vrsn = 2)
#     # precision matrix
#     Omega.k <- phenom.prec(theta.k[-mu.inds], n = n.y, W.p = W.plus, M.inv = M.inv, type = cov.type)
# 
#     resid.k <- (y.2011 - mu.k)
#     ig.a <- n.y / 2
#     ig.b <- as.numeric((t(resid.k) %*% (Omega.k %*% resid.k) / 2) * theta.k[s2.var.ind])
#     s2.var <- 1 / rgamma(n = 1, shape = ig.a, rate = ig.b)
# 
#     # replace structures with correct marginal variance parameter
#     theta.k[s2.var.ind] <- s2.var
#     ll.k <- phenomLikelihood(theta.k, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, 2)
# 
#     ## b. sample betas
# 
#     # propose new values
#     theta.star <- theta.k
#     theta.star[beta.inds] <- mvrnorm(n = 1, mu = theta.k[beta.inds], Sigma = adapt$beta$S)
#     ll.star <- phenomLikelihood(theta.star, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, 2)
# 
#     # compare MH ratio (remember: symmetric proposal distribution)
#     log.mh.beta <- (ll.star + beta.prior.ll(theta.star[beta.inds])) - 
#       (ll.k + beta.prior.ll(theta.k[beta.inds]))
#     
#     if (log(runif(1)) < log.mh.beta){
#       # accept betas
#       theta.k <- theta.star
#       ll.k <- ll.star
#       adapt$beta$n.jumps <- adapt$beta$n.jumps + 1
#     }
# 
#     ## c. sample alpha
#     theta.star <- theta.k
#     theta.star[alpha.ind] <- rnorm(n = 1, mean = theta.k[alpha.ind], sd = adapt$alpha$sd)
#     ll.star <- phenomLikelihood(theta.star, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, 2)
# 
#     # compare MH ratio (remember: symmetric proposal distribution)
#     log.mh.alpha <- (ll.star + alpha.prior.ll(theta.star[alpha.ind])) - 
#       (ll.k + alpha.prior.ll(theta.k[alpha.ind]))
# 
#     if (log(runif(1)) < log.mh.alpha){
#       # accept alpha
#       theta.k <- theta.star
#       ll.k <- ll.star
#       adapt$alpha$n.jumps <- adapt$alpha$n.jumps + 1
#     }
# 
#     ## d. sample s2.mu
# 
#     theta.star <- theta.k
#     theta.star[s2.mu.ind] <- rnorm(n = 1, mean = theta.k[s2.mu.ind], sd = adapt$s2.mu$sd)
#     ll.star <- phenomLikelihood(theta.star, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, 2)
# 
#     # compare MH ratio (remember: symmetric proposal distribution)  
#     log.mh.s2.mu <- (ll.star + s2.mu.prior.ll(theta.star[s2.mu.ind])) - 
#       (ll.k + s2.mu.prior.ll(theta.k[s2.mu.ind]))
# 
#     if (log(runif(1)) < log.mh.s2.mu){
#       # accept s2.mu
#       theta.k <- theta.star
#       ll.k <- ll.star
#       adapt$s2.mu$n.jumps <- adapt$s2.mu$n.jumps + 1
#     }
# 
# 
#     if (cov.type != "iid"){
#       ## e. sample rho
#       theta.star <- theta.k
#       theta.star[rho.ind] <- rnorm(n = 1, mean = theta.k[rho.ind], sd = adapt$rho$sd)
#       ll.star <- phenomLikelihood(theta.star, D.mat, y.2011, X.2011, u, v, W.plus, M.inv, cov.type, 2)
# 
#       # compare MH ratio (uniform prior (0,1), symmetric proposal distribution)  
#       log.mh.rho <- (ll.star - ll.k)
# 
#       if (log(runif(1)) < log.mh.rho){
#         # accept rho
#         theta.k <- theta.star
#         ll.k <- ll.star
#         adapt$rho$n.jumps <- adapt$rho$n.jumps + 1
#       }
#     }
# 
#     # save theta.k
#     theta.mat[k,] <- theta.k
#     ll.vals[k] <- ll.k
# 
#     ### update adaptive structures
#     if (k %% update == 0){
# 
#       # samples used for update
#       samples <- theta.mat[(adapt$time - 1) * update + 1:update, ]
#       
#       # update constants
#       gamma1 <- 1 / adapt$time^c1; gamma2 <- c0 * gamma1
#       
#       ## a. update beta adaptive tuning
#       S.hat <- var(samples[, beta.inds])
#       r.hat <- adapt$beta$n.jumps / update
#       log.var <- log(adapt$beta$var) + gamma2 * (r.hat - r.opt)
#       adapt$beta$var <- exp(log.var)
#       adapt$beta$S0 <- adapt$beta$S0 +
#         gamma1 * (S.hat - adapt$beta$S0)
#       adapt$beta$n.jumps <- 0
#       adapt$beta$S <- adapt$beta$var * adapt$beta$S0
#       
#       ## b. update alpha adaptive tuning
#       s.hat <- var(samples[, alpha.ind])
#       r.hat <- adapt$alpha$n.jumps / update
#       log.var <- log(adapt$alpha$var) + gamma2 * (r.hat - r.opt)
#       adapt$alpha$var <- exp(log.var)
#       adapt$alpha$s0 <- adapt$alpha$s0 +
#         gamma1 * (s.hat - adapt$alpha$s0)
#       adapt$alpha$n.jumps <- 0
#       adapt$alpha$sd <- sqrt(adapt$alpha$var * adapt$alpha$s0)
#       
#       ## c. update s2.mu adaptive tuning
#       s.hat <- var(samples[, s2.mu.ind])
#       r.hat <- adapt$s2.mu$n.jumps / update
#       log.var <- log(adapt$s2.mu$var) + gamma2 * (r.hat - r.opt)
#       adapt$s2.mu$var <- exp(log.var)
#       adapt$s2.mu$s0 <- adapt$s2.mu$s0 +
#         gamma1 * (s.hat - adapt$s2.mu$s0)
#       adapt$s2.mu$n.jumps <- 0
#       adapt$s2.mu$sd <- sqrt(adapt$s2.mu$var * adapt$s2.mu$s0)
#       
#       if (cov.type != "iid"){
#         ## d. update rho adaptive tuning
#         s.hat <- var(samples[, rho.ind])
#         r.hat <- adapt$rho$n.jumps / update
#         log.var <- log(adapt$rho$var) + gamma2 * (r.hat - r.opt)
#         adapt$rho$var <- exp(log.var)
#         adapt$rho$s0 <- adapt$rho$s0 +
#           gamma1 * (s.hat - adapt$rho$s0)
#         adapt$rho$n.jumps <- 0
#         adapt$rho$sd <- sqrt(adapt$rho$var * adapt$rho$s0)
#       }
# 
#       # update number of adaptive steps
#       adapt$time <- adapt$time + 1
#     }
# 
#     if (k %% print.iter == 0){
#       cat(paste("Iteration ", k, ",  ", 100 * k / n.mcmc, "% Complete.", "\n", sep = ""))
#       if (cov.type == "iid"){
#         cat("  Parameters: alpha, s2.mu, beta.0, beta.1, beta.2, s2.cov\n")
#         cat(paste(round(theta.k, 3))); cat("\n")
#       } else {
#         cat("  Parameters: alpha, s2.mu, beta.0, beta.1, beta.2, s2.cov, rho\n")
#         cat(paste(round(theta.k, 3))); cat("\n")
#       }
#       cat("--------------------------------------------------------------------------------------------\n")
#     }
#   }
#   ## return ll.vals (for DIC) and theta.mat (mcmc samples)
# 
#   if (cov.type == "iid"){
#     colnames(theta.mat) <- c("alpha", "s2.mu", "beta.0", "beta.1", "beta.2", "s2.cov")
#   } else {
#     colnames(theta.mat) <- c("alpha", "s2.mu", "beta.0", "beta.1", "beta.2", "s2.cov", "rho")
#   }
# 
#   return(list(samples = theta.mat, loglik = ll.vals, tuning = adapt))
# }
# 
# plotPhenom <- function(usa, theta, cov.type, mean.vrsn = 2, plot.type = "mean", ...){
#   # Create plots to assess fitted phenomenological model.
#   # Input: 
#   #   usa: data raster
#   #   theta: vector of parameter values
#   #   cov.type: type of covariance ("iid", "CAR", "SAR")
#   #   mean.vrsn: type of phenom. mean structure 
#   #     (1 = constant variance, 2 = variance scaled by power plant output)
#   #   plot.type: type of plot (either "mean" or "sim")
#   # Output: 
#   #   Creates a raster plot of the mean surface, or of a simulation from the 
#   #     phenomenological model.
# 
#   ## set-up
# 
#   # X and y values
#   y.2011 <- values(usa$so4$so4.2011)
#   X.2011 <- usa$X$X.2011
# 
#   # wind
#   u <- values(usa$wind.small$wind.2011$uwind)
#   v <- values(usa$wind.small$wind.2011$vwind)
# 
#   # distance
#   r.coords <- coordinates(usa$so4$so4.2011)
#   D.mat <- distm(r.coords, fun = distHaversine) / 1000  
# 
#   if (cov.type == "iid"){
#     # assume independent variance
#     W.plus = NA; M.inv = NA
#   } else {
#     # weights matrix for autoregressive covariances
#     adj <- raster::adjacent(usa$so4$so4.2011, cells = 1:length(y.2011), pairs = TRUE, directions = 4) 
#     W <- Matrix::sparseMatrix(i = adj[,1], j = adj[,2], x = 1)
#     M.inv <- Matrix::Diagonal(n = nrow(W), x = rowSums(W))
#     W.plus <- solve(M.inv, W)
#   }
# 
#   # split parameters
#   theta.mu <- theta[1:5]
#   theta.cov <- theta[-c(1:5)]
# 
#   # mean
#   mu <- phenom.mean(theta.mu, D = D.mat, x = X.2011, u = u, v = v, vrsn = mean.vrsn)  
# 
#   if (plot.type == "mean"){
# 
#     my.r <- usa$so4$so4.2011
#     values(my.r) <- mu
# 
#   } else {
#     # precision matrix
#     n.y <- length(y.2011)
#     Omega <- phenom.prec(theta.cov, n = n.y, W.p = W.plus, M.inv = M.inv, type = cov.type)
# 
#     my.r <- usa$so4$so4.2011
# 
#     # sample from full model
#     L <- chol(Omega)
#     n.y <- nrow(Omega)
#     values(my.r) <- mu + as.vector(solve(L, rnorm(n.y)))
#   }
# 
#   return(my.r)
# }
# 
# ################################################################################
# ### 5. Functions for fitting models with colored (Matern) driving noise
# ################################################################################
# 
# ll.inla <- function(theta, y, X, mats, m.reorder = list(FALSE)){
#   # Calculate the log-likelihood given sampled parameters for coupled 
#   #   SO2-SO4 time-averaged model (WITH WIND!).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the coupled SO2-SO4 model, with wind.
#   
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
#   
#   # create matrices
#   D <- mats$D
#   C <- mats$C
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#   
#   # constanst for log-likelihood calc.
#   Q <- Matrix::crossprod(A.1) / s2
#   mu <- rep(0, nrow(Q))
#   sources <- Matrix::solve(A.2, X * beta)
#   b <- Matrix::crossprod(sources, A.1) / s2
# 
#   if (m.reorder[[1]]){
# 		# log-likelihood
#   	ll <- inla.qsample(Q = Q,
#                      b = b,
#                      mu = mu,
#                      sample = y,
#                      logdens = TRUE,
#                      compute.mean = FALSE,
# 										 reordering = m.reorder[[2]])
# 	} else {
# 	  # log-likelihood
#   	ll <- inla.qsample(Q = Q,
#     	                 b = b,
#       	               mu = mu,
#         	             sample = y,
#           	           logdens = TRUE,
#             	         compute.mean = FALSE)
# 	}
# 
#   # return ll
#   return(ll$logdens)
# }
# 
# ll.data.Matern.W <- function(theta, y, X, mats){
#   # Calculate the log-likelihood given sampled parameters for coupled 
#   #   SO2-SO4 time-averaged model (WITH WIND!).
#   # Input:
#   #   theta: sampled parameters
#   #   y: observed data
#   #   X: vector of SOZ emissions
#   #   mats: matrix operators representing diffusion, advection, etc.
#   # Output:
#   #   Log-likelihood of the coupled SO2-SO4 model, with wind.
# 
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   alpha <- theta[4]
#   sigma <- theta[5]
#   rho <- theta[6]
#   delta <- theta[7]
# 
#   # transform rho and s2 to kappa and tau
#   kappa <- sqrt(8) / rho
#   tau <- sigma * kappa * sqrt(4 * pi)
# 
#   # create matrices
#   D <- mats$D
#   C <- mats$C
#   L <- mats$L
#   n.cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#   Q.half <- (L + Diagonal(n.cells, kappa^2)) / tau
#   Q <- Matrix::crossprod(Q.half)
# 
#   # constant for log-likelihood calc.
#   A1.y <- A.1 %*% y
#   sources <- solve(A.2, X * beta)
#   Q.sources <- Q %*% sources
#   sources.cp <- Matrix::crossprod(sources, Q.sources) 
#   log.det.A1 <- as.numeric(determinant(A.1)$modulus)
#   log.det.Q <- as.numeric(determinant(Q)$modulus)
#   log.det.S <- 2 * log.det.A1 + log.det.Q
#   
#   # log-likelihood
#   ll <- (-n.cells * log(2 * pi) / 2) + log.det.S / 2-
#     1 / 2 * (Matrix::crossprod(A1.y, (Q %*% A1.y))  - 2 * Matrix::crossprod(A1.y, Q.sources) + sources.cp)
# 
#   # return ll
#   return(as.numeric(ll))
# }
# 
# 
# pc.prior <- function(sigma, rho, hyper, log = TRUE){
#   # Penalized complexity prior for GRFs from Fuglstad et al (2019).
#   # Input:
#   #   sigma: sd of marginal variance for Matern GRF, nu = 1
#   #   rho: spatial range parameter
#   #   hyper: list with hyperparamters (alpha1, alpha2, rho0, sigma0)
#   #   log: Boolean, return log-density or density
#   # Output:
#   #   Prior density for Matern GRF (nu = 1) for (sigma, rho).
#   
#   lambda.1 <- -log(hyper$alpha1) * hyper$rho0
#   lambda.2 <- -log(hyper$alpha2) / hyper$sigma0
#   ll <- log(lambda.1) + log(lambda.2) - 2 * log(rho) - (lambda.1 / rho) - (lambda.2 * sigma)
# 
#   if (log){
#     return(ll)
#   } else {
#     return(exp(ll))
#   }
# }
# 
# coupledMCMC.Matern.W <- function(N, y, X, mats, params, priors, update, iter){
#   # MCMC sampler for the coupled SO2-SO4 model (WITH WIND AND Matern noise process!).
#   # Input:
#   #   N: number of MCMC samples.
#   #   y: observed data.
#   #   X: vector of SO2 outputs.
#   #   mats: matrix operators approximating physical processes.
#   #   params: vector of initial parameters.
#   #   priors: hyperparameters for prior distribution.
#   #   update: number of iterations at which to update MCMC proposal dist
#   #     (uses Shaby and Wells (2009) algorithm).
#   #   iter: when to print out iterations.
#   # Output:
#   #   List with MCMC samples, DIC calc, and adaptive structures used in MCMC 
#   # proposals.
#   
#   # number of observations
#   n.cells <- length(y)
# 
#   # initial parameter values
#   gamma.k <- params$gamma
#   xi.k <- params$xi
#   beta.k <- params$beta
#   alpha.k <- params$alpha
#   sigma.k <- params$sigma
#   rho.k <- params$rho
#   delta.k <- params$delta
#   
#   # all parameters in a vector
#   theta.k <- c(gamma.k, xi.k, beta.k, alpha.k, sigma.k, rho.k, delta.k)
#   
#   # matrix to hold parameters
#   samples <- matrix(NA_real_, nrow = N, ncol = length(theta.k))
#   colnames(samples) <- c("gamma", "xi", "beta", "alpha", "sigma", "rho", "delta")
#   
#   # vector to hold correction factor samples for DIC calculation
#   dic.vec <- rep(NA_real_, length = N)
#   
#   # calculate log-likelihood for current values
#   ll.curr <- ll.data.Matern.W(theta.k, y, X, mats)
#   
#   ### create adaptive structures
#   
#   # matrix to hold samples
#   update.samples <- matrix(NA_real_, nrow = update, ncol = length(theta.k))
#   
#   # default constants
#   c0 <- 1; c1 <- 0.8; r.opt <- 0.234; time <- 1
#   
#   # structure for gamma
#   ad.gamma <- list()
#   ad.gamma$s0 <- 1
#   ad.gamma$var <- 10
#   ad.gamma$n.jumps <- 0
#   ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#   
#   # structure for xi
#   ad.xi <- list()
#   ad.xi$s0 <- 1
#   ad.xi$var <- 1
#   ad.xi$n.jumps <- 0
#   ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#   
#   # structure for alpha
#   ad.alpha <- list()
#   ad.alpha$s0 <- 1
#   ad.alpha$var <- 1
#   ad.alpha$n.jumps <- 0
#   ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)
# 
#   # structure for sigma
#   ad.sigma <- list()
#   ad.sigma$s0 <- 1
#   ad.sigma$var <- 10
#   ad.sigma$n.jumps <- 0
#   ad.sigma$sd <- sqrt(ad.sigma$var * ad.sigma$s0)
# 
#   # structure for rho
#   ad.rho <- list()
#   ad.rho$s0 <- 1
#   ad.rho$var <- 10
#   ad.rho$n.jumps <- 0
#   ad.rho$sd <- sqrt(ad.rho$var * ad.rho$s0)
#   
#   ### generate mcmc samples
#   for (k in 1:N){
#     
#     # a. sample gamma
#     gamma.star <- rnorm(n = 1, mean = gamma.k, sd = ad.gamma$sd)
#     
#     if (gamma.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[1] <- gamma.star
#       ll.star <- ll.data.Matern.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(gamma.star, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(gamma.k, mean = 0, sd = priors$g.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new gamma
#         gamma.k <- gamma.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.gamma$n.jumps <- ad.gamma$n.jumps + 1
#       }
#     }
#     
#     # b. sample xi
#     xi.star <- rnorm(n = 1, mean = xi.k, sd = ad.xi$sd)
#     
#     if (xi.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[2] <- xi.star
#       ll.star <- ll.data.Matern.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dexp(xi.star, rate = priors$xi, log = TRUE)
#       mh.denom <- ll.curr + dexp(xi.k, rate = priors$xi, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new xi
#         xi.k <- xi.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.xi$n.jumps <- ad.xi$n.jumps + 1
#       }
#     }
#     
#     # sample alpha
#     alpha.star <- rnorm(n = 1, mean = alpha.k, sd = ad.alpha$sd)
#     
#     if (alpha.star > 0){
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[4] <- alpha.star
#       ll.star <- ll.data.Matern.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + dnorm(alpha.star, mean = 0, sd = priors$alpha.sd, log = TRUE)
#       mh.denom <- ll.curr + dnorm(alpha.k, mean = 0, sd = priors$alpha.sd, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new xi
#         alpha.k <- alpha.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.alpha$n.jumps <- ad.alpha$n.jumps + 1
#       }
#     }
#     
#     # c. sample sigma and rho
#     sigma.star <- rnorm(n = 1, mean = sigma.k, sd = ad.sigma$sd)
#     rho.star <- rnorm(n = 1, mean = rho.k, sd = ad.rho$sd)
#     if (sigma.star > 0 & rho.star > 0){
#       
#       # calculate new log-likelihood
#       theta.star <- theta.k
#       theta.star[5] <- sigma.star
#       theta.star[6] <- rho.star
#       ll.star <- ll.data.Matern.W(theta.star, y, X, mats)
#       
#       # determine log-MH ratio
#       mh.num <- ll.star + pc.prior(sigma.star, rho.star, hyper = priors$pc, log = TRUE)
#       mh.denom <- ll.curr + pc.prior(sigma.k, rho.k, hyper = priors$pc, log = TRUE)
#       mh.ratio <- mh.num - mh.denom
#       
#       if (log(runif(1)) < mh.ratio){
#         # accept new sigma and rho
#         sigma.k <- sigma.star
#         rho.k <- rho.star
#         theta.k <- theta.star
#         ll.curr <- ll.star
#         ad.sigma$n.jumps <- ad.sigma$n.jumps + 1
#         ad.rho$n.jumps <- ad.rho$n.jumps + 1
#       }
#     }
#     
#     # d. sample beta
#     D <- mats$D
#     C <- mats$C
#     L <- mats$L
#     n.cells <- nrow(D)
#     kappa.k <- sqrt(8) / rho.k
#     tau.k <- sigma.k * kappa.k * sqrt(4 * pi)
# 
#     A.1 <- gamma.k * D + alpha.k * C + Diagonal(n.cells, delta.k)
#     A.2 <- (gamma.k / xi.k) * D + (alpha.k / xi.k) * C + Diagonal(n.cells, 1)
#     Q.half <- L + Diagonal(n.cells, kappa.k) / tau.k
#     Q <- Matrix::crossprod(Q.half)
# 
#     Z.k <- solve(A.2, X)
#     Q.Z.k <- Q %*% Z.k
#     A1.y <- A.1 %*% y
#     beta.var <- 1 / (as.numeric(Matrix::crossprod(Z.k, Q.Z.k)) + (1 / priors$beta))
#     beta.mean <- beta.var * as.numeric(Matrix::crossprod(A1.y, Q.Z.k)) 
#     
#     beta.k <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = beta.mean, sd = sqrt(beta.var))
#     
#     # calculate new likelihood
#     theta.k[3] <- beta.k
#     ll.curr <- ll.data.Matern.W(theta.k, y, X, mats)
#     
#     # save ll for DIC calc
#     dic.vec[k] <- -2 * ll.curr
#     
#     # d. save theta.k
#     samples[k,] <- theta.k
#     
#     # e. update adaptive structures
#     
#     if (k %% update > 0){
#       # add theta.k to sample matrix
#       update.samples[k %% update, ] <- theta.k
#     } else {
#       # add theta.k to sample matrix
#       update.samples[update, ] <- theta.k
#       
#       # update constants
#       gamma1 <- 1 / time^c1; gamma2 <- c0 * gamma1
#       
#       ### update gamma adaptive structure
#       s.hat <- var(update.samples[,1])
#       r.hat <- ad.gamma$n.jumps / update
#       log.var <- log(ad.gamma$var) + gamma2 * (r.hat - r.opt)
#       ad.gamma$var <- exp(log.var)
#       ad.gamma$s0 <- ad.gamma$s0 +
#         gamma1 * (s.hat - ad.gamma$s0)
#       ad.gamma$n.jumps <- 0
#       ad.gamma$sd <- sqrt(ad.gamma$var * ad.gamma$s0)
#       
#       ### update xi adaptive structure
#       s.hat <- var(update.samples[,2])
#       r.hat <- ad.xi$n.jumps / update
#       log.var <- log(ad.xi$var) + gamma2 * (r.hat - r.opt)
#       ad.xi$var <- exp(log.var)
#       ad.xi$s0 <- ad.xi$s0 +
#         gamma1 * (s.hat - ad.xi$s0)
#       ad.xi$n.jumps <- 0
#       ad.xi$sd <- sqrt(ad.xi$var * ad.xi$s0)
#       
#       ### update alpha adaptive structure
#       s.hat <- var(update.samples[,4])
#       r.hat <- ad.alpha$n.jumps / update
#       log.var <- log(ad.alpha$var) + gamma2 * (r.hat - r.opt)
#       ad.alpha$var <- exp(log.var)
#       ad.alpha$s0 <- ad.alpha$s0 +
#         gamma1 * (s.hat - ad.alpha$s0)
#       ad.alpha$n.jumps <- 0
#       ad.alpha$sd <- sqrt(ad.alpha$var * ad.alpha$s0)
#       
#       ### update sigma adaptive structure
#       s.hat <- var(update.samples[,5])
#       r.hat <- ad.sigma$n.jumps / update
#       log.var <- log(ad.sigma$var) + gamma2 * (r.hat - r.opt)
#       ad.sigma$var <- exp(log.var)
#       ad.sigma$s0 <- ad.sigma$s0 +
#         gamma1 * (s.hat - ad.sigma$s0)
#       ad.sigma$n.jumps <- 0
#       ad.sigma$sd <- sqrt(ad.sigma$var * ad.sigma$s0)
# 
#       ### update rho adaptive structure
#       s.hat <- var(update.samples[,6])
#       r.hat <- ad.rho$n.jumps / update
#       log.var <- log(ad.rho$var) + gamma2 * (r.hat - r.opt)
#       ad.rho$var <- exp(log.var)
#       ad.rho$s0 <- ad.rho$s0 +
#         gamma1 * (s.hat - ad.rho$s0)
#       ad.rho$n.jumps <- 0
#       ad.rho$sd <- sqrt(ad.rho$var * ad.rho$s0)
# 
#       # update number of adaptive steps
#       time <- time + 1
#     }
#     
#     # f. print out the current iteration
#     if (iter > 0){
#       if (k %% iter == 0){
#         print(paste("iteration:", k))
#       }
#     }
#   }
#   
#   # save adaptive structures
#   ad.str <- list(gamma = ad.gamma,
#                  xi = ad.xi,
#                  alpha = ad.alpha,
#                  sigma = ad.sigma,
#                  rho = ad.rho)
#   
#   # return samples
#   results <- list(samples = samples, dic = dic.vec, ad.str = ad.str)
#   return(results)
# }
# 
# 
# matern.dic <- function(mcmc.res, burnin, y, X, mats){
#   # Calculate the DIC for the Matern noise model.
#   # Input:
#   #   mcmc.res: list with samples from MCMC
#   #   burnin: vector of burnin indices
#   #   y: response
#   #   X: covariates (e.g., power plant emissions data)
#   #   mats: list with matrices; used for numerical approximation
#   # Output: 
#   #   DIC for Matern noise model.
# 
#   mean.deviance <- mean(mcmc.res$dic[-burnin])
#   mean.theta <- colMeans(mcmc.res$samples[-burnin,])
#   ll.theta <- ll.data.Matern.W(mean.theta, y, X, mats)
#   deviance.theta.hat <- -2 * ll.theta
#   penalty <- mean.deviance - deviance.theta.hat
#   dic <- penalty + mean.deviance
#   return(dic)
# }
# 
# cov.sample <- function(theta, y, X, mats, matern = TRUE){
#   # Sample spatial surface with given covariance structure.
#   # Input:
#   #   theta: vector of parameters
#   #   y: response
#   #   X: covariates (e.g., power plant emissions data)
#   #   mats: list with matrices; used for numerical approximation
#   #   matern: boolean, matern or white noise
#   # Output: 
#   #   Sampled spatial surface.
# 
#   if (matern){
#     # grab parameters
#     gamma <- theta[1]
#     xi <- theta[2]
#     beta <- theta[3]
#     alpha <- theta[4]
#     sigma <- theta[5]
#     rho <- theta[6]
#     delta <- theta[7]
# 
#     # transform rho and s2 to kappa and tau
#     kappa <- sqrt(8) / rho
#     tau <- sigma * kappa * sqrt(4 * pi)
# 
#     # create matrices
#     D <- mats$D
#     C <- mats$C
#     L <- mats$L
#     n.cells <- nrow(D)
#     A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#     A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#     Q.half <- (L + Diagonal(n.cells, kappa^2)) / tau
#     Q <- Matrix::crossprod(Q.half)
#     S.inv <- t(A.1) %*% (Q %*% A.1)
# 
#   } else {
#     # grab parameters
#     gamma <- theta[1]
#     xi <- theta[2]
#     beta <- theta[3]
#     alpha <- theta[6]
#     s2 <- theta[4]
#     delta <- theta[5]
# 
#     # create matrices
#     D <- mats$D
#     C <- mats$C
#     L <- mats$L
#     n.cells <- nrow(D)
#     A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
# 
#     S.inv <- (t(A.1) %*% A.1) / s2
#   }
#   # constant for log-likelihood calc.
#   # A1.y <- A.1 %*% y
#   # Q.sources <- Q %*% sources
#   # sources.cp <- Matrix::crossprod(sources, Q.sources) 
#   # log.det.A1 <- as.numeric(determinant(A.1)$modulus)
#   # log.det.Q <- as.numeric(determinant(Q)$modulus)
#   # log.det.S <- 2 * log.det.A1 + log.det.Q
# 
#   inla.qsample(n = 1, Q = S.inv, b = rep(0, n.cells),
#                 mu = rep(0, n.cells))
# }
# 
# mean.calc <- function(theta, y, X, mats, matern = TRUE){
#   # Determine mean surface from given parameters.
#   # Input:
#   #   theta: vector of parameters
#   #   y: response
#   #   X: covariates (e.g., power plant emissions data)
#   #   mats: list with matrices; used for numerical approximation
#   #   matern: boolean, matern or white noise
#   # Output: 
#   #   Mean SO4 surface.
# 
#   if (matern){
#     # grab parameters
#     gamma <- theta[1]
#     xi <- theta[2]
#     beta <- theta[3]
#     alpha <- theta[4]
#     sigma <- theta[5]
#     rho <- theta[6]
#     delta <- theta[7]
# 
#     # transform rho and s2 to kappa and tau
#     kappa <- sqrt(8) / rho
#     tau <- sigma * kappa * sqrt(4 * pi)
# 
#     # create matrices
#     D <- mats$D
#     C <- mats$C
#     L <- mats$L
#     n.cells <- nrow(D)
#     A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#     A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#     # Q.half <- (L + Diagonal(n.cells, kappa^2)) / tau
#     # Q <- Matrix::crossprod(Q.half)
#   } else {
#     # grab parameters
#     gamma <- theta[1]
#     xi <- theta[2]
#     beta <- theta[3]
#     alpha <- theta[6]
#     s2 <- theta[4]
#     delta <- theta[5]
# 
#     # create matrices
#     D <- mats$D
#     C <- mats$C
#     L <- mats$L
#     n.cells <- nrow(D)
#     A.1 <- gamma * D + alpha * C + Diagonal(n.cells, delta)
#     A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n.cells, 1)
#   }
#   # constant for log-likelihood calc.
#   # A1.y <- A.1 %*% y
#   sources <- solve(A.2, X * beta)
#   mu <- solve(A.1, sources)
#   # Q.sources <- Q %*% sources
#   # sources.cp <- Matrix::crossprod(sources, Q.sources) 
#   # log.det.A1 <- as.numeric(determinant(A.1)$modulus)
#   # log.det.Q <- as.numeric(determinant(Q)$modulus)
#   # log.det.S <- 2 * log.det.A1 + log.det.Q
#   
#   return(list(so4 = mu, so2 = sources))
# }
# 
# 
# gCalc <- function(theta, mats, X_mat, P_mat, s_vec, key){
# 
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
#   b0 <- theta[7]
# 
#   D <- mats$D
#   C <- mats$C
# 
#   # create matrices
#   n_cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n_cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n_cells, 1)
#   # calculate mean (no intercept)
#   sources <- solve(A.2, X_mat * beta)
#   mu <- solve(A.1, sources)
# 
#   # determine source-receptor matris (T)
#   T_mat <- Matrix::crossprod(P_mat, mu)
# 
#   # calculate G^*
#   G_star <- T_mat %*% s_vec - T_mat[cbind(1:nrow(T_mat), key)] * s_vec[key]
# 
#   # calculate G from G^*
#   G <- (G_star / max(G_star))[,1]
# 
#   return(G)
# }
# 
# upwindCovar <- function(theta, mats, X_mat, P_mat, s_vec, key, covars){
# 
#   # grab parameters
#   gamma <- theta[1]
#   xi <- theta[2]
#   beta <- theta[3]
#   s2 <- theta[4]
#   delta <- theta[5]
#   alpha <- theta[6]
#   b0 <- theta[7]
# 
#   D <- mats$D
#   C <- mats$C
# 
#   # create matrices
#   n_cells <- nrow(D)
#   A.1 <- gamma * D + alpha * C + Diagonal(n_cells, delta)
#   A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n_cells, 1)
#   # calculate mean (no intercept)
#   sources <- solve(A.2, X_mat * beta)
#   mu <- solve(A.1, sources)
# 
#   # determine source-receptor matris (T)
#   T_mat <- Matrix::crossprod(P_mat, mu)
# 
#   # calculate upwind sum
#   upwind_covars <- T_mat %*% covars - T_mat[cbind(1:nrow(T_mat), key)] * covars[key, ]
# 
#   # standardize to [0,1]
#   upwind_max <- apply(upwind_covars, 2, max)
#   upwind_final <- upwind_covars / upwind_max[col(upwind_covars)] 
# 
#   return(upwind_final)
# }
# 
# 
# 
# createMexicoX <- function(r, emissions){
#   # Create the X design matrix for the pollution data.
#   # Input:
#   #   r: a single raster object
#   #   emissions: a data frame with power plant emissions data.
#   # Output:
#   #   A design vector, X, representing the total amount of emissions in 
#   #     each grid of the raster.
#   
#   # number of cells in rectangle
#   n.cells <- length(values(r))
#   
#   # create an SO2 variable using the emissions data
#   em.vals <- rep(0, n.cells)
#   em.inds <- cellFromXY(r, cbind(emissions$longitude, emissions$latitude))
#   
#   # make sure unique inds are
#   unique.inds <- unique(em.inds)
#   
#   for(curr.ind in unique.inds){
#     spec.points <- em.inds == curr.ind
#     em.vals[curr.ind] <- sum(emissions$total.SO2[spec.points]) / 12
#   }
#   
#   #em.vals[em.inds] <- emissions$totSO2emissions
#   
#   # combine all values into a single vector
#   X <- em.vals
#   
#   # return X
#   return(X)
# }
