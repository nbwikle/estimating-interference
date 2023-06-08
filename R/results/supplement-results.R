### supplement-results.R
### Nathan Wikle


####################################################################
### 1. HyADS vs OU model comparison
####################################################################

### i. data

# SO4 and emissions data
tx <- readRDS(here::here("data", "texas-2016-data.RDS"))

# FVM matrices
adv.mats <- readRDS(here::here("data", "advec-diff-mats.RDS"))

# gridded HyADS 2016 PM2.5 totals
grid.dat <- read.fst(
  here::here("data", "HyADS_grids_pm25_total_2016.fst"), 
  as.data.table = TRUE
)
grid.dat[ , year := 2016]

# fitted OU model
so4.results <- readRDS(here::here("output", "tx-2016-100k-int-mx.RDS"))


### ii. Compare estimated total SO4 due to power plants for both models

# OU Model Output

# posterior mean - sulfate model parameters
burnin <- 25000
theta.mu <- colMeans(so4.results$samples[-c(1:burnin), ])

# make all mexico emissions = 0
mx.vec <- tx$X$X.usa
mx.vec[mx.vec > 0] <- 0

# save surface as a raster
mu.hat <- meanSurface(
  theta = theta.mu,
  X_tx = tx$X$X.usa,
  X_mx = mx.vec,
  mats = adv.mats
)

# save results as raster
so4.ou <- tx$so4
values(so4.ou) <- mu.hat


# HyADS Model Output

# coordinate reference system projection string for spatial data
p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# dcast to get year columns
grid.dat.c <- data.table::dcast(grid.dat, x + y ~ year, value.var = 'vals.out')
# convert to raster
grid.dat.r <- raster::rasterFromXYZ(grid.dat.c, crs = p4s)

# get USA dataset
states48 <- datasets::state.name[!(state.name %in% c( 'Alaska', 'Hawaii'))]
usa.sub.sf <- USAboundaries::us_states() %>%
  st_transform(st_crs(p4s))
usa.sub.sf <- usa.sub.sf[usa.sub.sf$name %in% states48,]

# convert to correct raster
hyads.2016 <- raster::projectRaster(from = grid.dat.r$X2016, crs = crs(so4.ou))
so4.hyads <- raster::crop(extend(hyads.2016, so4.ou), so4.ou)


# Plot estimated surfaces
plotMeanSO4Comp(
  plot.name = here::here("output", "supplement", "sFIG1-total_so4.png"),
  ou.raster = so4.ou, 
  hyads.raster = so4.hyads
)


### iii. Compare relative contribution of power plant facility

# grab hyads unit-level contributions
hyads.units <- read.fst(
  here::here("data", "HyADS_grids_pm25_byunit_2016.fst"),
  as.data.table = TRUE
)
  
# OU Model: Facility 3497 ('Big Brown Power Plant')

# emissions totals
em.counterfactual <- tx$X$X.usa
em.counterfactual[-7126] <- 0

# save surface as a raster
mu.u3497 <- meanSurface(
  theta = theta.mu,
  X_tx = em.counterfactual,
  X_mx = mx.vec,
  mats = adv.mats
)

# calculate proportion of SO4 attributed to Facility 3497
prop.u3497.ou <- tx$so4
values(prop.u3497.ou) <- mu.u3497 / mu.hat

# HyADS: Facility 3497 ('Big Brown Power Plant')

# facility IDs (HyADS)
hyads.names <- colnames(hyads.units)[-c(1:2)]
hy.idx1 <- sub('.', '', unlist(lapply(
  base::strsplit(as.character(unique(hyads.names)), ".", TRUE), function(x) {
    x[1]
})))

# combine totals across units
fac.i <- which(hy.idx1 == "3497") + 2
hyads.u3497 <- hyads.units[ , .SD, .SDcols = c(1,2,fac.i)]
hyads.u3497[, vals.out := rowSums(.SD), .SDcols = 3:(2 + length(fac.i))]

# proportion of SO4 attributed to Facility 3497
hyads.prop <- hyads.u3497$vals.out / grid.dat$vals.out
hyads.prop.r <- grid.dat.r$X2016
values(hyads.prop.r) <- hyads.prop
# crop raster
prop.u3497.hyads <- raster::projectRaster(from = hyads.prop.r, crs = crs(so4.ou))
prop.u3497.hyads <- raster::crop(extend(prop.u3497.hyads, so4.ou), so4.ou)

# Plot proportion of SO4
plotPropSO4Comp(
  plot.name = here::here("output", "supplement", "sFIG2-prop_so4_u3497.png"),
  ou.raster = prop.u3497.ou,
  hyads.raster = prop.u3497.hyads
)


####################################################################
### 2. Covariate Balance and overlap
####################################################################

### i. data

# MCMC results
so4.mcmc <- readRDS(here::here("output", "tx-2016-100k-int-mx.RDS")) 
# outcome data
outcome.data <- readRDS(here::here("data", "outcome-data.RDS"))
out.df <- na.omit(outcome.data) 
# add black carbon data to outcome data
bc.data <- readRDS(here::here("data", "bc-2016-zcta.RDS"))
out.df$bc <- bc.data[,1]
# projection from raster to zip
proj.mat <- readRDS(here::here("data", "projection-mat.RDS")) 
# emissions matrix
em.mat <- readRDS(here::here("data", "emissions-mat.RDS")) 
# key-associated facilities
key.assoc <- readRDS(here::here("data", "key-assoc-vec.RDS"))
# scrubber status
scrubber.vec <- readRDS(here::here("data", "scrubber-vec.RDS"))
# advection-diffusion matrices
adv.mats <- readRDS(here::here("data", "advec-diff-mats.RDS"))
# facility-level data
fac.df <- readRDS(here::here("data", "facility-data-2016.RDS")) 

# add degree

burnin <- 25000
theta.mu <- colMeans(so4.mcmc$samples[-c(1:burnin), ])

# G.bar
G.bar <- gCalc(
  theta = theta.mu, 
  mats = adv.mats, 
  X_mat = em.mat, 
  P_mat = proj.mat, 
  s_vec = scrubber.vec, 
  key = key.assoc, 
  version = 3
)
out.df$G <- G.bar

# degree
upwind.covars <- upwindCovar(
  theta = theta.mu,
  mats = adv.mats,
  X_mat = em.mat,
  P_mat = proj.mat,
  s_vec = scrubber.vec,
  key = key.assoc,
  covars = fac.df[ ,"totHeatInput"],
  version = 2
)
out.df$degree <- upwind.covars[[1]]

### ii. covariate balance with Z

# restrict to Z and covariates
z_covars <- out.df %>%
  mutate(Z = scrubbed) %>%
  dplyr::select(
    "log_pop", "female", "ped_prop", "median_age", "white_prop",
    "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income",
    "move_rate", "ins_prop", "renter_housing", "urban_prop",
    "log_density", "tmin", "tmax", "prcp", "vp", "rel_humid",
    "smokerate", "log_heat", "log_optime", "pct_capacity",
    "dist_to_key", "degree", "Z"
  ) %>%
  st_drop_geometry() %>% drop_units()

# BART propensity scores
bart.fit <- dbarts::bart(
  x.train = z_covars[,-27], y.train = z_covars[,27], nthread = 1, 
  ntree = 200, ndpost = 500, keepevery = 5, nskip = 5000, keeptrees = TRUE
)
pred.vals <- predict(bart.fit, newdata = z_covars[,-27])
prob.z <- colMeans(pred.vals)

# BART weights
z.vals <- z_covars[,27]
bart.weights <- (z.vals) / prob.z + (1 - z.vals) / (1 - prob.z)

# CBPS propensity scores
cbps.est <- CBPS(Z ~ ., data = z_covars, ATT = 0)
prob.z2 <- cbps.est$fitted.values

# CBPS weights
cbps.weights <- (z.vals) / prob.z2 + (1 - z.vals) / (1 - prob.z2)
cbps.weights[which(cbps.weights > 20)] <- cbps.weights[1110]


# calculate absolute mean difference in treatment
diff.mat1 <- calc.balance(
  cov.matrix = z_covars[, -27],
  trt1 = which(z.vals == 1),
  trt2 = which(z.vals == 0),
  weights = bart.weights, d.type = 2
)

diff.mat2 <- calc.balance(
  cov.matrix = z_covars[, -27],
  trt1 = which(z.vals == 1),
  trt2 = which(z.vals == 0),
  weights = cbps.weights, d.type = 2
)

# combine into a single difference matrix
diff.mat <- cbind(diff.mat1, diff.mat2[,2])
colnames(diff.mat) <- c("Unweighted", "BART", "CBPS")

# plot Z balance results
plotZBalance(
  plot.name = here::here("output", "supplement", "sFIG3-z_balance.png"),
  diff.mat = diff.mat
)

# plot propensity score overlap (for Z)
plotZOverlap(
  plot.name = here::here("output", "supplement", "sFIG4-ps_overlap.png"),
  z = z.vals, prob.z = prob.z2
)

### iii. covariate balance with G

# restrict to G and relevant covariates
g_covars <- out.df %>%
  mutate(Z = scrubbed) %>%
  dplyr::select(
    "log_pop", "female", "ped_prop", "median_age", "white_prop",
    "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income",
    "move_rate", "ins_prop", "renter_housing", "urban_prop",
    "log_density", "tmin", "tmax", "prcp", "vp", "rel_humid",
    "smokerate", "log_heat", "log_optime", "pct_capacity",
    "dist_to_key", "degree", "G"
  ) %>%
  st_drop_geometry() %>% drop_units()

# unweighted correlations
unweighted.cor <- bal.tab(G ~ ., data = g_covars)$Balance[,2]

# CBPS-weighted correlations
w.cbps <- weightit(G ~ ., data = g_covars, method = "cbps", estimand = "ATE")
w.cbps <- WeightIt::trim(w.cbps, at = 10)
cbps.cor <- bal.tab(w.cbps)$Balance[,3]

# BART-weighted correlations
w.bart <- WeightIt::weightit(G ~ ., data = g_covars, method = "bart", estimand = "ATE")
w.bart <- WeightIt::trim(w.bart, at = 10)
bart.cor <- bal.tab(w.bart)$Balance[,3]

# combine correlations into single matrix
corr.mat <- cbind(
  abs(unweighted.cor[-27]),
  abs(bart.cor),
  abs(cbps.cor)
)
colnames(corr.mat) <- c("unweighted", "bart", "cbps")
rownames(corr.mat) <- rownames(diff.mat)

# plot G balance results
plotGBalance(
  plot.name = here::here("output", "supplement", "sFIG5-g_balance.png"),
  corr.mat = corr.mat
)

# plot generalized propensity score overlap

# colors and transparency
cols <- c("#434a60", "#c68f57", "#749e9b", "#89533e")
trans <- 1

# CBPS fit
propensity.g <- CBPS(G ~ ., data = g_covars)

# quantiles of g.bar
g.cutoff <- quantile(g_covars$G, probs = c(0.25, 0.5, 0.75))
g.q1 <- g_covars$G < g.cutoff[1]
g.q2 <- (g_covars$G >= g.cutoff[1]) & (g_covars$G < g.cutoff[2])
g.q3 <- (g_covars$G >= g.cutoff[2]) & (g_covars$G < g.cutoff[3])
g.q4 <- g_covars$G >= g.cutoff[3]


# plot G overlap

plotGOverlap(
  plot.name = here::here("output", "supplement", "sFIG6-gps_overlap.png"),
  g.cutoff = g.cutoff,
  g.propensity = propensity.g,
  g.quantiles = list(g.q1, g.q2, g.q3, g.q4),
  cols = cols, trans = trans
)



####################################################################
### 3. Stratified causal effect estimates
####################################################################


### i. Medicare analysis

# Medicare data

# log-linear BART
bart.plugin <- readRDS(here::here("output", "medicare-bart-plugin.RDS"))
bart.cut <- readRDS(here::here("output", "medicare-bart-cut.RDS"))

# Plot DE stratified by degree and distance to key-associated facility

# degree limits
y.lims.de.deg <- c(
  min(
    min(1000 * apply(bart.plugin$DE$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$DE$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

# dist limits
y.lims.de.dist <- c(
  min(
    min(1000 * apply(bart.plugin$DE$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$d4, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d4, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$DE$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$d4, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d4, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

y.lims.de.strat <- c(
  min(c(y.lims.de.deg, y.lims.de.dist)), 
  max(c(y.lims.de.deg, y.lims.de.dist))
)

plotStratEffects(
  plot.name = here::here("output", "supplement", "sFIG7a-medicare_bart_de_strat.png"),
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  y.lims = y.lims.de.strat,
  cols = c("#D3B1C2","#613659", "#E1C391", "#705446"),
  outcome.type = 1,
  ce.type = "DE"
)


# Plot IE0 stratified by degree and distance to key-associated facility

# degree limits
y.lims.ie0.deg <- c(
  min(
    min(1000 * apply(bart.plugin$IE0$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE0$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

# dist limits
y.lims.ie0.dist <- c(
  min(
    min(1000 * apply(bart.plugin$IE0$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$d4, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d4, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE0$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$d4, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d4, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

y.lims.ie0.strat <- c(
  min(c(y.lims.ie0.deg, y.lims.ie0.dist)), 
  max(c(y.lims.ie0.deg, y.lims.ie0.dist))
)

plotStratEffects(
  plot.name = here::here("output", "supplement", "sFIG7b-medicare_bart_ie0_strat.png"),
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  y.lims = y.lims.ie0.strat,
  cols = c("#D3B1C2","#613659", "#E1C391", "#705446"),
  outcome.type = 1,
  ce.type = "IE0"
)


# Plot IE1 stratified by degree and distance to key-associated facility

# degree limits
y.lims.ie1.deg <- c(
  min(
    min(1000 * apply(bart.plugin$IE1$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE1$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

# dist limits
y.lims.ie1.dist <- c(
  min(
    min(1000 * apply(bart.plugin$IE1$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$d4, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d4, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE1$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$d4, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d4, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

y.lims.ie1.strat <- c(
  min(c(y.lims.ie1.deg, y.lims.ie1.dist)), 
  max(c(y.lims.ie1.deg, y.lims.ie1.dist))
)

plotStratEffects(
  plot.name = here::here("output", "supplement", "sFIG7c-medicare_bart_ie1_strat.png"),
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  y.lims = y.lims.ie1.strat,
  cols = c("#D3B1C2","#613659", "#E1C391", "#705446"),
  outcome.type = 1,
  ce.type = "IE1"
)

# remove objects
rm(bart.plugin)
rm(bart.cut)


### ii. Asthma analysis

# Asthma data

# log-linear BART
bart.plugin <- readRDS(here::here("output", "asthma-bart-plugin.RDS"))
bart.cut <- readRDS(here::here("output", "asthma-bart-cut.RDS"))


# Plot DE stratified by degree and distance to key-associated facility

# degree limits
y.lims.de.deg <- c(
  min(
    min(1000 * apply(bart.plugin$DE$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$DE$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

# dist limits
y.lims.de.dist <- c(
  min(
    min(1000 * apply(bart.plugin$DE$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$DE$d4, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$DE$d4, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$DE$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$DE$d4, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$DE$d4, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

y.lims.de.strat <- c(
  min(c(y.lims.de.deg, y.lims.de.dist)), 
  max(c(y.lims.de.deg, y.lims.de.dist))
)

plotStratEffects(
  plot.name = here::here("output", "supplement", "sFIG8a-asthma_bart_de_strat.png"),
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  y.lims = y.lims.de.strat,
  cols = c("#D3B1C2","#613659", "#E1C391", "#705446"),
  outcome.type = 2,
  ce.type = "DE"
)


# Plot IE0 stratified by degree and distance to key-associated facility

# degree limits
y.lims.ie0.deg <- c(
  min(
    min(1000 * apply(bart.plugin$IE0$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE0$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

# dist limits
y.lims.ie0.dist <- c(
  min(
    min(1000 * apply(bart.plugin$IE0$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE0$d4, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE0$d4, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE0$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE0$d4, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE0$d4, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

y.lims.ie0.strat <- c(
  min(c(y.lims.ie0.deg, y.lims.ie0.dist)), 
  max(c(y.lims.ie0.deg, y.lims.ie0.dist))
)

plotStratEffects(
  plot.name = here::here("output", "supplement", "sFIG8b-asthma_bart_ie0_strat.png"),
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  y.lims = y.lims.ie0.strat,
  cols = c("#D3B1C2","#613659", "#E1C391", "#705446"),
  outcome.type = 2,
  ce.type = "IE0"
)


# Plot IE1 stratified by degree and distance to key-associated facility

# degree limits
y.lims.ie1.deg <- c(
  min(
    min(1000 * apply(bart.plugin$IE1$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$ldeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$mdeg, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$hdeg, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE1$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$ldeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$mdeg, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$hdeg, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

# dist limits
y.lims.ie1.dist <- c(
  min(
    min(1000 * apply(bart.plugin$IE1$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.plugin$IE1$d4, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d1, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d2, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d3, 2, quantile, probs = 0.025, na.rm = TRUE)),
    min(1000 * apply(bart.cut$IE1$d4, 2, quantile, probs = 0.025, na.rm = TRUE))
  ),
  max(
    max(1000 * apply(bart.plugin$IE1$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.plugin$IE1$d4, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d1, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d2, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d3, 2, quantile, probs = 0.975, na.rm = TRUE)),
    max(1000 * apply(bart.cut$IE1$d4, 2, quantile, probs = 0.975, na.rm = TRUE))
  )
)

y.lims.ie1.strat <- c(
  min(c(y.lims.ie1.deg, y.lims.ie1.dist)), 
  max(c(y.lims.ie1.deg, y.lims.ie1.dist))
)

plotStratEffects(
  plot.name = here::here("output", "supplement", "sFIG8c-asthma_bart_ie1_strat.png"),
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  y.lims = y.lims.ie1.strat,
  cols = c("#D3B1C2","#613659", "#E1C391", "#705446"),
  outcome.type = 2,
  ce.type = "IE1"
)

# remove objects
rm(bart.plugin)
rm(bart.cut)


####################################################################
### 3. Simulation Study: Major Findings
####################################################################


### i. Coverage and Variance - PM1

# load results
ss.results <- readRDS(here::here("output", "simstudy-results.RDS"))

# coverage results

# Poisson outcome model
p.plugin <- ss.results$coverage$PM1$pois$plugin[2:7,]
p.cut <- ss.results$coverage$PM1$pois$cut[2:7,]
# BART outcome model
b.plugin <- ss.results$coverage$PM1$bart$plugin[2:7,]
b.cut <- ss.results$coverage$PM1$bart$cut[2:7,]

# proportion of variance

# Poisson outcome model
p.var <- matrix(0, nrow = 6, ncol = 3)
p.var[,1] <- ss.results$variance$PM1$pois$DE[2:7,2] / 
  ss.results$variance$PM1$pois$DE[2:7,3]
p.var[,2] <- ss.results$variance$PM1$pois$IE0[2:7,2] / 
  ss.results$variance$PM1$pois$IE0[2:7,3]
p.var[,3] <- ss.results$variance$PM1$pois$IE1[2:7,2] / 
  ss.results$variance$PM1$pois$IE1[2:7,3]
# BART outcome model
b.var <- matrix(0, nrow = 6, ncol = 3)
b.var[,1] <- ss.results$variance$PM1$bart$DE[2:7,2] / 
  ss.results$variance$PM1$bart$DE[2:7,3]
b.var[,2] <- ss.results$variance$PM1$bart$IE0[2:7,2] / 
  ss.results$variance$PM1$bart$IE0[2:7,3]
b.var[,3] <- ss.results$variance$PM1$bart$IE1[2:7,2] / 
  ss.results$variance$PM1$bart$IE1[2:7,3]

# plot results (Figure 10 in Supplement)
plotPM1Results(
  plot.name = here::here("output", "supplement", "sFIG10-coverage_pm1.png"),
  p.plugin = p.plugin,
  p.cut = p.cut,
  b.plugin = b.plugin,
  b.cut = b.cut,
  p.var = p.var,
  b.var = b.var
)



### ii. Coverage and Variance - PM3

# coverage results

# Poisson outcome model
p.plugin <- ss.results$coverage$PM3$pois$plugin[2:7,]
p.cut <- ss.results$coverage$PM3$pois$cut[2:7,]
# BART outcome model
b.plugin <- ss.results$coverage$PM3$bart$plugin[2:7,]
b.cut <- ss.results$coverage$PM3$bart$cut[2:7,]

# proportion of variance

# Poisson outcome model
p.var <- matrix(0, nrow = 6, ncol = 3)
p.var[,1] <- ss.results$variance$PM3$pois$DE[2:7,2] / 
  ss.results$variance$PM3$pois$DE[2:7,3]
p.var[,2] <- ss.results$variance$PM3$pois$IE0[2:7,2] / 
  ss.results$variance$PM3$pois$IE0[2:7,3]
p.var[,3] <- ss.results$variance$PM3$pois$IE1[2:7,2] / 
  ss.results$variance$PM3$pois$IE1[2:7,3]
# BART outcome model
b.var <- matrix(0, nrow = 6, ncol = 3)
b.var[,1] <- ss.results$variance$PM3$bart$DE[2:7,2] / 
  ss.results$variance$PM3$bart$DE[2:7,3]
b.var[,2] <- ss.results$variance$PM3$bart$IE0[2:7,2] / 
  ss.results$variance$PM3$bart$IE0[2:7,3]
b.var[,3] <- ss.results$variance$PM3$bart$IE1[2:7,2] / 
  ss.results$variance$PM3$bart$IE1[2:7,3]

# plot results (Figure 10 in Supplement)
plotPM3Results(
  plot.name = here::here("output", "supplement", "sFIG11-coverage_pm3.png"),
  p.plugin = p.plugin,
  p.cut = p.cut,
  b.plugin = b.plugin,
  b.cut = b.cut,
  p.var = p.var,
  b.var = b.var
)
  
### iii. Compare change in coverage with variance results

# data frame showing change in coverage and % variance for each simulation study
delta.comps <- deltaCoverage(ss.results)

# plot change in coverage vs % variance
plotDeltaComps(
  plot.name = here::here("output", "supplement", "sFIG12-delta_coverage"),
  delta.comps = delta.comps
)



