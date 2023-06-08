### asthma-bart-plugin.R
### Nathan Wikle

########################################################################
# 1. Load libraries and data
########################################################################

# air pollution transport data
so4.mcmc <- readRDS(here::here("output", "tx-2016-100k-int-mx.RDS")) # MCMC results

# asthma data
out.df <- readRDS(here::here("data", "outcome-data.RDS"))
out.df <- out.df %>% drop_units() 

# remove NA rows (also removes asthma rate, as it is not used later)
out.df <- na.omit(out.df) # 142 zcta's are dropped because of NAs

# add black carbon data to outcome data
bc.data <- readRDS(here::here("data", "bc-2016-zcta.RDS"))
out.df$bc <- bc.data[,1]

# facility-level data
fac.df <- readRDS(here::here("data", "facility-data-2016.RDS")) 


########################################################################
### 2. Estimate G
########################################################################

### A. data structures used in the analysis

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


### B. Calculate G.bar, upwind heat input, and degree from E(theta | sulfate)

# g value
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

# degree and upwind heat input
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
out.df$log_up_heat <- upwind.covars[[2]]

# remove ZCTAs with 0 children...
zero.peds <- which(out.df$Ped_16 == 0)
out.final <- out.df[-zero.peds,]


########################################################################
# 3. Fit log-linear BART regression
########################################################################

# covariate matrix
X.mat <- as.matrix(out.final %>%
  st_drop_geometry() %>%
  mutate(Z = scrubbed, Z.G = scrubbed * G) %>%
  dplyr::select(c(
    "Ped_16", "female", "ped_prop", "median_age", "white_prop",
    "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income",
    "move_rate", "ins_prop", "renter_housing", "urban_prop",
    "log_density", "smokerate", "tmin", "tmax", "prcp", "vp", "rel_humid", "bc",
    "log_heat", "log_optime", "pct_capacity", "dist_to_key", "log_up_heat", "degree", 
    "Z", "G", "Z.G"
)) %>% drop_units())

# response
y <- out.final %>% dplyr::pull(asthma)

# offset
offset.vec <- log(out.final %>% dplyr::pull(Ped_16))

# determine leaf hyperparameters
rates <- y / exp(offset.vec)
r_star <- quantile(rates, probs = 0.975)
r_mean <- mean(rates)
a_0 <- 0.5 * (log(r_star) - log(r_mean))


# BART MCMC samples
bart.nounc <- count_bart(
  y = y,
  x = X.mat,
  offset = exp(offset.vec + log(r_mean)),
  nburn = 10000,
  nsim = 1000,
  nthin = 5,
  update_interval = 1000,
  ntree = 200,
  a0 = a_0,
  sd = 2 * sd(y),
  base = 0.95,
  power = 2,
  nu = 3,
  lambda = NULL,
  sigq = .9,
  sighat = NULL, # prior specification for continuous y... only used as 'ghost' inputs
  kappa_a = 5, kappa_b = 3, # shape parameters for kappa prior (beta prime distribution)
  count_model = "poisson", # must be either 'poisson', 'nb', 'zipoisson', or 'zinb'
  debug = FALSE,
  # used as input... but not relevant for count models...?
  randeff_design = matrix(1),
  randeff_variance_component_design = matrix(1),
  randeff_scales = 1,
  randeff_df = 3
)


########################################################################
# 6. Estimate causal effects
########################################################################

bart.plugin <- bartEstimates(
  fit = bart.nounc, n_post = 1000, r_mean = r_mean, 
  X = X.mat[, -c(ncol(X.mat) - 1, ncol(X.mat))], 
  degree_idx = which(colnames(X.mat) == "degree"), 
  dist_idx = which(colnames(X.mat) == "dist_to_key"),
  up_heat_idx = which(colnames(X.mat) == "log_up_heat"),
  heat_idx = which(colnames(X.mat) == "log_heat")
)

########################################################################
# 7. Save results
########################################################################
saveRDS(bart.plugin, here::here("output", "asthma-bart-plugin_synth.RDS"))

