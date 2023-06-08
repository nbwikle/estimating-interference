### asthma-pois-plugin.R
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
# 3. Fit Poisson regression
########################################################################

# separate into X and y matrices

# covariate matrix
X.mat <- as.matrix(out.final %>%
  st_drop_geometry() %>%
  mutate(Z = scrubbed, Z.G = scrubbed * G) %>%
  dplyr::select(c(
    "female", "ped_prop", "median_age", "white_prop",
    "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income",
    "move_rate", "ins_prop", "renter_housing", "urban_prop",
    "log_density", "smokerate", "tmin", "tmax", "prcp", "vp", "rel_humid", "bc",
    "log_heat", "log_optime", "pct_capacity", "dist_to_key", "log_up_heat", "degree", 
    "Z", "G", "Z.G"
)) %>% drop_units())

# response vector
y <- out.final %>% dplyr::pull(asthma)

# offset vector
offset.vec <- out.final %>% dplyr::pull(Ped_16)

# MCMC parameterization
n.mcmc = 750000
n.burn = 250000
n.thin = 10

# number of posterior samples to retain
n.post <- 10000

# GLM estimate
glm.mean <- stats::glm(
  y ~ . + offset(log(denom)) - denom,
  family = stats::poisson, 
  data = data.frame(y, denom = offset.vec, X.mat)
)

# Poisson regression 
pois.res <- pois_reg_cpp(
  y_vector = y,
  X_matrix = X.mat,
  offset = log(offset.vec),
  theta_0 = glm.mean$coefficients,
  n_mcmc = n.mcmc,
  thin = n.thin,
  burnin = n.burn,
  n_adapt = 5000,
  Sigma_0 = tuningMatFull()[-23,-23], 
  keep_burnin = FALSE,
  prior_var = 1
)

# acceptance rate
pois.res[[2]]


########################################################################
# 6. Estimate causal effects
########################################################################

pois.plugin <- poisEstimates(
  fit = pois.res,
  n_post = n.post,
  X = X.mat[, -c((ncol(X.mat) - 2):ncol(X.mat))],
  degree_idx = which(colnames(X.mat) == "degree"),
  dist_idx = which(colnames(X.mat) == "dist_to_key"),
  up_heat_idx = which(colnames(X.mat) == "log_up_heat"),
  heat_idx = which(colnames(X.mat) == "log_heat")
)


########################################################################
# 7. Save results
########################################################################
saveRDS(pois.plugin, here::here("output", "asthma-pois-plugin_synth.RDS"))


