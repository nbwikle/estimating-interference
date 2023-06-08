### sensitivity-ntrees.R
### Nathan Wikle

###############################################################################
### 1. Sensitivity study set-up
###############################################################################

### parallel set-up

library(parallel) # note: 'parallel' only works on POSIX systems!
workers <- detectCores()

### load data

# outcome data
out.df <- readRDS(here::here("data", "outcome-data.RDS"))
out.df <- out.df %>% drop_units() 
out.df <- na.omit(out.df) 
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
# air pollution transport data
so4.mcmc <- readRDS(here::here("output", "tx-2016-100k-int-mx.RDS"))
# copula of theta posterior values used in the analysis
theta_post <- readRDS(here::here("output", "copula-posterior.RDS"))


### create covariate matrix and z, g values

# covariate matrix
X.sim <- data.frame(out.df %>%
  st_drop_geometry() %>%
  dplyr::select(c(
    "log_pop", "black_prop", "log_income", "log_density", "smokerate", 
    "prcp", "log_optime", 
  )))

# Z vector
z.vec <- out.df %>% dplyr::pull(scrubbed)

# G_bar
burnin <- 25000
theta.mu <- colMeans(so4.mcmc$samples[-c(1:burnin), ])

# G.bar
g.mean <- gCalc(
  theta = theta.mu, 
  mats = adv.mats, 
  X_mat = em.mat, 
  P_mat = proj.mat, 
  s_vec = scrubber.vec, 
  key = key.assoc, 
  version = 3
)


###############################################################################
### 2. Perform sensitivity study
###############################################################################

### i. m = 100 trees
m100 <- bartSensitivity(
  n_chains = 50,
  n_tree = 100,
  base_val = 0.95,
  power_val = 2,
  theta_posterior = theta_post,
  data_list = list(
    proj_mat = proj.mat,
    em_mat = em.mat,
    adv_mats = adv.mats,
    s_vec = scrubber.vec,
    key_assoc = key.assoc,
    x = X.sim,
    z = z.vec,
    g_mean = g.mean
  ),
  n_workers = workers
)

saveRDS(m100, here::here("output", "sensitivity-m100.RDS"))

### ii. m = 200 trees
m200 <- bartSensitivity(
  n_chains = 50,
  n_tree = 200,
  base_val = 0.95,
  power_val = 2,
  theta_posterior = theta_post,
  data_list = list(
    proj_mat = proj.mat,
    em_mat = em.mat,
    adv_mats = adv.mats,
    s_vec = scrubber.vec,
    key_assoc = key.assoc,
    x = X.sim,
    z = z.vec,
    g_mean = g.mean
  ),
  n_workers = workers
)

saveRDS(m200, here::here("output", "sensitivity-m200.RDS"))

### iii. m = 300 trees
m300 <- bartSensitivity(
  n_chains = 50,
  n_tree = 300,
  base_val = 0.95,
  power_val = 2,
  theta_posterior = theta_post,
  data_list = list(
    proj_mat = proj.mat,
    em_mat = em.mat,
    adv_mats = adv.mats,
    s_vec = scrubber.vec,
    key_assoc = key.assoc,
    x = X.sim,
    z = z.vec,
    g_mean = g.mean
  ),
  n_workers = workers
)

saveRDS(m300, here::here("output", "sensitivity-m300.RDS"))

### iv. m = 400 trees
m400 <- bartSensitivity(
  n_chains = 50,
  n_tree = 400,
  base_val = 0.95, 
  power_val = 2,
  theta_posterior = theta_post,
  data_list = list(
    proj_mat = proj.mat,
    em_mat = em.mat,
    adv_mats = adv.mats,
    s_vec = scrubber.vec,
    key_assoc = key.assoc,
    x = X.sim,
    z = z.vec,
    g_mean = g.mean
  ),
  n_workers = workers
)

saveRDS(m400, here::here("output", "sensitivity-m400.RDS"))



