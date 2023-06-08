### CM2-BART.R
### Nathan Wikle


###############################################################################
### 1. Simulation study set-up
###############################################################################

### parallel set-up

# run in parallel, if possible
library(parallel) # note: 'parallel' only works on POSIX systems!
workers <- detectCores()

# number of outcome simulations 
n_sims <- 100
# number of multiple imputation draws
n_gsamples <- 100


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
### 2. Perform simulation study
###############################################################################

ss.results.BART <- mclapply(
  1:n_sims, 
  function(x) CM_SimStudy_BART(
    ss_num = x,
    n_gsamples = n_gsamples,
    data_list = list(
      so4_mats = adv.mats,
      em_mat = em.mat,
      proj_mat = proj.mat,
      s_vec = scrubber.vec,
      key_assoc = key.assoc,
      x = X.sim,
      z = z.vec,
      g_mean = g.mean
    ),
    theta_posterior = theta_post, 
    model = 2
  ),
  mc.cores = workers
)

###############################################################################
### 3. Save results
###############################################################################

# save result
saveRDS(ss.results.BART, file = here::here("output", "simstudy-cm2-bart.RDS"))


