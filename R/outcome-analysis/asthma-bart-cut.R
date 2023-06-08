### asthma-bart-cut.R
### Nathan Wikle

### 1. set up cluster

library(parallel) # note: 'parallel' only works on POSIX systems!

# number of multiple imputation samples
n.sims <- 250
# number of cores to run in parallel
workers <- detectCores()

### 2. data structures used in the analysis

# air pollution transport data
so4.mcmc <- readRDS(here::here("output", "tx-2016-100k-int-mx.RDS")) 
# asthma data
out.df <- readRDS(here::here("data", "outcome-data.RDS"))
out.df <- out.df %>% drop_units() 
out.df <- na.omit(out.df) 
# add black carbon data to outcome data
bc.data <- readRDS(here::here("data", "bc-2016-zcta.RDS"))
out.df$bc <- bc.data[,1]
# facility-level data
fac.df <- readRDS(here::here("data", "facility-data-2016.RDS")) 
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


### 3. Fit Poisson regression with cut feedback

# sample from pi(theta | sulfate data)
set.seed(329)
burnin <- 25000
so4.samples <- so4.mcmc$samples[-c(1:burnin), ]
theta.samples <- so4.samples[sample(nrow(so4.samples), size = n.sims), ]

# estimate DE(g), IE(g) in parallel
results.cut <- mclapply(
  1:n.sims, 
  function(x) guncBART(
    sim.k = x, n_burn = 5000, n_post = 250, n_thin = 5
  ),
  mc.cores = workers
)

### 4. Save results
saveRDS(results.cut, here::here("output", "asthma-bart-cut_synth.RDS"))








