### sensitivity-results.R
### Nathan Wikle
###
### Processes the sensitivity study results.


#####################################################################
### 1. Load relevant data, set-up sensitivity study objects
#####################################################################

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

# G_mean
burnin <- 25000
theta.mu <- colMeans(so4.mcmc$samples[-c(1:burnin), ])

g.mean <- gCalc(
  theta = theta.mu, 
  mats = adv.mats, 
  X_mat = em.mat, 
  P_mat = proj.mat, 
  s_vec = scrubber.vec, 
  key = key.assoc, 
  version = 3
)

#####################################################################
### 2. Number of Trees
#####################################################################

### load results
m100 <- readRDS(here::here("output", "sensitivity-m100.RDS"))
m200 <- readRDS(here::here("output", "sensitivity-m200.RDS"))
m300 <- readRDS(here::here("output", "sensitivity-m300.RDS"))
m400 <- readRDS(here::here("output", "sensitivity-m400.RDS"))

### marginal standard deviation

# estimates
mtree.sd <- mtreeStatCalc(
  m100 = m100, 
  m200 = m200,
  m300 = m300,
  m400 = m400,
  X = X.sim,
  res.type = "sd"
)

# plot results
plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG31a-sd_DE.png"),
  mtree = mtree.sd,
  res.type = "sd",
  ce.type = "DE"
)

plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG31b-sd_IE0.png"),
  mtree = mtree.sd,
  res.type = "sd",
  ce.type = "IE0"
)

plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG31c-sd_IE1.png"),
  mtree = mtree.sd,
  res.type = "sd",
  ce.type = "IE1"
)

### bias

# estimates
mtree.bias <- mtreeStatCalc(
  m100 = m100, 
  m200 = m200,
  m300 = m300,
  m400 = m400,
  X = X.sim,
  res.type = "bias"
)

# plot results
plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG32a-bias_DE.png"),
  mtree = mtree.bias,
  res.type = "bias",
  ce.type = "DE"
)

plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG32b-bias_IE0.png"),
  mtree = mtree.bias,
  res.type = "bias",
  ce.type = "IE0"
)

plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG32c-bias_IE1.png"),
  mtree = mtree.bias,
  res.type = "bias",
  ce.type = "IE1"
)

### rmse

# estimates
mtree.rmse <- mtreeStatCalc(
  m100 = m100, 
  m200 = m200,
  m300 = m300,
  m400 = m400,
  X = X.sim,
  res.type = "rmse"
)

# plot results
plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG33a-rmse_DE.png"),
  mtree = mtree.rmse,
  res.type = "rmse",
  ce.type = "DE"
)

plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG33b-rmse_IE0.png"),
  mtree = mtree.rmse,
  res.type = "rmse",
  ce.type = "IE0"
)

plotNTreesRes(
  plot.name = here::here("output", "supplement", "sFIG33c-rmse_IE1.png"),
  mtree = mtree.rmse,
  res.type = "rmse",
  ce.type = "IE1"
)



#####################################################################
### 3. Tree Depth
#####################################################################

### tree depth results
pow05 <- readRDS(here::here("output", "sensitivity-power05.RDS"))
pow1 <- readRDS(here::here("output", "sensitivity-power1.RDS"))
pow15 <- readRDS(here::here("output", "sensitivity-power15.RDS"))
pow2 <- readRDS(here::here("output", "sensitivity-power2.RDS"))
pow25 <- readRDS(here::here("output", "sensitivity-power25.RDS"))
pow3 <- readRDS(here::here("output", "sensitivity-power3.RDS"))
pow4 <- readRDS(here::here("output", "sensitivity-power4.RDS"))
pow5 <- readRDS(here::here("output", "sensitivity-power5.RDS"))


### marginal standard deviation

# estimates
power.sd <- depthStatCalc(
  p05 = pow05, 
  p1 = pow1,
  p15 = pow15,
  p2 = pow2,
  p25 = pow25,
  p3 = pow3,
  p4 = pow4,
  p5 = pow5,
  X = X.sim,
  res.type = "sd"
)

# plot results
plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG36a-sd_DE.png"),
  power = power.sd,
  res.type = "sd",
  ce.type = "DE"
)

plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG36b-sd_IE0.png"),
  power = power.sd,
  res.type = "sd",
  ce.type = "IE0"
)

plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG36c-sd_IE1.png"),
  power = power.sd,
  res.type = "sd",
  ce.type = "IE1"
)


### bias

# estimates
power.bias <- depthStatCalc(
  p05 = pow05, 
  p1 = pow1,
  p15 = pow15,
  p2 = pow2,
  p25 = pow25,
  p3 = pow3,
  p4 = pow4,
  p5 = pow5,
  X = X.sim,
  res.type = "bias"
)

# plot results
plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG37a-bias_DE.png"),
  power = power.bias,
  res.type = "bias",
  ce.type = "DE"
)

plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG37b-bias_IE0.png"),
  power = power.bias,
  res.type = "bias",
  ce.type = "IE0"
)

plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG37c-bias_IE1.png"),
  power = power.sd,
  res.type = "bias",
  ce.type = "IE1"
)


### rmse

# estimates
power.rmse <- depthStatCalc(
  p05 = pow05, 
  p1 = pow1,
  p15 = pow15,
  p2 = pow2,
  p25 = pow25,
  p3 = pow3,
  p4 = pow4,
  p5 = pow5,
  X = X.sim,
  res.type = "rmse"
)

# plot results
plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG38a-rmse_DE.png"),
  power = power.rmse,
  res.type = "rmse",
  ce.type = "DE"
)

plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG38b-rmse_IE0.png"),
  power = power.rmse,
  res.type = "rmse",
  ce.type = "IE0"
)

plotDepthRes(
  plot.name = here::here("output", "supplement", "sFIG38c-rmse_IE1.png"),
  power = power.rmse,
  res.type = "rmse",
  ce.type = "IE1"
)

### Influence of (alpha, beta) on # terminal nodes

set.seed(35)

# simulate prior probabilities for number of terminal nodes
tree_sim1 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 0.5)
tree_sim2 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 1)
tree_sim3 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 1.5)
tree_sim4 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 2)
tree_sim5 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 2.5)
tree_sim6 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 3)
tree_sim7 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 4)
tree_sim8 <- simTermNodes(n_sims = 100000, alpha = 0.95, beta = 5)

# plot number of terminal nodes

plotTermNodes(
  plot.name = here::here("output", "supplement", "sFIG34-tuning_nodes.png"),
  s1 = tree_sim1,
  s2 = tree_sim2,
  s3 = tree_sim3,
  s4 = tree_sim4,
  s5 = tree_sim5,
  s6 = tree_sim6,
  s7 = tree_sim7,
  s8 = tree_sim8
)

### Influence of (alpha, beta) on tree depth

# simulate prior probability for max tree depth
depth_sim1 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 0.5)
depth_sim2 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 1)
depth_sim3 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 1.5)
depth_sim4 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 2)
depth_sim5 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 2.5)
depth_sim6 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 3)
depth_sim7 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 4)
depth_sim8 <- simMaxDepth(n_sims = 10000, alpha = 0.95, beta = 5)

# plot tree depth
plotDepthDist(
  plot.name = here::here("output", "supplement", "sFIG35-tuning_tree_depth.png"),
  s1 = depth_sim1,
  s2 = depth_sim2,
  s3 = depth_sim3,
  s4 = depth_sim4,
  s5 = depth_sim5,
  s6 = depth_sim6,
  s7 = depth_sim7,
  s8 = depth_sim8
)

