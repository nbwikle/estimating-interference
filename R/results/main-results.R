### main-results.R
### Nathan Wikle


### Create a new folder for output in main text
dir.create(here::here("output", "main"), showWarnings = FALSE)

### Load data

# texas so4 data
tx.so4 <- readRDS(here::here("data", "texas-2016-data.RDS"))
# air pollution transport data
so4.mcmc <- readRDS(here::here("output", "tx-2016-100k-int-mx.RDS")) 
# outcome data
out.df <- readRDS(here::here("data", "outcome-data.RDS"))
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
# Texas shapefile 
tx.sf <- readRDS(here::here("data", "tx-state-sf.RDS"))


######################################################################
### 1. Scrubber Locations and Asthma Rates (Figure 1)
######################################################################

### A. Scrubber Locations (Figure 1a)

# facility locations
fac.locs <- fac.df[, 20:19]
names(fac.locs) <- c("long", "lat")

# convert to proper projection
fac.coords <- fac.locs %>%
  st_as_sf(
    coords = c("long", "lat"),
    crs = st_crs(out.df)
  )

scrubberPlot(
  plot.name = here::here("output", "main", "FIG1a-fac_locations.png"),
  coordinates = fac.coords,
  scrubbers = scrubber.vec
)

# ### B. Rate of Asthma ED visits in 2016 (Figure 1b)
# 
# # Note: Due to the private nature of the pediatric asthma data, it cannot 
# #   be recreated in this repository. However, I have included the relevant code. 
# 
# # create sf object with asthma rate per ZCTA
# ped.asthma <- out.df %>% 
#   mutate(asthma_rate = asthma / Ped_16) %>%
#   dplyr::select(c(
#     "Ped_16", "asthma", "asthma_rate"
#   ))
# # convert Inf values to NA (zero children in ZCTA)
# ped.asthma$asthma_rate[ped.asthma$asthma_rate == Inf] <- NA
# 
# # create figure
# asthmaPlot(
#   plot.name = here::here("output", "main", "FIG1b-asthma_rate.png"),
#   shp_file = ped.asthma, tx_shp = tx.sf, type = "rates", geo = "zcta"
# )



######################################################################
### 2. Sulfate Ananlysis (Figure 2)
######################################################################


### A. Average sulfate concentrations (Figure 2a)

sulfateDataPlot(
  plot.name = here::here("output", "main", "FIG2a-so4_wind_comb.png"),
  tx = tx.so4
)

### B. Estimated SO4 due to power plant emissions (Figure 2b)

expectedSO4Plot(
  plot.name = here::here("output", "main", "FIG2b-mean_so4.png"),
  mats = adv.mats,
  so4.data = tx.so4,
  so4.res = so4.mcmc
)


######################################################################
### 3. Key-associated and Upwind Treatments (Figure 3)
######################################################################

# remove NA rows from out.df

out.df <- readRDS(here::here("data", "outcome-data.RDS"))
out.df <- out.df %>% drop_units() 
out.df <- na.omit(out.df) 

# estimate G.bar and sigma(G)
burnin <- 25000
theta.mu <- colMeans(so4.mcmc$samples[-c(1:burnin), ])

G_mu <- gCalc(
  theta = theta.mu, 
  mats = adv.mats, 
  X_mat = em.mat, 
  P_mat = proj.mat, 
  s_vec = scrubber.vec, 
  key = key.assoc, 
  version = 3
)

out.df$G <- G_mu



### A. Key-associated treatments (Figure 3a)

plotTreatment(
  plot_name = here::here("output", "main", "FIG3a-z_plot.png"),
  shp_file = out.df %>% mutate(Z = scrubbed),
  tx_shp = tx.sf,
  trt = "Z",
  title = ""
)

### B. Mean upwind treatments (Figure 3b)

plotTreatment(
  plot_name = here::here("output", "main", "FIG3b-g_plot.png"),
  shp_file = out.df,
  tx_shp = tx.sf,
  trt = "G",
  title = ""
)

### C. Upwind treatment standard deviation (Figure 3c)

# sample theta values
set.seed(23)
theta_samples <- so4.mcmc$samples[
  25000 + sample(nrow(so4.mcmc$samples) - 25000, 1000),
]

# estimate g_i for each theta
g_vals <- apply(
  theta_samples, 1, gCalc,
  mats = adv.mats,
  X_mat = em.mat, 
  P_mat = proj.mat,
  s_vec = scrubber.vec, 
  key = key.assoc,
  version = 3
)
# marginal standard deviation of g_i
g_sd <- apply(g_vals, 1, sd)
out.df$G_sd <- g_sd

# generate plot
plotTreatment(
  plot_name = here::here("output", "main", "FIG3c-g_sd_plot.png"),
  shp_file = out.df,
  tx_shp = tx.sf,
  trt = "G sd",
  title = ""
)




######################################################################
### 4. Medicare Results (Figure 4)
######################################################################

### Medicare Data

# log-linear BART
bart.plugin <- readRDS(here::here("output", "medicare-bart-plugin.RDS"))
bart.cut <- readRDS(here::here("output", "medicare-bart-cut.RDS"))

# Poisson regression
pois.plugin <- readRDS(here::here("output", "medicare-pois-plugin.RDS"))
pois.cut <- readRDS(here::here("output", "medicare-pois-cut.RDS"))


### Estimate proportion of variance attributed to uncertainty propagation

# BART results
mult <- 1000
bart.var <- list()
bart.var$DE <- varComponents(res = mult * bart.cut$DE$all, n_samps = 260)
bart.var$IE0 <- varComponents(res = mult * bart.cut$IE0$all, n_samps = 260)
bart.var$IE1 <- varComponents(res = mult * bart.cut$IE1$all, n_samps = 260)

# Poisson regression results
mult <- 1000
pois.var <- list()
pois.var$DE <- varComponents(res = mult * pois.cut$DE$all, n_samps = 250)
pois.var$IE0 <- varComponents(res = mult * pois.cut$IE0$all, n_samps = 250)
pois.var$IE1 <- varComponents(res = mult * pois.cut$IE1$all, n_samps = 250)


# compare % variance due to interference across g, estimation strategy 
prop.var.DE <- cbind(
  seq(0.25, 0.9, by = 0.01),
  pois.var$DE$B / pois.var$DE$T,
  bart.var$DE$B / bart.var$DE$T
)
colnames(prop.var.DE) <- c("g", "Poisson", "BART")

prop.var.IE0 <- cbind(
  seq(0.25, 0.9, by = 0.01),
  pois.var$IE0$B / pois.var$IE0$T,
  bart.var$IE0$B / bart.var$IE0$T
)
colnames(prop.var.IE0) <- c("g", "Poisson", "BART")

prop.var.IE1 <- cbind(
  seq(0.25, 0.9, by = 0.01),
  pois.var$IE1$B / pois.var$IE1$T,
  bart.var$IE1$B / bart.var$IE1$T
)
colnames(prop.var.IE1) <- c("g", "Poisson", "BART")


### Plot estimated effects

# colors
cols <- c("#D3B1C2","#613659", "#E1C391", "#705446")

# y-axis limits


y.lims.de <- c(
  min(
    min(1000 * apply(bart.plugin$DE$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(bart.cut$DE$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(pois.plugin$DE$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(pois.cut$DE$all, 2, quantile, probs = 0.025))
  ),
  max(
    max(1000 * apply(bart.plugin$DE$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(bart.cut$DE$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(pois.plugin$DE$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(pois.cut$DE$all, 2, quantile, probs = 0.975))
  )
)

y.lims.ie0 <- c(
  min(
    min(1000 * apply(bart.plugin$IE0$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(bart.cut$IE0$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(pois.plugin$IE0$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(pois.cut$IE0$all, 2, quantile, probs = 0.025))
  ),
  max(
    max(1000 * apply(bart.plugin$IE0$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(bart.cut$IE0$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(pois.plugin$IE0$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(pois.cut$IE0$all, 2, quantile, probs = 0.975))
  )
)


y.lims.ie1 <- c(
  min(
    min(1000 * apply(bart.plugin$IE1$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(bart.cut$IE1$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(pois.plugin$IE1$all, 2, quantile, probs = 0.025)),
    min(1000 * apply(pois.cut$IE1$all, 2, quantile, probs = 0.025))
  ),
  max(
    max(1000 * apply(bart.plugin$IE1$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(bart.cut$IE1$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(pois.plugin$IE1$all, 2, quantile, probs = 0.975)),
    max(1000 * apply(pois.cut$IE1$all, 2, quantile, probs = 0.975))
  )
)



# match space between text across plots
de.ratio <- 1.15 / 14.7622 * (y.lims.de[2] - y.lims.de[1])
ie0.ratio <- de.ratio / (y.lims.de[2] - y.lims.de[1]) * (y.lims.ie0[2] - y.lims.ie0[1])
ie1.ratio <- de.ratio / (y.lims.de[2] - y.lims.de[1]) * (y.lims.ie1[2] - y.lims.ie1[1])


plotMedicareEffects(
  plot.name = here::here("output", "main", "FIG4-medicare_results.png"), 
  pois.plugin = pois.plugin, 
  pois.cut = pois.cut,
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  prop.var.DE = prop.var.DE,
  prop.var.IE0 = prop.var.IE0,
  prop.var.IE1 = prop.var.IE1,
  y.lims.de = y.lims.de, 
  y.lims.ie0 = y.lims.ie0,
  y.lims.ie1 = y.lims.ie1,
  cols = cols
)


######################################################################
### 5. Asthma Results (Figure 5)
######################################################################

### Asthma Data

# log-linear BART
bart.plugin <- readRDS(here::here("output", "asthma-bart-plugin.RDS"))
bart.cut <- readRDS(here::here("output", "asthma-bart-cut.RDS"))

# Poisson regression
pois.plugin <- readRDS(here::here("output", "asthma-pois-plugin.RDS"))
pois.cut <- readRDS(here::here("output", "asthma-pois-cut.RDS"))

### Estimate proportion of variance attributed to uncertainty propagation

# BART results
mult <- 1000
bart.var <- list()
bart.var$DE <- varComponents(res = mult * bart.cut$DE$all, n_samps = 250)
bart.var$IE0 <- varComponents(res = mult * bart.cut$IE0$all, n_samps = 250)
bart.var$IE1 <- varComponents(res = mult * bart.cut$IE1$all, n_samps = 250)

# Poisson regression results
mult <- 1000
pois.var <- list()
pois.var$DE <- varComponents(res = mult * pois.cut$DE$all, n_samps = 250)
pois.var$IE0 <- varComponents(res = mult * pois.cut$IE0$all, n_samps = 250)
pois.var$IE1 <- varComponents(res = mult * pois.cut$IE1$all, n_samps = 250)

# compare % variance due to interference across g, estimation strategy 
prop.var.DE <- cbind(
  seq(0.25, 0.9, by = 0.01),
  pois.var$DE$B / pois.var$DE$T,
  bart.var$DE$B / bart.var$DE$T
)
colnames(prop.var.DE) <- c("g", "Poisson", "BART")

prop.var.IE0 <- cbind(
  seq(0.25, 0.9, by = 0.01),
  pois.var$IE0$B / pois.var$IE0$T,
  bart.var$IE0$B / bart.var$IE0$T
)
colnames(prop.var.IE0) <- c("g", "Poisson", "BART")

prop.var.IE1 <- cbind(
  seq(0.25, 0.9, by = 0.01),
  pois.var$IE1$B / pois.var$IE1$T,
  bart.var$IE1$B / bart.var$IE1$T
)
colnames(prop.var.IE1) <- c("g", "Poisson", "BART")


### Plot estimated causal effects (DE(g), IE(z,g))

# colors
cols <- c("#D3B1C2", "#613659", "#E1C391", "#705446")

y.lims.de <- c(-10, 4.7622)
y.lims.ie0 <- c(-5.6852, 15)
y.lims.ie1 <- c(-6.2551, 15)

de.ratio <- 1.15
ie0.ratio <- de.ratio / (y.lims.de[2] - y.lims.de[1]) * (y.lims.ie0[2] - y.lims.ie0[1])
ie1.ratio <- de.ratio / (y.lims.de[2] - y.lims.de[1]) * (y.lims.ie1[2] - y.lims.ie1[1])


plotAsthmaEffects(
  plot.name = here::here("output", "main", "FIG5-asthma_results.png"), 
  pois.plugin = pois.plugin, 
  pois.cut = pois.cut,
  bart.plugin = bart.plugin, 
  bart.cut = bart.cut,
  prop.var.DE = prop.var.DE,
  prop.var.IE0 = prop.var.IE0,
  prop.var.IE1 = prop.var.IE1,
  y.lims.de = y.lims.de, 
  y.lims.ie0 = y.lims.ie0,
  y.lims.ie1 = y.lims.ie1,
  cols = cols
)






