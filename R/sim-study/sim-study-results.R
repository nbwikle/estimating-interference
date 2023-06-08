### sim-study-results.R
### Nathan Wikle
###
### Processes the simulation study results.


#####################################################################
### 1. Load relevant data, set-up sim. study objects
#####################################################################

### create a new folder for sim study figures
dir.create(here::here("output", "supplement"), showWarnings = FALSE)
dir.create(here::here("output", "supplement", "sim-study"), showWarnings = FALSE)

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
### 2. CM1
#####################################################################

### i. Load data

# linear regression 
lm.cm1 <- readRDS(here::here("output", "simstudy-cm1-lm.RDS"))

# BART regression
bart.cm1 <- readRDS(here::here("output", "simstudy-cm1-bart.RDS"))

### ii. save bias, variance, and coverage results

# bias
bias <- list(
  CM1 = list(),
  CM2 = list(),
  CM3 = list()
)

# variance
variance <- list(
  CM1 = list(),
  CM2 = list(),
  CM3 = list()
)

# coverage
coverage <- list(
  CM1 = list(),
  CM2 = list(),
  CM3 = list()
)

### iii. bias/variance/coverage using linear model

bias$CM1$lm <- biasCalc(
  results = lm.cm1, 
  te = trueCM1(
    X = as.matrix(X.sim[,-1]), 
    params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = FALSE
)

coverage$CM1$lm <- coverageCalc(
  results = lm.cm1, 
  te = trueCM1(
    X = as.matrix(X.sim[,-1]), 
    params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = FALSE
)

variance$CM1$lm <- varianceCalc(
  results = lm.cm1,
  loglin = FALSE
)

### iv. bias/variance/coverage using BART model

bias$CM1$bart <- biasCalc(
  results = bart.cm1, 
  te = trueCM1(
    X = as.matrix(X.sim[,-1]), 
    params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = FALSE
)

coverage$CM1$bart <- coverageCalc(
  results = bart.cm1, 
  te = trueCM1(
    X = as.matrix(X.sim[,-1]), 
    params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = FALSE
)

variance$CM1$bart <- varianceCalc(
  results = bart.cm1,
  loglin = FALSE
)


### v. supplementary material plots

te.cm1 <- trueCM1(
  X = as.matrix(X.sim[,-1]),
  params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7)
)

k <- 100

# DE(g)
png(
  file = here::here("output", "supplement", "sFIG13a-lm_cm1_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm1[[k]]$plugin$DE,
  unc = lm.cm1[[k]]$cut$DE,
  true = te.cm1$DE,
  title = "Linear Regression:  CM1",
  y_title = "DE(g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG13b-bart_cm1_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm1[[k]]$plugin$DE,
  unc = bart.cm1[[k]]$cut$DE,
  true = te.cm1$DE,
  title = "BART:  CM1",
  y_title = "DE(g)",
  loglin = FALSE
)
dev.off()

# IE(0,g)
png(
  file = here::here("output", "supplement", "sFIG14a-lm_cm1_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm1[[k]]$plugin$IE0,
  unc = lm.cm1[[k]]$cut$IE0,
  true = te.cm1$IE0,
  title = "Linear Regression:  CM1",
  y_title = "IE(0,g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG14b-bart_cm1_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm1[[k]]$plugin$IE0,
  unc = bart.cm1[[k]]$cut$IE0,
  true = te.cm1$IE0,
  title = "BART:  CM1",
  y_title = "IE(0,g)",
  loglin = FALSE
)
dev.off()

# IE(1,g)
png(
  file = here::here("output", "supplement", "sFIG15a-lm_cm1_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm1[[k]]$plugin$IE1,
  unc = lm.cm1[[k]]$cut$IE1,
  true = te.cm1$IE1,
  title = "Linear Regression:  CM1",
  y_title = "IE(1,g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG15b-bart_cm1_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm1[[k]]$plugin$IE1,
  unc = bart.cm1[[k]]$cut$IE1,
  true = te.cm1$IE1,
  title = "BART:  CM1",
  y_title = "IE(1,g)",
  loglin = FALSE
)
dev.off()

### vi. remove data structures
rm(lm.cm1)
rm(bart.cm1)


#####################################################################
### 3. CM2
#####################################################################

### i. Load data

# linear regression 
lm.cm2 <- readRDS(here::here("output", "simstudy-cm2-lm.RDS"))

# BART regression
bart.cm2 <- readRDS(here::here("output", "simstudy-cm2-bart.RDS"))

### ii. bias/variance/coverage using linear model

bias$CM2$lm <- biasCalc(
  results = lm.cm2, 
  te = trueCM2(
    X = as.matrix(X.sim[,-1]), 
    params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = FALSE
)

coverage$CM2$lm <- coverageCalc(
  results = lm.cm2, 
  te = trueCM2(
    X = as.matrix(X.sim[,-1]), 
    params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = FALSE
)

variance$CM2$lm <- varianceCalc(
  results = lm.cm2,
  loglin = FALSE
)

### iii. bias/variance/coverage using BART model

bias$CM2$bart <- biasCalc(
  results = bart.cm2, 
  te = trueCM2(
    X = as.matrix(X.sim[,-1]), 
    params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = FALSE
)

coverage$CM2$bart <- coverageCalc(
  results = bart.cm2, 
  te = trueCM2(
    X = as.matrix(X.sim[,-1]), 
    params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = FALSE
)

variance$CM2$bart <- varianceCalc(
  results = bart.cm2,
  loglin = FALSE
)

### iv. supplementary material plots

te.cm2 <- trueCM2(
  X = as.matrix(X.sim[,-1]),
  params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6)
)

k <- 100

# DE(g)
png(
  file = here::here("output", "supplement", "sFIG16a-lm_cm2_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm2[[k]]$plugin$DE,
  unc = lm.cm2[[k]]$cut$DE,
  true = te.cm2$DE,
  title = "Linear Regression:  CM2",
  y_title = "DE(g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG16b-bart_cm2_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm2[[k]]$plugin$DE,
  unc = bart.cm2[[k]]$cut$DE,
  true = te.cm2$DE,
  title = "BART:  CM2",
  y_title = "DE(g)",
  loglin = FALSE
)
dev.off()

# IE(0,g)
png(
  file = here::here("output", "supplement", "sFIG17a-lm_cm2_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm2[[k]]$plugin$IE0,
  unc = lm.cm2[[k]]$cut$IE0,
  true = te.cm2$IE0,
  title = "Linear Regression:  CM2",
  y_title = "IE(0,g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG17b-bart_cm2_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm2[[k]]$plugin$IE0,
  unc = bart.cm2[[k]]$cut$IE0,
  true = te.cm2$IE0,
  title = "BART:  CM2",
  y_title = "IE(0,g)",
  loglin = FALSE
)
dev.off()

# IE(1,g)
png(
  file = here::here("output", "supplement", "sFIG18a-lm_cm2_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm2[[k]]$plugin$IE1,
  unc = lm.cm2[[k]]$cut$IE1,
  true = te.cm2$IE1,
  title = "Linear Regression:  CM2",
  y_title = "IE(1,g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG18b-bart_cm2_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm2[[k]]$plugin$IE1,
  unc = bart.cm2[[k]]$cut$IE1,
  true = te.cm2$IE1,
  title = "BART:  CM2",
  y_title = "IE(1,g)",
  loglin = FALSE
)
dev.off()

### v. remove data structures
rm(lm.cm2)
rm(bart.cm2)


#####################################################################
### 4. CM3
#####################################################################

### i. Load data

# linear regression 
lm.cm3 <- readRDS(here::here("output", "simstudy-cm3-lm.RDS"))

# BART regression
bart.cm3 <- readRDS(here::here("output", "simstudy-cm3-bart.RDS"))

### ii. bias/variance/coverage using linear model

bias$CM3$lm <- biasCalc(
  results = lm.cm3, 
  te = trueCM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = FALSE
)

coverage$CM3$lm <- coverageCalc(
  results = lm.cm3, 
  te = trueCM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = FALSE
)

variance$CM3$lm <- varianceCalc(
  results = lm.cm3,
  loglin = FALSE
)

### iii. bias/variance/coverage using BART model

bias$CM3$bart <- biasCalc(
  results = bart.cm3, 
  te = trueCM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = FALSE
)

coverage$CM3$bart <- coverageCalc(
  results = bart.cm3, 
  te = trueCM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = FALSE
)

variance$CM3$bart <- varianceCalc(
  results = bart.cm3,
  loglin = FALSE
)

### iv. supplementary material plots

te.cm3 <- trueCM3(
  X = as.matrix(X.sim[,-1]),
)

k <- 100

# DE(g)
png(
  file = here::here("output", "supplement", "sFIG19a-lm_cm3_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm3[[k]]$plugin$DE,
  unc = lm.cm3[[k]]$cut$DE,
  true = te.cm3$DE,
  title = "Linear Regression:  CM3",
  y_title = "DE(g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG19b-bart_cm3_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm3[[k]]$plugin$DE,
  unc = bart.cm3[[k]]$cut$DE,
  true = te.cm3$DE,
  title = "BART:  CM3",
  y_title = "DE(g)",
  loglin = FALSE
)
dev.off()

# IE(0,g)
png(
  file = here::here("output", "supplement", "sFIG20a-lm_cm3_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm3[[k]]$plugin$IE0,
  unc = lm.cm3[[k]]$cut$IE0,
  true = te.cm3$IE0,
  title = "Linear Regression:  CM3",
  y_title = "IE(0,g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG20b-bart_cm3_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm3[[k]]$plugin$IE0,
  unc = bart.cm3[[k]]$cut$IE0,
  true = te.cm3$IE0,
  title = "BART:  CM3",
  y_title = "IE(0,g)",
  loglin = FALSE
)
dev.off()

# IE(1,g)
png(
  file = here::here("output", "supplement", "sFIG21a-lm_cm3_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = lm.cm3[[k]]$plugin$IE1,
  unc = lm.cm3[[k]]$cut$IE1,
  true = te.cm3$IE1,
  title = "Linear Regression:  CM3",
  y_title = "IE(1,g)",
  loglin = FALSE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG21b-bart_cm3_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.cm3[[k]]$plugin$IE1,
  unc = bart.cm3[[k]]$cut$IE1,
  true = te.cm3$IE1,
  title = "BART:  CM3",
  y_title = "IE(1,g)",
  loglin = FALSE
)
dev.off()

### v. remove data structures
rm(lm.cm3)
rm(bart.cm3)


#####################################################################
### 5. PM1
#####################################################################

### i. Load data

# Poisson regression
pois.pm1 <- readRDS(here::here("output", "simstudy-pm1-pois.RDS"))

# log-linear BART regression
bart.pm1 <- readRDS(here::here("output", "simstudy-pm1-bart.RDS"))

### ii. bias/variance/coverage using linear model

bias$PM1$pois <- biasCalc(
  results = pois.pm1, 
  te = truePM1(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = TRUE
)

coverage$PM1$pois <- coverageCalc(
  results = pois.pm1, 
  te = truePM1(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = TRUE
)

variance$PM1$pois <- varianceCalc(
  results = pois.pm1,
  loglin = TRUE
)

### iii. bias/variance/coverage using BART model

bias$PM1$bart <- biasCalc(
  results = bart.pm1, 
  te = truePM1(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = TRUE
)

coverage$PM1$bart <- coverageCalc(
  results = bart.pm1, 
  te = truePM1(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = TRUE
)

variance$PM1$bart <- varianceCalc(
  results = bart.pm1,
  loglin = TRUE
)

### iv. supplementary material plots

# true effects for PM1
te.pm1 <- truePM1(X = as.matrix(X.sim[,-1]))
k <- 4

# DE(g)
png(
  file = here::here("output", "supplement", "sFIG22a-pois_pm1_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm1[[k]]$plugin$DE,
  unc = pois.pm1[[k]]$cut$DE,
  true = te.pm1$DE,
  title = "Poisson Regression:  PM1",
  y_title = "DE(g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG22b-bart_pm1_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm1[[k]]$plugin$DE,
  unc = bart.pm1[[k]]$cut$DE,
  true = te.pm1$DE,
  title = "Log-linear BART:  PM1",
  y_title = "DE(g)",
  loglin = TRUE
)
dev.off()

# IE(0,g)
png(
  file = here::here("output", "supplement", "sFIG23a-pois_pm1_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm1[[k]]$plugin$IE0,
  unc = pois.pm1[[k]]$cut$IE0,
  true = te.pm1$IE0,
  title = "Poisson Regression:  PM1",
  y_title = "IE(0,g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG23b-bart_pm1_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm1[[k]]$plugin$IE0,
  unc = bart.pm1[[k]]$cut$IE0,
  true = te.pm1$IE0,
  title = "Log-linear BART:  PM1",
  y_title = "IE(0,g)",
  loglin = TRUE
)
dev.off()

# IE(1,g)
png(
  file = here::here("output", "supplement", "sFIG24a-pois_pm1_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm1[[k]]$plugin$IE1,
  unc = pois.pm1[[k]]$cut$IE1,
  true = te.pm1$IE1,
  title = "Poisson Regression:  PM1",
  y_title = "IE(1,g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG24b-bart_pm1_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm1[[k]]$plugin$IE1,
  unc = bart.pm1[[k]]$cut$IE1,
  true = te.pm1$IE1,
  title = "Log-linear BART:  PM1",
  y_title = "IE(1,g)",
  loglin = TRUE
)
dev.off()

### v. remove data structures
rm(pois.pm1)
rm(bart.pm1)


#####################################################################
### 6. PM2
#####################################################################

### i. Load data

# Poisson regression
pois.pm2 <- readRDS(here::here("output", "simstudy-pm2-pois.RDS"))

# log-linear BART regression
bart.pm2 <- readRDS(here::here("output", "simstudy-pm2-bart.RDS"))

### ii. bias/variance/coverage using linear model

bias$PM2$pois <- biasCalc(
  results = pois.pm2, 
  te = truePM2(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = TRUE
)

coverage$PM2$pois <- coverageCalc(
  results = pois.pm2, 
  te = truePM2(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = TRUE
)

variance$PM2$pois <- varianceCalc(
  results = pois.pm2,
  loglin = TRUE
)

### iii. bias/variance/coverage using BART model

bias$PM2$bart <- biasCalc(
  results = bart.pm2, 
  te = truePM2(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = TRUE
)

coverage$PM2$bart <- coverageCalc(
  results = bart.pm2, 
  te = truePM2(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = TRUE
)

variance$PM2$bart <- varianceCalc(
  results = bart.pm2,
  loglin = TRUE
)

### iv. supplementary material plots

# true effects for PM1
te.pm2 <- truePM2(X = as.matrix(X.sim[,-1]))
k <- 4

# DE(g)
png(
  file = here::here("output", "supplement", "sFIG25a-bart_pm2_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm2[[k]]$plugin$DE,
  unc = pois.pm2[[k]]$cut$DE,
  true = te.pm2$DE,
  title = "Poisson Regression:  PM2",
  y_title = "DE(g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG25b-bart_pm2_DE.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm2[[k]]$plugin$DE,
  unc = bart.pm2[[k]]$cut$DE,
  true = te.pm2$DE,
  title = "Log-linear BART:  PM2",
  y_title = "DE(g)",
  loglin = TRUE
)
dev.off()

# IE(0,g)
png(
  file = here::here("output", "supplement", "sFIG26b-pois_pm2_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm2[[k]]$plugin$IE0,
  unc = pois.pm2[[k]]$cut$IE0,
  true = te.pm2$IE0,
  title = "Poisson Regression:  PM2",
  y_title = "IE(0,g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG26b-bart_pm2_IE0.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm2[[k]]$plugin$IE0,
  unc = bart.pm2[[k]]$cut$IE0,
  true = te.pm2$IE0,
  title = "Log-linear BART:  PM2",
  y_title = "IE(0,g)",
  loglin = TRUE
)
dev.off()

# IE(1,g)
png(
  file = here::here("output", "supplement", "sFIG27a-pois_pm2_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm2[[k]]$plugin$IE1,
  unc = pois.pm2[[k]]$cut$IE1,
  true = te.pm2$IE1,
  title = "Poisson Regression:  PM2",
  y_title = "IE(1,g)",
  loglin = TRUE
)
dev.off()


png(
  file = here::here("output", "supplement", "sFIG27b-bart_pm2_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm2[[k]]$plugin$IE1,
  unc = bart.pm2[[k]]$cut$IE1,
  true = te.pm2$IE1,
  title = "Log-linear BART:  PM2",
  y_title = "IE(1,g)",
  loglin = TRUE
)
dev.off()

### v. remove data structures
rm(pois.pm2)
rm(bart.pm2)


#####################################################################
### 6. PM3
#####################################################################

### i. Load data

# Poisson regression
pois.pm3 <- readRDS(here::here("output", "simstudy-pm3-pois.RDS"))

# log-linear BART regression
bart.pm3 <- readRDS(here::here("output", "simstudy-pm3-bart.RDS"))

### ii. bias/variance/coverage using linear model

bias$PM3$pois <- biasCalc(
  results = pois.pm3, 
  te = truePM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = TRUE
)

coverage$PM3$pois <- coverageCalc(
  results = pois.pm3, 
  te = truePM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = TRUE
)

variance$PM3$pois <- varianceCalc(
  results = pois.pm3,
  loglin = TRUE
)

### iii. bias/variance/coverage using BART model

bias$PM3$bart <- biasCalc(
  results = bart.pm3, 
  te = truePM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ), 
  loglin = TRUE
)

coverage$PM3$bart <- coverageCalc(
  results = bart.pm3, 
  te = truePM3(
    X = as.matrix(X.sim[,-1]), 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  ),
  loglin = TRUE
)

variance$PM3$bart <- varianceCalc(
  results = bart.pm3,
  loglin = TRUE
)

### iv. supplementary material plots

# true effects for PM3
te.pm3 <- truePM3(X = as.matrix(X.sim[,-1]))
k <- 4

# DE(g)
png(
  file = here::here("output", "supplement", "sFIG28a-pois_pm3_DE.png"), 
    width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm3[[k]]$plugin$DE,
  unc = pois.pm3[[k]]$cut$DE,
  true = te.pm3$DE,
  title = "Poisson Regression:  PM3",
  y_title = "DE(g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG28b-bart_pm3_DE.png"), 
    width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm3[[k]]$plugin$DE,
  unc = bart.pm3[[k]]$cut$DE,
  true = te.pm3$DE,
  title = "Log-linear BART:  PM3",
  y_title = "DE(g)",
  loglin = TRUE
)
dev.off()

# IE(0,g)
png(
  file = here::here("output", "supplement", "sFIG29a-pois_pm3_IE0.png"), 
    width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm3[[k]]$plugin$IE0,
  unc = pois.pm3[[k]]$cut$IE0,
  true = te.pm3$IE0,
  title = "Poisson Regression:  PM3",
  y_title = "IE(0,g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG29b-bart_pm3_IE0.png"), 
    width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm3[[k]]$plugin$IE0,
  unc = bart.pm3[[k]]$cut$IE0,
  true = te.pm3$IE0,
  title = "Log-linear BART:  PM3",
  y_title = "IE(0,g)",
  loglin = TRUE
)
dev.off()

# IE(1,g)
png(
  file = here::here("output", "supplement", "sFIG30a-pois_pm3_IE1.png"), 
    width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = pois.pm3[[k]]$plugin$IE1,
  unc = pois.pm3[[k]]$cut$IE1,
  true = te.pm3$IE1,
  title = "Poisson Regression:  PM3",
  y_title = "IE(1,g)",
  loglin = TRUE
)
dev.off()

png(
  file = here::here("output", "supplement", "sFIG30b-bart_pm3_IE1.png"), 
  width = 8, height = 5, units = "in", res = 300, bg = "white"
)
plotSimStudyEst(
  no_unc = bart.pm3[[k]]$plugin$IE1,
  unc = bart.pm3[[k]]$cut$IE1,
  true = te.pm3$IE1,
  title = "Log-linear BART:  PM3",
  y_title = "IE(1,g)",
  loglin = TRUE
)
dev.off()

### v. remove data structures
rm(pois.pm3)
rm(bart.pm3)


#####################################################################
### 7. Save Bias, Coverage, and Variance results
#####################################################################

# combine into a list
ss.results <- list(
  bias = bias,
  coverage = coverage,
  variance = variance
)

# save
saveRDS(ss.results, file = here::here("output", "simstudy-results.RDS"))



