### so4-model.R
### Nathan Wikle

########################################################################
### 1. Set up MCMC sampler 
########################################################################

# load data
so4.ls <- readRDS(here::here("data", "texas-2016-data.RDS"))

# create D matrix
dims <- dim(so4.ls$so4)
inds.d <- FVMInds(dims, "insulated")
D <- diffFVM(dims, inds.d, "insulated")

# create C matrix
inds.w <- FVMInds(dims, "periodic")
C <- createC(
  big = so4.ls$wind.big$wind.2016,
  small = so4.ls$wind.small$wind.2016,
  dims, inds = inds.w
)

# X and y values
y.2016 <- values(so4.ls$so4)
X.2016 <- as.vector(so4.ls$X$X.usa + so4.ls$X$X.mexico)

# initialize parameters
gamma.init <- 1500
xi.init <- 0.45
alpha.init <- 0.45
beta.init <- 4
s2.init <- 25000
delta.init <- 50
b0.init <- 0.5 

params.init <- list(
  gamma = gamma.init,
  xi = xi.init,
  beta = beta.init,
  s2 = s2.init,
  delta = delta.init,
  alpha = alpha.init,
  b0 = b0.init
)

# generate MCMC samples !!! CAUTION: THIS TAKES A LONG TIME (~110 hrs) !!!
mcmc.res <- coupledMCMC.Wind(
  N = 100000, y = y.2016, X = X.2016,
  mats = list(D = D, C = C),
  params = params.init,
  priors = list(g.sd = 10000, xi = 0.01, s2 = 0.001, beta = 10, alpha.sd = 1000, b0 = 5),
  update = 50, iter = 100,
  n.y = 50, burnin = 10000, intercept = TRUE
)

# save results
saveRDS(mcmc.res, here::here("output", "tx-2016-100k-int-mx.RDS"))
