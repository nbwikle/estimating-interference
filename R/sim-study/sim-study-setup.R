### sim-study-setup.R
### Nathan Wikle

#######################################################################
### 1. Compile 'lm.stan' for repeated use in simulation study
#######################################################################

# load package
library(rstan)
# avoids recompiling unchanged stan code
rstan_options(auto_write = TRUE)
# compile stan LM model from source code; saved as 'lm.rds' in 'src'
stan_model(file = here::here("src", "lm.stan"), model_name = "lm", save_dso = TRUE)

#######################################################################
### 2. Create copula model approximating sulfate model's posterior 
#######################################################################


### Approximate sulfae model posterior with copula

# load package
library(copula)
library(truncnorm)

# gaussian copula structure with correlation matching actual SO4 param posterior
norm_cop <- normalCopula(
  c(0.25, 0, -0.15),
  dim = 3,
  dispstr = "un"
)

# marginal parameterization
cop_params <- list(
  gamma_mean = 2600,
  gamma_sd = 600,
  xi_mean = 175, 
  xi_var = 40000,
  alpha_mean = 25,
  alpha_var = 750
)
# log-normal mean and variance
cop_params$log_var = log(cop_params$xi_var / cop_params$xi_mean^2 + 1)
cop_params$log_mu = log(cop_params$xi_mean) - cop_params$log_var / 2
# gamma shape and rate parameters
cop_params$alpha_a = cop_params$alpha_mean^2 / cop_params$alpha_var
cop_params$alpha_b = cop_params$alpha_mean / cop_params$alpha_var

# specify marginal densities, create copula
copula_post <- mvdc(
  norm_cop,
  c("truncnorm", "lnorm", "gamma"), 
  list(
    list(a = 0, b = Inf, mean = cop_params$gamma_mean, sd = cop_params$gamma_sd),  
    list(meanlog = cop_params$log_mu, sdlog = sqrt(cop_params$log_var)), 
    list(shape = cop_params$alpha_a, rate = cop_params$alpha_b)
  )
)

# save copula
saveRDS(copula_post, file = here::here("output", "copula-posterior.RDS"))

