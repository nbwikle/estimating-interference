### sim-study-functions.R
### Nathan Wikle
###
### A collection of functions used in the simulation study.


#####################################################################
### 1. General functions; possibly used for all simulations
#####################################################################

gCalc <- function(theta, mats, X_mat, P_mat, s_vec, key, version = 3){
  # Calculate G_i for a given sample of theta.
  # Input:
  #   theta:   Vector of SO4 model parameters.
  #   mats:    List containing advection and diffusion FVM matrices: C and D.
  #   X_mat:   Matrix of counterfactual emissions scenarios. 
  #   P_mat:   Projection matrix mapping raster elements to ZCTAs.
  #   s_vec:   Vector listing treatment status of key-associated power plants.
  #   key:     Vector of key-associated power plants.
  #   version: The version of G_i to calculate (default = 3; version used in paper).
  # Output:
  #   A vector of G_i values, one for each ZCTA in the analysis.
  
  # SO4 model parameters
  gamma <- theta[1]
  xi <- theta[2]
  beta <- theta[3]
  s2 <- theta[4]
  delta <- theta[5]
  alpha <- theta[6]
  # b0 <- theta[7]
    
  # matrix structures used to solve advection-diffusion equation
  D <- mats$D
  C <- mats$C
    
  # solve for mu | theta
  n_cells <- nrow(D)
  A.1 <- gamma * D + alpha * C + Matrix::Diagonal(n_cells, delta)
  A.2 <- (gamma/xi) * D + (alpha/xi) * C + Matrix::Diagonal(n_cells, 1)
  sources <- Matrix::solve(A.2, X_mat * beta)
  mu <- Matrix::solve(A.1, sources)
  
  # convert to source-receptor matrix
  T_mat <- Matrix::crossprod(P_mat, mu)
    
  ### calculate G_i
  
  # weighted sum of treated facilities
  G_star <- T_mat %*% s_vec - T_mat[cbind(1:nrow(T_mat), key)] * s_vec[key]
  # weighted degree
  sr_total <- Matrix::rowSums(T_mat) - T_mat[cbind(1:nrow(T_mat), key)]

  if (version == 1) {
    G <- (G_star / max(G_star))[, 1]
  } else if (version == 2) {
    G <- as.vector(G_star)
  } else if (version == 3) {
    # version used in the paper!
    G <- (G_star / sr_total)[, 1]
  } else if (version == 4) {
    G <- sr_total
  }
    
  return(G)
}

ieEst <- function(mu){
  # Estimate indirect effects (IE(z,g)) using a given mu(z,g) estimate.
  # Input:
  #   mu: Estimated mu(z,g) values.
  # Output:
  #   IE(z,g) estimates.
  t(apply(mu, 1, FUN = function(vec){vec - vec[1]})) 
}

tuningMat <- function() {
  # Creates an initial 9x9 covariance matrix for the M-H proposal distribution.
  matrix(
    c(
      3.408334e-03, 3.675509e-04, -2.848031e-04, -1.695028e-05, -1.050344e-03, -1.206524e-07, -1.745315e-06, -9.912764e-05, 7.910869e-04,
      3.675509e-04, 2.953484e-03, 2.028150e-04, -3.585426e-05, -6.657295e-04, -1.574000e-07, -3.301334e-04, 1.566296e-04, 5.068728e-04,
      -2.848031e-04, 2.028150e-04, 1.280350e-04, 1.036525e-06, 9.895577e-05, -4.625312e-08, -1.549233e-04, -3.314847e-05, -1.004646e-05,
      -1.695028e-05, -3.585426e-05, 1.036525e-06, 1.333538e-05, 3.237167e-05, -5.227321e-09, -1.005353e-05, -2.812927e-06, -2.846182e-05,
      -1.050344e-03, -6.657295e-04, 9.895577e-05, 3.237167e-05, 2.495804e-03, -4.027825e-08, -8.372018e-05, -1.016680e-04, -1.122048e-03,
      -1.206524e-07, -1.574000e-07, -4.625312e-08, -5.227321e-09, -4.027825e-08, 3.182481e-10, 5.553762e-08, 5.310186e-08, -3.089040e-08,
      -1.745315e-06, -3.301334e-04, -1.549233e-04, -1.005353e-05, -8.372018e-05, 5.553762e-08, 2.469013e-04, 5.071728e-05, -5.130964e-05,
      -9.912764e-05, 1.566296e-04, -3.314847e-05, -2.812927e-06, -1.016680e-04, 5.310186e-08, 5.071728e-05, 2.116932e-04, 7.434809e-05,
      7.910869e-04, 5.068728e-04, -1.004646e-05, -2.846182e-05, -1.122048e-03, -3.089040e-08, -5.130964e-05, 7.434809e-05, 1.924868e-03
    ),
    nrow = 9, ncol = 9, byrow = TRUE
  )
}


#####################################################################
### 2. Continuous Outcome Models
#####################################################################

simCM1 <- function(
    X, Z, G, 
    params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7), 
    sd = 1
){
  # Simulates continuous outcomes (i.e., 'CM1') from a given set of covariates 
  #   and parameter values.
  # Input:
  #   X: Matrix of covariates.
  #   Z: Vector of Z (direct treatment) values.
  #   G: Vector of G (indirect treatment) values.
  #   params: Vector of "true" parameter values 
  #     (default = values used in the paper for CM1).
  #   sd: Measurement effor (default = 1).
  # Output:
  #   A vector of simulated outcomes.

  # Create design matrix
  X_new <- cbind(rep(1, nrow(X)), scale(X), Z, G, Z * G)
  colnames(X_new)[ncol(X_new)] <- "Z.G"
  # Mean function
  mu <- X_new %*% params
  # Sample outcomes
  y_vals <- rnorm(n = length(mu), mean = mu, sd = sd)
  return(y_vals)
}

simCM2 <- function(
    X, Z, G, 
    params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
    sd = 1
){
  # Simulates continuous outcomes (i.e., 'CM2') from a given set of covariates 
  #   and parameter values.
  # Input:
  #   X: Matrix of covariates.
  #   Z: Vector of Z (direct treatment) values.
  #   G: Vector of G (indirect treatment) values.
  #   params: Vector of "true" parameter values 
  #     (default = values used in the paper for CM2).
  #   sd: Measurement effor (default = 1).
  # Output:
  #   A vector of simulated outcomes.
  
  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  W <- matrix(0.5, nrow = nrow(X_new), ncol = ncol(X_new))
  center_vals <- attr(X_new, which = 'scaled:center')
  scale_vals <- attr(X_new, which = 'scaled:scale')
  
  X_beta <- X_new %*% params
  XW_beta <- (X_new - W * G) %*% params
  
  # means for Y(0) vs Y(1)
  mu0 <- exp(XW_beta)
  mu1 <- X_beta - G^2
  
  # mean function
  mu <- mu0 * (1 - Z) + mu1 * Z
  
  # sample from N(mu, 1) 
  y_vals <- rnorm(n = length(mu), mean = mu, sd = sd)
  
  return(y_vals)
}

simCM3 <- function(
    X, Z, G, 
    params = c(0.2, -0.1, 1.3, -1.1, 2.3, 1.3, 0.9), 
    params_v2 = c(-1.5, -0.7, -0.3, 0.4, -0.5, 2.0, 1.0),
    sd = 1
){
  # Simulates continuous outcomes (i.e., 'CM3') from a given set of covariates 
  #   and parameter values.
  # Input:
  #   X: Matrix of covariates.
  #   Z: Vector of Z (direct treatment) values.
  #   G: Vector of G (indirect treatment) values.
  #   params: Vector of "true" parameter values 
  #     (default = values used in the paper for CM3).
  #   params_v2: Vector of "true" parameter values when G < 0.65
  #     (default = values used in the paper for CM3).
  #   sd: Measurement effor (default = 1).
  # Output:
  #   A vector of simulated outcomes.
  
  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  center_vals <- attr(X_new, which = 'scaled:center')
  scale_vals <- attr(X_new, which = 'scaled:scale')  
  
  # stratification
  g1 <- (G < 0.65)
  g3 <- (G >= 0.75)
  g2 <- !g1 & !g3
  
  # X * beta
  X_beta <- X_new %*% params
  
  # mu1
  mu1 <- X_beta + (Z * X_new) %*% params_v2 - 0.3 * Z
  mu1[G < 0.4] <- mu1[G < 0.4] + 1.7
  
  # mu2
  Xb_mod <- X_beta
  Xb_mod[X_beta < -3] <- -3  
  mu2 <- exp(-Xb_mod) - (0.25 * Z) - (1 + Z) * sin(40 * (G - 0.7))
  
  # mu3
  mu3 <- X_beta - (0.5 * Z) - 0.5 * exp(2 * G)
  
  # combine into a single mean function
  mu <- rep(0, nrow(X))
  mu[g1] <- mu1[g1]
  mu[g2] <- mu2[g2]
  mu[g3] <- mu3[g3]

  # sample from N(mu, 1) 
  y_vals <- rnorm(n = length(mu), mean = mu, sd = sd)
  
  return(y_vals)
}

CM_SimStudy_BART <- function(
    ss_num, n_gsamples, data_list, theta_posterior, model = 1, 
    gseq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
) {
  # Performs a simulation study with (continuous) BART regression. This includes
  #   simulation of the outcomes, and inference both with and without 
  #   incorporation of uncertainty in G.
  # Input:
  #   ss_num: A scalar value for setting the seed of the random number generator.
  #   theta_0: Posterior sample used to generate g
  #   n_gsamples: The number of independent MCMC samplers to run when 
  #       incorporating uncertainty in G.
  #   data_list: A list with necessary data for the simulation study 
  #       (so4_mats, em_mat, proj_mat, s_vec, key_assoc, x, z, g_mean).
  #   theta_posterior: A copula model, representing the theta posterior from 
  #       an analysis of the SO4 data.
  #   model: Specifies the continuous outcome model to simulate from 
  #     (i.e., 1 = CM1, 2 = CM2, 3 = CM3; default = 1).
  #   gseq: G values at which to estimate (marginal) causal effects.
  # Output: 
  #   A list with i) the 'true' theta vector, ii) MCMC samples using E(G | *), 
  #     iii) MCMC samples incorporating uncertainty in G, and the posterior 
  #     coverage probability of delta using the posterior iv) without uncertaity
  #     in G and v) with uncertainty in G.
  
  # set seed
  set.seed(ss_num)
  
  ## i. simulate outcome
  
  # sample theta
  s0 <- copula::rMvdc(1, theta_posterior)
  theta_0 <- c(s0[1:2], 2.6, 49600, 50, s0[3], 0.65)
  
  # generate G
  g_0 <- gCalc(
    theta = theta_0,
    mats = data_list$so4_mats,
    X_mat = data_list$em_mat,
    P_mat = data_list$proj_mat,
    s_vec = data_list$s_vec,
    key = data_list$key_assoc, 
    version = 3
  )
  
  # simulate outcome
  if (model == 1){
    y <- simCM1(
      X = data_list$x[,-1],
      Z = data_list$z,
      G = g_0,
      params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7), 
      sd = 1
    )
  } else if (model == 2){
    y <- simCM2(
      X = data_list$x[,-1],
      Z = data_list$z,
      G = g_0,
      params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
      sd = 1
    )
  } else if (model == 3){
    y <- simCM3(
      X = data_list$x[,-1],
      Z = data_list$z,
      G = g_0, 
      params = c(0.2, -0.1, 1.3, -1.1, 2.3, 1.3, 0.9), 
      params_v2 = c(-1.5, -0.7, -0.3, 0.4, -0.5, 2.0, 1.0),
      sd = 1
    )
  }
  
  ## ii. estimate using two-step approach

  # create data frame with x, z, and g_bar
  X.mu <- data.frame(data_list$x[,-1], data_list$z, data_list$g_mean)
  colnames(X.mu)[c(ncol(X.mu) - 1, ncol(X.mu))] <- c("Z", "G")
  
  # fit BART model
  bart.plugin <- gbart(
    x.train = X.mu, y.train = y, 
    ndpost = 1000, keepevery = 1, nskip = 500, ntree = 200
  )
  
  # estimate DE, IE0, IE1
  eff.plugin <- bartContEffects(
    fit = bart.plugin, 
    X = data_list$x[,-1], 
    n_post = nrow(bart.plugin$yhat.train),
    g_seq = gseq
  )

  ## iii. estimate using cutting feedback approach

  # save results in a list
  eff.k <- list()

  for (k in 1:n_gsamples) {
    
    cat(paste("\n Iteration ", k, "\n", sep = ""))

    # sample theta from posterior
    post_sample <- copula::rMvdc(1, theta_posterior)
    theta_k <- c(post_sample[1:2], 2.6, 49600, 50, post_sample[3], 0.65)

    cat("   theta: ")
    cat(theta_k)

    # create G
    g_k <- gCalc(
      theta = theta_k,
      mats = data_list$so4_mats,
      X_mat = data_list$em_mat,
      P_mat = data_list$proj_mat,
      s_vec = data_list$s_vec,
      key = data_list$key_assoc,
      version = 3
    )

    # create data frame with x, z, and g
    X.df <- data.frame(data_list$x[,-1], data_list$z, g_k)
    colnames(X.df)[c(ncol(X.df) - 1, ncol(X.df))] <- c("Z", "G")

    # fit BART model
    bart.k <- gbart(
      x.train = X.df, y.train = y, 
      ndpost = 1000, keepevery = 1, nskip = 500, ntree = 200
    )
    
    # estimate DE(g), IE(z,g)
    eff.k[[k]] <- bartContEffects(
      fit = bart.k, 
      X = data_list$x[,-1], 
      n_post = nrow(bart.k$yhat.train),
      g_seq = gseq
    )
  }

  ## process results

  # combine into single matrices

  eff.cut <- list()

  eff.cut$mu0 <- do.call(rbind, lapply(eff.k, function(res){res$mu0}))
  eff.cut$mu1 <- do.call(rbind, lapply(eff.k, function(res){res$mu1}))
  eff.cut$DE <- do.call(rbind, lapply(eff.k, function(res){res$DE}))
  eff.cut$IE0 <- do.call(rbind, lapply(eff.k, function(res){res$IE0}))
  eff.cut$IE1 <- do.call(rbind, lapply(eff.k, function(res){res$IE1}))

  # compile all results in a list
  results <- list(
    theta_sim = theta_0,
    plugin = eff.plugin, 
    cut = eff.cut
  )
  
  # return results
  return(results)
}

bartContEffects <- function(
    fit, X, n_post, 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
){
  # Estimate the effects (DE(g), IE(z,g) from a fitted (continuous) BART model.
  # Input:
  #   fit: Fitted BART regression (output from 'gbart').
  #   X: Covariate matrix.
  #   n_post: Number of posterior samples.
  #   g_seq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1.
  
  # number of g values used in prediction
  n_g <- length(g_seq)
  
  # estimate mu0
  z <- 0
  mu0 <- matrix(NA_real_, nrow = n_post, ncol = n_g)
  for (k in 1:n_g){
    g_val <- g_seq[k]
    X_m <- cbind(X, z, g_val)
    colnames(X_m)[(ncol(X_m) - 1):ncol(X_m)] <- c("Z", "G")
    p0 <- predict(fit, X_m)
    mu0[, k] <- rowMeans(p0)
  }
  
  # estimate mu1
  mu1 <- matrix(NA_real_, nrow = n_post, ncol = n_g)
  z <- 1
  for (k in 1:n_g){
    g_val <- g_seq[k]
    X_m <- cbind(X, z, g_val)
    colnames(X_m)[(ncol(X_m) - 1):ncol(X_m)] <- c("Z", "G")
    p1 <- predict(fit, X_m)
    mu1[, k] <- rowMeans(p1)
  }
  
  # DE(g) = mu(1,g) - mu(0,g)
  DE <- mu1 - mu0
  # IE(0,g)
  IE0 <- ieEst(mu0)
  # IE(1,g)
  IE1 <- ieEst(mu1)
  
  # return results
  results <- list(
    mu0 = mu0,
    mu1 = mu1,
    DE = DE,
    IE0 = IE0,
    IE1 = IE1
  )
  
  return(results)
}

CM_SimStudy_LM <- function(
    ss_num, n_gsamples, data_list, theta_posterior, stan_loc, interact = FALSE, model = 1, 
    gseq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
) {
  # Performs a simulation study with (Bayesian) linear regression. This includes
  #   simulation of the outcomes, and inference both with and without 
  #   incorporation of uncertainty in G.
  # Output:
  #   ss_num: A scalar value for setting the seed of the random number generator.
  #   n_gsamples: The number of independent MCMC samplers to run when 
  #       incorporating uncertainty in G.
  #   data_list: A list with necessary data for the simulation study 
  #       (so4_mats, em_mat, proj_mat, s_vec, key_assoc, x, z, g_mean).
  #   theta_posterior: A copula model, representing the theta posterior from 
  #       an analysis of the SO4 data.
  #   stan_loc: Stan 'lm.rds' object.
  #   interact: Boolean indicating if the correct Z*G interaction should be 
  #       included in the linear regression outcome model.
  #   model: Specifies the continuous outcome model to simulate from 
  #     (i.e., 1 = CM1, 2 = CM2, 3 = CM3; default = 1).
  #   gseq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list with i) the 'true' theta vector, ii) MCMC samples using E(G | *), 
  #     iii) MCMC samples incorporating uncertainty in G, and the posterior 
  #     coverage probability of delta using the posterior iv) without uncertaity
  #     in G and v) with uncertainty in G.
  
  # set seed
  set.seed(ss_num)
  
  # load compiled stan DSO
  lm.mod <- readRDS(stan_loc)
  
  ## i. simulate outcome
  
  # sample theta
  s0 <- copula::rMvdc(1, theta_posterior)
  theta_0 <- c(s0[1:2], 2.6, 49600, 50, s0[3], 0.65)
  
  # generate G
  g_0 <- gCalc(
    theta = theta_0,
    mats = data_list$so4_mats,
    X_mat = data_list$em_mat,
    P_mat = data_list$proj_mat,
    s_vec = data_list$s_vec,
    key = data_list$key_assoc, 
    version = 3
  )
  
  # simulate outcome
  if (model == 1){
    y <- simCM1(
      X = data_list$x[,-1],
      Z = data_list$z,
      G = g_0,
      params = c(2.5,  1.6, -2.3,  0.8, -1.1, -1.5, -0.5, -1.2, -0.4, -1.7), 
      sd = 1
    )
  } else if (model == 2){
    y <- simCM2(
      X = data_list$x[,-1],
      Z = data_list$z,
      G = g_0,
      params = c(-0.3, 0.3, -0.3, -0.3, 0.2, -0.4, 0.6),
      sd = 1
    )
  } else if (model == 3){
    y <- simCM3(
      X = data_list$x[,-1],
      Z = data_list$z,
      G = g_0, 
      params = c(0.2, -0.1, 1.3, -1.1, 2.3, 1.3, 0.9), 
      params_v2 = c(-1.5, -0.7, -0.3, 0.4, -0.5, 2.0, 1.0),
      sd = 1
    )
  }
  
  ## ii. estimate using two-step approach
  
  if (interact){
    # set up data frame for stan
    X.mu <- data.frame(
      1, data_list$x[,-1], data_list$z, data_list$g_mean, 
      data_list$z * data_list$g_mean
    )
    colnames(X.mu)[(ncol(X.mu) - 2):ncol(X.mu)] <- c("Z", "G", "Z.G")
  } else {
    ## use same data frame as BART
    
    # set up data frame for stan
    X.mu <- data.frame(1, data_list$x[,-1], data_list$z, data_list$g_mean)
    colnames(X.mu)[(ncol(X.mu) - 1):ncol(X.mu)] <- c("Z", "G")
  }
  
  # sample from lm.mod
  lm.2s <- rstan::sampling(
    lm.mod, 
    data = list(
      N = length(y),
      Y = y,
      K = ncol(X.mu),
      X = X.mu,
      prior_only = FALSE
    ),
    chains = 1, 
    iter = 5000
  )
  
  # estimate DE(g), IE(z,g)
  eff.plugin <- stanLMEffects(
    fit = extract(lm.2s), 
    X = data_list$x[,-1], 
    n.post = 1000,
    model = model, 
    g_seq = gseq,
    interact = interact
  )
  
  ## iii. estimate using cutting feedback approach
  
  # save results in a list
  eff.k <- list()
  
  for (k in 1:n_gsamples) {
    
    cat(paste("\n Iteration ", k, "\n", sep = ""))
    
    # sample theta from posterior
    post_sample <- copula::rMvdc(1, theta_posterior)
    theta_k <- c(post_sample[1:2], 2.6, 49600, 50, post_sample[3], 0.65)
    
    cat("   theta: ")
    cat(theta_k)
    
    # create G
    g_k <- gCalc(
      theta = theta_k,
      mats = data_list$so4_mats,
      X_mat = data_list$em_mat,
      P_mat = data_list$proj_mat,
      s_vec = data_list$s_vec,
      key = data_list$key_assoc,
      version = 3
    )
    
    if (interact){
      ## use correct design matrix
      
      # set up data frame for stan
      X.df <- data.frame(1, data_list$x[,-1], data_list$z, g_k, data_list$z * g_k)
      colnames(X.df)[(ncol(X.df) - 2):ncol(X.df)] <- c("Z", "G", "Z.G")
        
    } else {
      ## use same data frame as BART
      
      # set up data frame for stan
      X.df <- data.frame(1, data_list$x[,-1], data_list$z, g_k)
      colnames(X.df)[(ncol(X.df) - 1):ncol(X.df)] <- c("Z", "G")
    }
    
    # sample from lm.mod
    lm.k <- rstan::sampling(
      lm.mod, 
      data = list(
        N = length(y),
        Y = y,
        K = ncol(X.df),
        X = X.df,
        prior_only = FALSE
      ),
      chains = 1, 
      iter = 5000
    )
    
    # estimate DE(g), IE(z,g)
    eff.k[[k]] <- stanLMEffects(
      fit = extract(lm.k), 
      X = data_list$x[,-1], 
      n.post = 1000,
      model = model, 
      g_seq = gseq,
      interact = interact
    )
  }
  
  ## process results
  
  # combine into single matrices
  
  eff.cut <- list()
  
  eff.cut$mu0 <- do.call(rbind, lapply(eff.k, function(res){res$mu0}))
  eff.cut$mu1 <- do.call(rbind, lapply(eff.k, function(res){res$mu1}))
  eff.cut$DE <- do.call(rbind, lapply(eff.k, function(res){res$DE}))
  eff.cut$IE0 <- do.call(rbind, lapply(eff.k, function(res){res$IE0}))
  eff.cut$IE1 <- do.call(rbind, lapply(eff.k, function(res){res$IE1}))
  
  # save all results
  results <- list(
    theta_sim = theta_0,
    plugin = eff.plugin, 
    cut = eff.cut
  )
  
  # return results
  return(results)
}

stanLMEffects <- function(
    fit, X, n.post, model = 0, 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1)), 
    interact = FALSE
){
  # Estimates causal effects (DE(g), IE(z,g)) using output from a 
  #   Bayesian linear regression model fitted with STAN.
  # Input:
  #   fit: Fitted BART regression (output from 'gbart').
  #   X: Covariate matrix.
  #   n.post: Number of posterior samples.
  #   model: CM model used to generate output.
  #   g_seq: G values at which to estimate (marginal) causal effects.
  #   interact: Boolean indicator of whether outcome model includes Z.G interaction.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1.

  # number of observations
  n.obs <- nrow(X)
  
  # number of monte carlo samples
  n.mcmc <- length(fit$b_Intercept)
  
  # samples kept for effect estimates
  post.samps <- sample(n.mcmc, n.post)
  theta <- cbind(fit$b_Intercept[post.samps], fit$b[post.samps,])
  
  if (interact){
    ## estimate mu(z,g)
    
    # mu(0,g)
    mu0 <- sapply(g_seq, FUN = function(g_val){
      X.m <- cbind(
        rep(1, n.obs), 
        as.matrix(X),
        rep(0, n.obs), 
        rep(g_val, n.obs), 
        rep(0 * g_val, n.obs)
      )
      mu0.g <- X.m %*% t(theta)
      colMeans(mu0.g)
    })
    
    # Z = 1
    mu1 <- sapply(g_seq, FUN = function(g_val){
      X.m <- cbind(
        rep(1, n.obs), 
        as.matrix(X),
        rep(1, n.obs), 
        rep(g_val, n.obs), 
        rep(1 * g_val, n.obs)
      )
      mu1.g <- X.m %*% t(theta)
      colMeans(mu1.g)
    })
    
  } else {
    ## estimate mu(z,g)
    
    # mu(0,g)
    mu0 <- sapply(g_seq, FUN = function(g_val){
      X.m <- cbind(
        rep(1, n.obs), 
        as.matrix(X),
        rep(0, n.obs), 
        rep(g_val, n.obs)
      )
      mu0.g <- X.m %*% t(theta)
      colMeans(mu0.g)
    })
    
    # Z = 1
    mu1 <- sapply(g_seq, FUN = function(g_val){
      X.m <- cbind(
        rep(1, n.obs), 
        as.matrix(X),
        rep(1, n.obs), 
        rep(g_val, n.obs) 
      )
      mu1.g <- X.m %*% t(theta)
      colMeans(mu1.g)
    })
    
  }

  # DE(g) = mu1 - mu0
  DE <- mu1 - mu0
  # IE(0,g)
  IE0 <- ieEst(mu0)
  # IE(1,g)
  IE1 <- ieEst(mu1)
  
  # return results
  results <- list(
    mu0 = mu0,
    mu1 = mu1,
    DE = DE,
    IE0 = IE0,
    IE1 = IE1
  )
  
  return(results)
}


#####################################################################
### 3. Poisson Outcome Models
#####################################################################

simPM1 <- function(
  log_pop, X, Z, G, 
  params = list(
    beta_0 = 1, 
    beta_vec = c(1.5, -0.5, 0.05, 5, -0.0001, -0.3),
    tau = -0.25, 
    gamma = -1.25
  )
) {
  # Simulates Poisson outcomes (i.e., 'PM1') from a given set of covariates 
  #   and parameter values.
  # Input:
  #   log_pop: Vector with log(ZCTA population) values, used as an offset.
  #   X: Matrix of covariates.
  #   Z: Vector of Z (direct treatment) values.
  #   G: Vector of G (indirect treatment) values.
  #   params: List of parameters used to simulate outcomes (default = PM1).
  # Output:
  #   A vector of simulated outcomes from PM1.
  
  # rate function
  mu <- exp(log_pop + params$beta_0 + X %*% params$beta_vec + params$tau * Z + params$gamma * G)
  # sample from Pois(mu_i)
  stats::rpois(n = length(mu), lambda = mu)
}

xMatPM2 <- function(X, z, g){
  # Creates the design matrix used in PM2.
  # Input:
  #   X: Matrix of covariates.
  #   z: Vector of Z (direct treatment) values.
  #   g: Vector of G (indirect treatment) values.
  # Output:
  #   Design matrix used to simulate outcomes from PM2.

  n_obs <- nrow(X)
  
  intercept <- rep(1, n_obs) # intercept
  x1 <- 1 / (1 + exp(-5 * X[,1])) # black proportion
  x2 <- X[,2]
  x3 <- X[,3]
  x4 <- X[,4]^2
  x5 <- X[,4]
  x6 <- exp(-X[,5] / 10000)
  x7 <- X[,6]
  x8 <- z
  x9 <- z * X[,1]
  x10 <- g
  x11 <- g * X[,1] * X[,4]
  
  X_v2 <- cbind(
    intercept,
    x1, x2, x3, x4, x5, 
    x6, x7, x8, x9, x10, x11
  )
  
  return(X_v2)
}

simPM2 <- function(
  log_pop, X, Z, G,
  params = c(0.75, 1.5, -0.5, 0.05, 4, 2, -1, 0.1, -0.25, -1.5, -0.5, -0.75)
) {
  # Simulates Poisson outcomes (i.e., 'PM2') from a given set of covariates 
  #   and parameter values.
  # Input:
  #   log_pop: Vector with log(ZCTA population) values, used as an offset.
  #   X: Matrix of covariates.
  #   Z: Vector of Z (direct treatment) values.
  #   G: Vector of G (indirect treatment) values.
  #   params: Vector of parameters used to generate outcomes
  #     (default = values used in the paper for PM2).
  # Output:
  #   A vector of simulated outcomes from PM2.
  
  # Create design matrix
  X_new <- xMatPM2(X, Z, G)
  # log(rate)
  log_rate <- c(X_new %*% params)
  # rate (accounts for offset)
  mu <- exp(log_pop + log_rate)
  # simulate outcomes
  stats::rpois(n = length(mu), lambda = mu)
}

simPM3 <- function(
    log_pop, X, Z, G, 
    params = c(0.050, -0.025, 0.325, -0.275, 0.575, 0.325, 0.225) 
  ) {
  # Simulates Poisson outcomes (i.e., 'PM3') from a given set of covariates 
  #   and parameter values.
  # Input:
  #   log_pop: Vector with log(ZCTA population) values, used as an offset.
  #   X: Matrix of covariates.
  #   Z: Vector of Z (direct treatment) values.
  #   G: Vector of G (indirect treatment) values.
  #   params: Vector of "true" parameter values 
  #     (default = values used in the paper for PM3).
  # Output:
  #   A vector of simulated outcomes from PM3.

  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  center_vals <- attr(X_new, which = 'scaled:center')
  scale_vals <- attr(X_new, which = 'scaled:scale')  
  
  # stratification
  g1 <- (G < 0.65)
  g3 <- (G >= 0.75)
  g2 <- !g1 & !g3
  
  # X * beta
  X_beta <- X_new %*% params
  
  # mu1
  log.r1 <- X_beta - 0.3 * Z # (Z * X_new) %*% params_v2 - 0.3 * Z
  log.r1[G < 0.4] <- log.r1[G < 0.4] + 1.7
  
  # mu2
  Xb_mod <- X_beta
  Xb_mod[X_beta < -1] <- -1  
  log.r2 <- exp(-Xb_mod) - (0.25 * Z) - 0.1 * (1 + Z) * sin(40 * (G - 0.7))
  
  # mu3
  log.r3 <- X_beta - (0.5 * Z) - 0.5 * exp(2 * G)
  
  # combine into a single mean function
  log.rate <- rep(0, nrow(X))
  log.rate[g1] <- log.r1[g1]
  log.rate[g2] <- log.r2[g2]
  log.rate[g3] <- log.r3[g3]
  log.rate <- log.rate - 5
  
  # sample from Poisson(mu_{0i} * rate_i)
  mu <- exp(log_pop + log.rate)
  y_vals  <- stats::rpois(n = length(mu), lambda = mu)
  
  return(y_vals)
}

PM_SimStudy_Pois <- function(
    ss_num, n_gsamples, data_list, theta_posterior, 
    param_ls, model = 1
) {
  # Performs a simulation study with Bayesian Poisson regression. This includes
  #   simulation of the outcomes, and inference both with and without 
  #   incorporation of uncertainty in G.
  # Input:
  #   ss_num: A scalar value for setting the seed of the random number generator.
  #   n_gsamples: The number of independent MCMC samplers to run when 
  #       incorporating uncertainty in G.
  #   data_list: A list with necessary data for the simulation study 
  #       (so4_mats, em_mat, proj_mat, s_vec, key_assoc, x, z, g_mean).
  #   theta_posterior: A copula model, representing the theta posterior from 
  #       an analysis of the SO4 data.
  #   param_ls: Parameters used to generate PM outcomes.
  #   model: Specifies the Poisson outcome model to simulate from 
  #     (i.e., 1 = PM1, 2 = PM2, 3 = PM3; default = 1).
  # Output: 
  #   A list with i) the 'true' theta vector, ii) MCMC samples using 
  #     G = E(G | *), and iii) MCMC samples incorporating uncertainty in G.
  
  # set seed
  set.seed(ss_num)
  
  ## simulate outcome
  
  # sample theta
  s0 <- copula::rMvdc(1, theta_posterior)
  theta_0 <- c(s0[1:2], 2.6, 49600, 50, s0[3], 0.65)
  
  # generate G
  g_0 <- gCalc(
    theta = theta_0,
    mats = data_list$so4_mats,
    X_mat = data_list$em_mat,
    P_mat = data_list$proj_mat,
    s_vec = data_list$s_vec,
    key = data_list$key_assoc, 
    version = 3
  )
  
  # simulate outcomes
  if (model == 1){
    # simulate from PM1
    y <- simPM1(
      log_pop = data_list$x[, 1],
      X = as.matrix(data_list$x[, -1]),
      Z = data_list$z, G = g_0,
      params = param_ls
    )
  } else if (model == 2){
    # simulate from PM2
    y <- simPM2(
      log_pop = data_list$x[, 1],
      X = as.matrix(data_list$x[, -1]),
      Z = data_list$z, G = g_0,
      params = param_ls
    )
  } else if (model == 3){
    # simulate from PM3
    y <- simPM3(
      log_pop = data_list$x[, 1],
      X = as.matrix(data_list$x[, -1]),
      Z = data_list$z, G = g_0,
      params = param_ls[[1]]
    )
  }
  
  ## estimate model without uncertainty
  
  # create data frame with x, z, and g
  X_mean <- data.frame(y, data_list$x, data_list$z, data_list$g_mean)
  # use log_pop as offset
  off <- X_mean$log_pop
  # fit GLM to initialize MCMC
  glm_mean <- stats::glm(
    y ~ . + offset(log_pop) - log_pop,
    family = stats::poisson, data = X_mean
  )
    
  # mcmc results with no uncertainty in G
  mcmc_meanG <- pois_reg_cpp(
    y_vector = y,
    X_matrix = as.matrix(X_mean[, -c(1:2)]),
    offset = off,
    theta_0 = glm_mean$coefficients,
    n_mcmc = 750000,
    thin = 100,
    burnin = 250000,
    n_adapt = 5000,
    Sigma_0 = tuningMat(),
    keep_burnin = FALSE,
    prior_var = 1
  )
  
  eff.plugin <- poisEffects(
    fit = mcmc_meanG[[1]], 
    X = as.matrix(data_list$x[,-1]),
    n_post = 1000
  )
    
  ## Estimate model WITH uncertainty
  
  # save results in a list
  effects_k <- list()
  
  for (k in 1:n_gsamples) {
    
    cat(paste("\n Iteration ", k, "\n", sep = ""))
    
    # sample theta from posterior
    post_sample <- copula::rMvdc(1, theta_posterior)
    theta_k <- c(post_sample[1:2], 2.6, 49600, 50, post_sample[3], 0.65)
    
    cat("   theta: ")
    cat(theta_k)
    
    # create G
    g_k <- gCalc(
      theta = theta_k,
      mats = data_list$so4_mats,
      X_mat = data_list$em_mat, 
      P_mat = data_list$proj_mat,
      s_vec = data_list$s_vec, 
      key = data_list$key_assoc, 
      version = 3
    )
    
    # create data frame with x, z, and g
    X_df <- data.frame(y, data_list$x, data_list$z, g_k)
    # use log_pop as offset
    off <- X_df$log_pop
      
    # fit GLM to initialize MCMC
    glm_fit_k <- stats::glm(
      y ~ . + offset(log_pop) - log_pop, 
      family = stats::poisson, data = X_df
    )
      
    # generate MCMC samples
    g_uncertainty <- pois_reg_cpp(
      y_vector = y,
      X_matrix = as.matrix(X_df[, -c(1:2)]),
      offset = off,
      theta_0 = glm_fit_k$coefficients,
      n_mcmc = 750000,
      thin = 100,
      burnin = 250000,
      n_adapt = 5000,
      Sigma_0 = tuningMat(),
      keep_burnin = FALSE,
      prior_var = 1
    )
  
    effects_k[[k]] <- poisEffects(
      fit = g_uncertainty[[1]], 
      X = as.matrix(data_list$x[,-1]),
      n_post = 100
    )
  }
  
  ## process results
  
  # combine into single matrices
  
  eff.cut <- list()
  
  eff.cut$mu0 <- do.call(rbind, lapply(effects_k, function(res){res$mu0}))
  eff.cut$mu1 <- do.call(rbind, lapply(effects_k, function(res){res$mu1}))
  eff.cut$DE <- do.call(rbind, lapply(effects_k, function(res){res$DE}))
  eff.cut$IE0 <- do.call(rbind, lapply(effects_k, function(res){res$IE0}))
  eff.cut$IE1 <- do.call(rbind, lapply(effects_k, function(res){res$IE1}))
  
  # save all results
  results <- list(
    theta_sim = theta_0,
    plugin = eff.plugin, 
    cut = eff.cut
  )
  
  # return results
  return(results)
}

poisEffects <- function(
  fit, X, n_post, 
  g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
){
  # Estimates causal effects (DE(g), IE(z,g)) from a fitted Poisson 
  #   regression model.
  # Input:
  #   fit: Fitted Poisson regression (output from 'pois_reg_cpp').
  #   X: Covariate matrix.
  #   n_post: Number of posterior samples.
  #   g_seq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1
  
  # subsample from posterior, if necessary
  if (n_post < nrow(fit)){
    samps <- sample(nrow(fit), n_post)
  } else {
    samps <- 1:nrow(fit)
  }
  
  # add intercept
  X.m <- cbind(int = rep(1, nrow(X)), X)
  
  # number of g values 
  n_g <- length(g_seq)
  
  # estimate mu0
  z <- 0
  mu0 <- matrix(NA_real_, nrow = n_post, ncol = n_g)
  for (k in 1:n_g){
    g_val <- g_seq[k]
    X.full <- cbind(X.m, z, g_val)
    log.rate <- X.full %*% t(fit[samps,])
    mu0[,k] <- colMeans(exp(log.rate))
  }
  
  # estimate mu1
  z <- 1
  mu1 <- matrix(NA_real_, nrow = n_post, ncol = n_g)
  for (k in 1:n_g){
    g_val <- g_seq[k]
    X.full <- cbind(X.m, z, g_val)
    log.rate <- X.full %*% t(fit[samps,])
    mu1[,k] <- colMeans(exp(log.rate))
  }
  
  # DE(g) = mu(1,g) - mu(0,g)
  DE <- mu1 - mu0
  # IE(0,g)
  IE0 <- ieEst(mu0)
  # IE(1,g)
  IE1 <- ieEst(mu1)
  
  # return results
  results <- list(
    mu0 = mu0,
    mu1 = mu1,
    DE = DE,
    IE0 = IE0,
    IE1 = IE1
  )
  
  return(results)
}

PM_SimStudy_BART <- function(
    ss_num, n_gsamples, data_list, theta_posterior, burnin, param_ls, model = 1
) {
  # Performs a simulation study with log-linear BART regression. This includes
  #   simulation of the outcomes, and inference both with and without 
  #   incorporation of uncertainty in G.
  # Input:
  #   ss_num: A scalar value for setting the seed of the random number generator.
  #   n_gsamples: The number of independent MCMC samplers to run when 
  #       incorporating uncertainty in G.
  #   data_list: A list with necessary data for the simulation study 
  #       (so4_mats, em_mat, proj_mat, s_vec, key_assoc, x, z, g_mean).
  #   theta_posterior: A copula model, representing the theta posterior from 
  #       an analysis of the SO4 data.
  #   burnin: Burnin period for BART model.
  #   param_ls: Parameters used to generate PM outcomes.
  #   model: Specifies the Poisson outcome model to simulate from 
  #     (i.e., 1 = PM1, 2 = PM2, 3 = PM3; default = 1).
  # Output: 
  #   A list with i) the 'true' theta vector, ii) MCMC samples using 
  #     G = E(G | *), and iii) MCMC samples incorporating uncertainty in G.
  
  # set seed
  set.seed(ss_num)
  
  ## simulate outcome
  
  # sample theta
  s0 <- copula::rMvdc(1, theta_posterior)
  theta_0 <- c(s0[1:2], 2.6, 49600, 50, s0[3], 0.65)
  
  # generate G
  g_0 <- gCalc(
    theta = theta_0,
    mats = data_list$so4_mats,
    X_mat = data_list$em_mat,
    P_mat = data_list$proj_mat,
    s_vec = data_list$s_vec,
    key = data_list$key_assoc, 
    version = 3
  )
  
  # simulate outcomes
  if (model == 1){
    # simulate from PM1
    y <- simPM1(
      log_pop = data_list$x[, 1],
      X = as.matrix(data_list$x[, -1]),
      Z = data_list$z, G = g_0,
      params = param_ls
    )
  } else if (model == 2){
    # simulate from PM2
    y <- simPM2(
      log_pop = data_list$x[, 1],
      X = as.matrix(data_list$x[, -1]),
      Z = data_list$z, G = g_0,
      params = param_ls
    )
  } else if (model == 3){
    # simulate from PM3
    y <- simPM3(
      log_pop = data_list$x[, 1],
      X = as.matrix(data_list$x[, -1]),
      Z = data_list$z, G = g_0,
      params = param_ls[[1]]
    )
  }
  
  ## estimate model without uncertainty
  
  # create data frame with x, z, and g
  X_mean <- data.frame(y, data_list$x, data_list$z, data_list$g_mean)
  colnames(X_mean)[c(ncol(X_mean) - 1, ncol(X_mean))] <- c("Z", "G")
  # offset 
  offset_vec <- exp(X_mean$log_pop)
  
  # determine leaf prior hyperparameter
  rates <- y / offset_vec
  r_star <- quantile(rates, probs = 0.975)
  r_mean <- mean(rates)
  a_0 <- 0.5 * (log(r_star) - log(r_mean))
  
  # fit model 
  bart.plugin <- count_bart(
    y = y,
    x = as.matrix(X_mean[,-c(1)]),
    offset = exp(log(offset_vec) + log(r_mean)),
    nburn = burnin,
    nsim = 500,
    nthin = 5,
    update_interval = 100,
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
  
  # estimate tau, gamma0, gamma1
  eff.plugin <- bartPoisEffects(
    fit = bart.plugin, 
    X = data_list$x, 
    r_mean = r_mean, 
    n_post = 500
  )

  ## Estimate model WITH uncertainty
  
  # save results in a list
  effects_k <- list()
  
  for (k in 1:n_gsamples) {
    cat(paste("\n Iteration ", k, "\n", sep = ""))
    
    # sample theta from posterior
    post_sample <- copula::rMvdc(1, theta_posterior)
    theta_k <- c(post_sample[1:2], 2.6, 49600, 50, post_sample[3], 0.65)
    
    cat("   theta: ")
    cat(theta_k)
    
    # create G
    g_k <- gCalc(
      theta = theta_k,
      mats = data_list$so4_mats,
      X_mat = data_list$em_mat,
      P_mat = data_list$proj_mat,
      s_vec = data_list$s_vec,
      key = data_list$key_assoc,
      version = 3
    )
    
    # create data frame with x, z, and g
    X_df <- data.frame(y, data_list$x, data_list$z, g_k)
    colnames(X_df)[c(ncol(X_df) - 1, ncol(X_df))] <- c("Z", "G")
    # offset
    offset_vec <- exp(X_df$log_pop)
    
    # determine leaf hyperparameter value
    rates <- y / offset_vec
    r_star <- quantile(rates, probs = 0.975)
    r_mean <- mean(rates)
    a_0 <- 0.5 * (log(r_star) - log(r_mean))
    
    # generate MCMC samples
    bart_k <- count_bart(
      y = y,
      x = as.matrix(X_df[, -1]),
      offset = exp(log(offset_vec) + log(r_mean)),
      nburn = burnin,
      nsim = 100,
      nthin = 5,
      update_interval = 100,
      ntree = 200,
      a0 = a_0,
      sd = 2 * sd(y),
      base = 0.95,
      power = 2,
      nu = 3,
      lambda = NULL,
      sigq = .9,
      sighat = NULL, 
      kappa_a = 5, kappa_b = 3,
      count_model = "poisson", 
      debug = FALSE,
      # used as input... but not relevant for count models...?
      randeff_design = matrix(1),
      randeff_variance_component_design = matrix(1),
      randeff_scales = 1,
      randeff_df = 3
    )
    
    # estimate tau, gamma0, gamma1
    effects_k[[k]] <- bartPoisEffects(
      fit = bart_k,
      X = data_list$x,
      r_mean = r_mean,
      n_post = 100
    )
  }
  
  ## process results
  
  # combine into single matrices
  
  eff.cut <- list()
  
  eff.cut$mu0 <- do.call(rbind, lapply(effects_k, function(res){res$mu0}))
  eff.cut$mu1 <- do.call(rbind, lapply(effects_k, function(res){res$mu1}))
  eff.cut$DE <- do.call(rbind, lapply(effects_k, function(res){res$DE}))
  eff.cut$IE0 <- do.call(rbind, lapply(effects_k, function(res){res$IE0}))
  eff.cut$IE1 <- do.call(rbind, lapply(effects_k, function(res){res$IE1}))
  
  # save all results
  results <- list(
    theta_sim = theta_0,
    plugin = eff.plugin, 
    cut = eff.cut
  )
  
  # return results
  return(results)
}

bartPoisEffects <- function(
    fit, X, r_mean, n_post, 
    g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
){
  # Estimates casual effects using output from log-linear BART.
  # Input:
  #   fit: Fitted log-linear BART model (from 'count_bart').
  #   X: Covariate matrix.
  #   r_mean: Rate standardization value used in log-linear BART.
  #   n_post: Number of posterior samples to generate.
  #   g_seq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1.
  
  # number of g values used in prediction
  n_g <- length(g_seq)

  # estimate mu0
  z <- 0
  mu0 <- matrix(NA_real_, nrow = n_post, ncol = n_g)
  for (k in 1:n_g){
    g_val <- g_seq[k]
    X_m <- cbind(X, z, g_val)
    colnames(X_m)[(ncol(X_m) - 1):ncol(X_m)] <- c("Z", "G")
    out <- fit$tree_fit$tree_samples$coefs(as.matrix(t(X_m)))
    f_hat <- exp(out[1, , ] + log(r_mean))
    mu0[, k] <- colMeans(f_hat)
  }

  # estimate mu1
  mu1 <- matrix(NA_real_, nrow = n_post, ncol = n_g)
  z <- 1
  for (k in 1:n_g){
    g_val <- g_seq[k]
    X_m <- cbind(X, z, g_val)
    colnames(X_m)[(ncol(X_m) - 1):ncol(X_m)] <- c("Z", "G")
    out <- fit$tree_fit$tree_samples$coefs(as.matrix(t(X_m)))
    f_hat <- exp(out[1, , ] + log(r_mean))
    mu1[, k] <- colMeans(f_hat)
  }

  # DE(g) = mu(1,g) - mu(0,g)
  DE <- mu1 - mu0
  # IE(0,g)
  IE0 <- ieEst(mu0)
  # IE(1,g)
  IE1 <- ieEst(mu1)

  # return results
  results <- list(
    mu0 = mu0,
    mu1 = mu1,
    DE = DE,
    IE0 = IE0,
    IE1 = IE1
  )
  
  return(results)
}


#####################################################################
### 4. Processing the results
#####################################################################

trueCM1 <- function(X, params, g_seq = seq(from = 0.25, to = 0.9, by = 0.01)){
  # Generates the "true" causal effects from CM1.
  # Input:
  #   X: covariate matrix.
  #   params: parameters used in the simulation study.
  #   g_seq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1.
  
  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  n.obs <- nrow(X_new)
  
  ### estimate (mu0, mu1) for different levels of (z,g)
  
  # Z = 0
  mu0 <- sapply(g_seq, FUN = function(g_val){
    X_m <- cbind(X_new, rep(0, n.obs), rep(g_val, n.obs), rep(0 * g_val, n.obs))
    mu0.g <- X_m %*% params
    mean(mu0.g)
  })
  
  # Z = 1
  mu1 <- sapply(g_seq, FUN = function(g_val){
    X_m <- cbind(X_new, rep(1, n.obs), rep(g_val, n.obs), rep(1 * g_val, n.obs))
    mu1.g <- X_m %*% params
    mean(mu1.g)
  })
  
  list(
    mu0 = mu0,
    mu1 = mu1,
    DE = mu1 - mu0,
    IE0 = mu0 - mu0[1],
    IE1 = mu1 - mu1[1]
  )
}

trueCM2 <- function(X, params, g_seq = seq(from = 0.25, to = 0.9, by = 0.01)){
  # Generates the "true" causal effects from CM2.
  # Input:
  #   X: covariate matrix.
  #   params: parameters used in the simulation study.
  #   g_seq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1.
  
  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  W <- matrix(0.5, nrow = nrow(X_new), ncol = ncol(X_new))
  center_vals <- attr(X_new, which = 'scaled:center')
  scale_vals <- attr(X_new, which = 'scaled:scale')
  X_beta <- X_new %*% params
  
  ### estimate (mu0, mu1) for different levels of (z,g)
  
  # Z = 0
  mu0 <- sapply(g_seq, FUN = function(g_val){
    XW_beta <- (X_new - W * g_val) %*% params
    mu0_g <- exp(XW_beta)
    mean(mu0_g)
  })
  
  mu1 <- sapply(g_seq, FUN = function(g_val){
    mu1_g <- X_beta - g_val^2
    mean(mu1_g)
  })
  
  list(
    mu0 = mu0,
    mu1 = mu1,
    DE = mu1 - mu0,
    IE0 = mu0 - mu0[1],
    IE1 = mu1 - mu1[1]
  )
}

meanCM3 <- function(
    X, Z, G, 
    params = c(0.2, -0.1, 1.3, -1.1, 2.3, 1.3, 0.9),
    params_v2 = c(-1.5, -0.7, -0.3, 0.4, -0.5, 2.0, 1.0)
){
  # Mean function used in CM3.
  # Input:
  #   X: covariate matrix.
  #   Z: direct treatment values.
  #   G: indirect treatment values.
  #   params: first set of parameters used in CM3.
  #   params_v2: second set of parameters used in CM3.
  # Output:
  #   Mean function given X, Z, G, and params.
  
  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  center_vals <- attr(scale(X), which = 'scaled:center')
  scale_vals <- attr(scale(X), which = 'scaled:scale')  
  
  # stratification
  g1 <- (G < 0.65)
  g3 <- (G >= 0.75)
  g2 <- !g1 & !g3
  
  # X * beta
  X_beta <- X_new %*% params
  
  # mu1
  mu1 <- X_beta + (Z * X_new) %*% params_v2 - 0.3 * Z
  mu1[G < 0.4] <- mu1[G < 0.4] + 1.7
  
  # mu2
  Xb_mod <- X_beta
  Xb_mod[X_beta < -3] <- -3  
  mu2 <- exp(-Xb_mod) - (0.25 * Z) - (1 + Z) * sin(40 * (G - 0.7))
  
  # mu3
  mu3 <- X_beta - (0.5 * Z) - 0.5 * exp(2 * G)
  
  # combine into a single mean function
  mu <- rep(0, nrow(X))
  mu[g1] <- mu1[g1]
  mu[g2] <- mu2[g2]
  mu[g3] <- mu3[g3]
  
  return(mu)
}

trueCM3 <- function(X, g_seq = seq(from = 0.25, to = 0.9, by = 0.01)){
  # Generates the "true" causal effects from CM3.
  # Input:
  #   X: covariate matrix
  #   g_seq: G values at which to estimate (marginal) causal effects.
  # Output:
  #   A list containing mu0, mu1, DE, IE0, and IE1.
  
  ### estimate (mu0, mu1) for different levels of (z,g)
  
  # Z = 0
  mu0 <- sapply(g_seq, FUN = function(g_val){
    mu0.g <- meanCM3(
      X = X, Z = rep(0, nrow(X)), G = rep(g_val, nrow(X)) 
    )
    mean(mu0.g)
  })
  
  # Z = 1
  mu1 <- sapply(g_seq, FUN = function(g_val){
    mu1.g <- meanCM3(
      X = X, Z = rep(1, nrow(X)), G = rep(g_val, nrow(X)) 
    )
    mean(mu1.g)
  })
  
  # Results
  list(
    mu0 = mu0,
    mu1 = mu1,
    DE = mu1 - mu0,
    IE0 = mu0 - mu0[1],
    IE1 = mu1 - mu1[1]
  )
}

biasCalc <- function(results, te, loglin = TRUE){
  # Estimate the bias of the causal effect estimates.
  # Input:
  #   results: simulation study results
  #   te: "true effect" for a given simulation study
  #   loglin: Boolean indicator if the data are from a continuous
  #     or count (e.g., log-linear BART) study.
  # Output:
  #   Matrix with estimated bias for DE, IE0, and IE1.
  
  bias.plugin <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  bias.cut <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  if (loglin){
    mult <- 1000
  } else {
    mult <- 1
  }
  
  for (k in 1:length(results)){
    
    bias.plugin$DE[k,] <- colMeans(mult * results[[k]]$plugin$DE) - te$DE
    bias.plugin$IE0[k,] <- colMeans(mult * results[[k]]$plugin$IE0) - te$IE0
    bias.plugin$IE1[k,] <- colMeans(mult * results[[k]]$plugin$IE1) - te$IE1
    
    bias.cut$DE[k,] <- colMeans(mult * results[[k]]$cut$DE) - te$DE
    bias.cut$IE0[k,] <- colMeans(mult * results[[k]]$cut$IE0) - te$IE0
    bias.cut$IE1[k,] <- colMeans(mult * results[[k]]$cut$IE1) - te$IE1
  }
  
  bias <- list(
    plugin = cbind(
      colMeans(bias.plugin$DE),
      colMeans(bias.plugin$IE0),
      colMeans(bias.plugin$IE1)
    ),
    cut = cbind(
      colMeans(bias.cut$DE),
      colMeans(bias.cut$IE0),
      colMeans(bias.cut$IE1)
    )
  )
  
  colnames(bias$plugin) <- c("DE", "IE0", "IE1")
  colnames(bias$cut) <- c("DE", "IE0", "IE1")
  
  return(bias)
}

coverageCalc <- function(results, te, loglin = TRUE){
  # Estimates the 95% credible interval coverage rates.
  # Input:
  #   results: simulation study results
  #   te: "true effect" for a given simulation study
  #   loglin: Boolean indicator if the data are from a continuous
  #     or count (e.g., log-linear BART) study
  # Output:
  #   Matrix with credible interval coverage rates for DE, IE0, and IE1
  
  cover.plugin <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  cover.cut <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  if (loglin){
    mult <- 1000
  } else {
    mult <- 1
  }
  
  for (k in 1:length(results)){
    cover.plugin$DE[k,] <- coverageCheck(mult * results[[k]]$plugin$DE, te$DE)
    cover.plugin$IE0[k,] <- coverageCheck(mult * results[[k]]$plugin$IE0, te$IE0) 
    cover.plugin$IE1[k,] <- coverageCheck(mult * results[[k]]$plugin$IE1, te$IE1)
    
    cover.cut$DE[k,] <- coverageCheck(mult * results[[k]]$cut$DE, te$DE)
    cover.cut$IE0[k,] <- coverageCheck(mult * results[[k]]$cut$IE0, te$IE0) 
    cover.cut$IE1[k,] <- coverageCheck(mult * results[[k]]$cut$IE1, te$IE1)
  }
  
  coverage <- list(
    plugin = cbind(
      colMeans(cover.plugin$DE),
      colMeans(cover.plugin$IE0),
      colMeans(cover.plugin$IE1)
    ),
    cut = cbind(
      colMeans(cover.cut$DE),
      colMeans(cover.cut$IE0),
      colMeans(cover.cut$IE1)
    )
  )
  
  colnames(coverage$plugin) <- c("DE", "IE0", "IE1")
  colnames(coverage$cut) <- c("DE", "IE0", "IE1")
  
  return(coverage)
}

coverageCheck <- function(res, true){
  # Check if 'true' causal effect is within 95% credible interval.
  # Input:
  #   res: results vector
  #   true: true effect
  # Output:
  #   Percent of estimates with 95% credible intervals containing true param
  quants <- rowMeans(apply(res, 1, FUN = function(est){true < est}))
  (quants >= 0.025 & quants <= 0.975)
}

varComponents <- function(res, n_samps){
  # Compute the between, within, and total variance of an MI estimate.
  # Input:
  #   res: results
  #   n_samps: number of multiple imputation samples
  # Output:
  #   List containing total (T), within (U), and between (B) variance.
  
  n <- nrow(res)
  samp_idx <- n / n_samps
  start_idx <- seq(from = 1, to = n, by = samp_idx)
  stop_idx <- seq(from = samp_idx, to = n, by = samp_idx)
  
  u_vals <- matrix(NA_real_, nrow = n_samps, ncol = ncol(res))
  q_vals <- matrix(NA_real_, nrow = n_samps, ncol = ncol(res))
  q_bar <- colMeans(res)
  
  for (k in 1:n_samps){
    idx <- start_idx[k]:stop_idx[k]
    u_vals[k,] <- apply(res[idx,], 2, var)
    q_vals[k,] <- colMeans(res[idx,])
  }
  
  u_bar <- colMeans(u_vals)
  b <- apply(q_vals, 2, var)
  total <- u_bar + (1 + 1/n_samps) * b
  
  list(
    T = total, 
    U = u_bar, 
    B = b
  )
}

varianceCalc <- function(results, loglin = TRUE){
  # Estimate the within, between, and total variance for a given set of 
  #   causal effect estimates.
  # Input:
  #   results: results from a simulation study
  #   loglin: Boolean indicator if the data are from a continuous
  #     or count (e.g., log-linear BART) study
  # Output:
  #   List with DE, IE0, and IE1 variance estimates.
  
  var_within <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  var_between <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  var_total <- list(
    DE = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE0 = matrix(NA_real_, nrow = length(results), ncol = 8),
    IE1 = matrix(NA_real_, nrow = length(results), ncol = 8)
  )
  
  if (loglin){
    mult <- 1000
  } else {
    mult <- 1
  }
  
  for (s in 1:length(results)){
    var.DE <- varComponents(res = mult * results[[s]]$cut$DE, n_samps = 100)
    var.IE0 <- varComponents(res = mult * results[[s]]$cut$IE0, n_samps = 100)
    var.IE1 <- varComponents(res = mult * results[[s]]$cut$IE1, n_samps = 100)
    
    # DE  
    var_within$DE[s,] <- var.DE$U
    var_between$DE[s,] <- var.DE$B
    var_total$DE[s,] <- var.DE$T
    
    # IE0
    var_within$IE0[s,] <- var.IE0$U
    var_between$IE0[s,] <- var.IE0$B
    var_total$IE0[s,] <- var.IE0$T
    
    # IE1
    var_within$IE1[s,] <- var.IE1$U
    var_between$IE1[s,] <- var.IE1$B
    var_total$IE1[s,] <- var.IE1$T
  }
  
  variance <- list(
    DE = cbind(
      colMeans(var_within$DE),
      colMeans(var_between$DE),
      colMeans(var_total$DE)
    ),
    IE0 = cbind(
      colMeans(var_within$IE0),
      colMeans(var_between$IE0),
      colMeans(var_total$IE0)
    ),
    IE1 = cbind(
      colMeans(var_within$IE1),
      colMeans(var_between$IE1),
      colMeans(var_total$IE1)
    )
  )
  
  colnames(variance$DE) <- c("Within", "Between", "Total")
  colnames(variance$IE0) <- c("Within", "Between", "Total")
  colnames(variance$IE1) <- c("Within", "Between", "Total")
  
  return(variance)
}


#####################################################################
### 5. Plotting the results
#####################################################################

plotSimStudyEst <- function(no_unc, unc, true, title = "", y_title = "", loglin = TRUE){
  # Plots the estimated causal effect curves with and without uncertainty
  #   propagation.
  # Input:
  #   no_unc: results without uncertainty propagation
  #   unc: results with uncertainty propagation
  #   true: true effects
  #   title: main figure title
  #   y_title: y-axis title
  #   loglin: indicates if results are from count data simulation
  # Output:
  #   Plot with CE curves, both with and without uncertainty propagation
  
  col1 <- "#04647E"
  col2 <- "#8C4646"
  col_true <- "red"
        
  # c("#588C7E", "#F2E394", "#F2AE72", "#D96459", "#8C4646")
  g_seq <- c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
  # c(0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  g_seq2 <- seq(from = 0.25, to = 0.9, by = 0.01)
      
      
  if (loglin){
    mu_no <- colMeans(1000 * no_unc)
    ub_no <- apply(1000 * no_unc, 2, quantile, probs = 0.975)
    lb_no <- apply(1000 * no_unc, 2, quantile, probs = 0.025)
    
    mu_unc <- colMeans(1000 * unc)
    ub_unc <- apply(1000 * unc, 2, quantile, probs = 0.975)
    lb_unc <- apply(1000 * unc, 2, quantile, probs = 0.025)
  } else {
    mu_no <- colMeans(no_unc)
    ub_no <- apply(no_unc, 2, quantile, probs = 0.975)
    lb_no <- apply(no_unc, 2, quantile, probs = 0.025)
    
    mu_unc <- colMeans(unc)
    ub_unc <- apply(unc, 2, quantile, probs = 0.975)
    lb_unc <- apply(unc, 2, quantile, probs = 0.025)
  }
      
  y_lims <- c(min(lb_no, lb_unc, true), max(ub_no, ub_unc, true))
      
  # create plot
  par(
    mfrow = c(1,2), 
    mar = c(4,4,2,1),
    oma = c(1,1,2,1)
  )
      
  # plug-in plot
  plot(g_seq, mu_no, type = "l", ylim = y_lims, bty = "l",
       xlab = "g", ylab = y_title, col = col1, lwd = 2)
  polygon(
    c(g_seq, rev(g_seq)), c(ub_no, rev(lb_no)),
    col = alpha(col1, alpha = 0.2), border = FALSE
  )
  lines(g_seq2, true, col = col_true, lty = 3, lwd = 2)
  title(main = "A.  plug-in", cex.main = 1, line = 0, adj = 0.05)
      
  # cut feedback plot
  plot(g_seq, mu_unc, type = "l", ylim = y_lims, bty = "l",
       xlab = "g", ylab = y_title, col = col2, lwd = 2)
  polygon(
    c(g_seq, rev(g_seq)), c(ub_unc, rev(lb_unc)),
    col = alpha(col2, alpha = 0.2), border = FALSE
  )
  lines(g_seq2, true, col = col_true, lty = 3, lwd = 2)
  title(main = "B.  cut feedback", cex.main = 1, line = 0, adj = 0.05)
    
  title(main = title, line = 0, adj = 0.05, outer = TRUE)
}

truePM1 <- function(
    X, params = c(1, 1.5, -0.5, 0.05, 5, -0.0001, -0.3, -0.25, -1.25),
    g_seq = seq(from = 0.25, to = 0.9, by = 0.01)
){
  # Generates the "true" causal effect for a sim study PM1.
  # Input:
  #   X: covariate matrix
  #   params: parameters used in the simulation study
  #   g_seq: g values at which to return the true effect
  # Output:
  #   List with the true mu0, mu1, DE, IE0, and IE1 values
  
  n.obs <- nrow(X)
  
  ### estimate (mu0, mu1) for different levels of (z,g)
  
  # Z = 0
  mu0 <- sapply(g_seq, FUN = function(g_val){
    X.m <- cbind(rep(1, n.obs), X, rep(0, n.obs), rep(g_val, n.obs))
    mu0.g <- exp(X.m %*% params)
    mean(mu0.g) * 1000
  })
  
  # Z = 1
  mu1 <- sapply(g_seq, FUN = function(g_val){
    X.m <- cbind(rep(1, n.obs), X, rep(1, n.obs), rep(g_val, n.obs))
    mu1.g <- exp(X.m %*% params)
    mean(mu1.g) * 1000
  })
  
  list(
    mu0 = mu0,
    mu1 = mu1,
    DE = mu1 - mu0,
    IE0 = mu0 - mu0[1],
    IE1 = mu1 - mu1[1]
  )
}

truePM2 <- function(
    X, params = c(0.75, 1.5, -0.5, 0.05, 4, 2, -1, 0.1, -0.25, -1.5, -0.5, -0.75),
    g_seq = seq(from = 0.25, to = 0.9, by = 0.01)
){
  # Generates the "true" causal effect for a sim study PM2.
  # Input:
  #   X: covariate matrix
  #   params: parameters used in the simulation study
  #   g_seq: g values at which to return the true effect
  # Output:
  #   List with the true mu0, mu1, DE, IE0, and IE1 values

  ### estimate (mu0, mu1) for different levels of (z,g)
  
  # Z = 0
  mu0 <- sapply(g_seq, FUN = function(g_val){
    X.m <- xMatPM2(X = X, z = 0, g = g_val)
    mu0.g <- exp(X.m %*% params)
    mean(mu0.g) * 1000
  })
  
  # Z = 1
  mu1 <- sapply(g_seq, FUN = function(g_val){
    X.m <- xMatPM2(X = X, z = 1, g = g_val)
    mu1.g <- exp(X.m %*% params)
    mean(mu1.g) * 1000
  })
  
  list(
    mu0 = mu0,
    mu1 = mu1,
    DE = mu1 - mu0,
    IE0 = mu0 - mu0[1],
    IE1 = mu1 - mu1[1]
  )
}

meanPM3 <- function(
    X, Z, G, 
    params = c(0.050, -0.025, 0.325, -0.275, 0.575, 0.325, 0.225)
){
  # Returns the true rate function used in PM3.
  # Input:
  #   X: covariate matrix
  #   Z: z values
  #   G: g values
  #   params: parameters used in the simulation study
  # Output:
  #   Rate function for given Z, G, and parameter values
  
  # standardization
  X_new <- cbind(rep(1, nrow(X)), scale(X))
  center_vals <- attr(X_new, which = 'scaled:center')
  scale_vals <- attr(X_new, which = 'scaled:scale')  
  
  # stratification
  g1 <- (G < 0.65)
  g3 <- (G >= 0.75)
  g2 <- !g1 & !g3
  
  # X * beta
  X_beta <- X_new %*% params
  
  # mu1
  log.r1 <- X_beta - 0.3 * Z # X_beta + (Z * X_new) %*% params_v2 - 0.3 * Z
  log.r1[G < 0.4] <- log.r1[G < 0.4] + 1.7
  
  # mu2
  Xb_mod <- X_beta
  Xb_mod[X_beta < -1] <- -1  
  log.r2 <- exp(-Xb_mod) - (0.25 * Z) - 0.1 * (1 + Z) * sin(40 * (G - 0.7))
  
  # mu3
  log.r3 <- X_beta - (0.5 * Z) - 0.5 * exp(2 * G)
  
  # combine into a single mean function
  log.rate <- rep(0, nrow(X))
  log.rate[g1] <- log.r1[g1]
  log.rate[g2] <- log.r2[g2]
  log.rate[g3] <- log.r3[g3]
  log.rate <- log.rate - 5
  
  # return the rate / 1000
  rate <- exp(log.rate) * 1000
  
  return(rate)
}

truePM3 <- function(X, g_seq = seq(from = 0.25, to = 0.9, by = 0.01)){
  # Generates the "true" causal effect for a sim study PM3.
  # Input:
  #   X: covariate matrix
  #   g_seq: g values at which to return the true effect
  # Output:
  #   List with the true mu0, mu1, DE, IE0, and IE1 values
  
  ### estimate (mu0, mu1) for different levels of (z,g)
  
  # Z = 0
  mu0 <- sapply(g_seq, FUN = function(g_val){
    mu0.g <- meanPM3(
      X = X, Z = rep(0, nrow(X)), G = rep(g_val, nrow(X)) 
    )
    mean(mu0.g) 
  })
  
  # Z = 1
  mu1 <- sapply(g_seq, FUN = function(g_val){
    mu1.g <- meanPM3(
      X = X, Z = rep(1, nrow(X)), G = rep(g_val, nrow(X)) 
    )
    mean(mu1.g) 
  })
  
  list(
    mu0 = mu0,
    mu1 = mu1,
    DE = mu1 - mu0,
    IE0 = mu0 - mu0[1],
    IE1 = mu1 - mu1[1]
  )
}


#####################################################################
### 6. Log-linear BART Sensitivity Analysis
#####################################################################

bartSensitivity <- function(
  n_chains = 50,
  n_tree = 200,
  base_val = 0.95,
  power_val = 2,
  theta_posterior, 
  data_list,
  n_workers = 1
){
  # Performs a sensitivity study with log-linear BART regression. Repeated
  #   log-linear BART MCMC chains are fit to simulated data, using the specified
  #   BART priors. The resulting causal estimates are returned in a list
  # Input:
  #   n_chains: number of independent MCMC chains to fit
  #   n_tree: number of trees to use for log-linear BART prior
  #     (default = 200)
  #   base_val: parameter controlling distribution of terminal nodes 
  #     (default = 0.95)
  #   power_val: parameter controlling distribution of terminal nodes
  #     (default = 2)
  #   theta_posterior: copula model, representing the theta posterior from 
  #       an analysis of the SO4 data.
  #   data_list: list with necessary data for the simulation study 
  #       (adv_mats, em_mat, proj_mat, s_vec, key_assoc, x, z, g_mean).
  #   n_workers: number of chains to fit in parallel
  # Output: 
  #   A list with DE, IE0, and IE1 posterior estimates for each MCMC chain

  # set seed
  set.seed(29)
  
  # sample theta
  s0 <- copula::rMvdc(1, theta_posterior)
  theta_0 <- c(s0[1:2], 2.6, 49600, 50, s0[3], 0.65)
  
  # generate G
  g_0 <- gCalc(
    theta = theta_0,
    mats = data_list$adv_mats,
    X_mat = data_list$em_mat,
    P_mat = data_list$proj_mat,
    s_vec = data_list$s_vec,
    key = data_list$key_assoc, 
    version = 3
  )
  
  # sample y
  y_sim <- simPM2(
    log_pop = data_list$x[,1], 
    X = data_list$x[,-1], 
    Z = data_list$z, 
    G = g_0,
    params = c(0.75, 1.5, -0.5, 0.05, 4, 2, -1, 0.1, -0.25, -1.5, -0.5, -0.75)
  )
  
  # create data frame with x, z, and g
  X_mean <- data.frame(y_sim, data_list$x, data_list$z, data_list$g_mean)
  colnames(X_mean)[1] <- "y"
  colnames(X_mean)[c(ncol(X_mean) - 1, ncol(X_mean))] <- c("Z", "G")
  
  # offset 
  offset_vec <- exp(X_mean$log_pop)
  
  # determine leaf prior hyperparameter
  rates <- y_sim / offset_vec
  r_star <- quantile(rates, probs = 0.975)
  r_mean <- mean(rates)
  a_0 <- 0.5 * (log(r_star) - log(r_mean))
  
  # fit model n_chains times
  ntree.effects <- mclapply(
    1:n_chains, 
    FUN = function(x){
      # fit log-linear BART model
      bf <- count_bart(
        y = y_sim,
        x = as.matrix(X_mean[,-c(1)]),
        offset = exp(log(offset_vec) + log(r_mean)),
        nburn = 10000,
        nsim = 500,
        nthin = 5,
        update_interval = 100,
        ntree = n_tree,
        a0 = a_0,
        sd = 2 * sd(y_sim),
        base = base_val,
        power = power_val,
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
      
      # return causal effect estimates
      bartPoisEffects(
        fit = bf, 
        X = data_list$x, 
        r_mean = r_mean, 
        n_post = 500
      )
    },
    mc.cores = n_workers
  )
  
  # return results
  return(ntree.effects)
}

mtreeStatCalc <- function(m100, m200, m300, m400, X, res.type = "sd"){
  # Estimate results (marginal standard deviation, bias, and RMSE) from
  #   the sensitivity analysis (number of trees).
  # Input:
  #   m100: results from m = 100 analysis
  #   m200: results from m = 200 analysis
  #   m300: results from m = 300 analysis
  #   m400: results from m = 400 analysis
  #   X: covariate matrix
  #   res.type: type of result to calculate ("sd", "bias", "rmse")
  # Output:
  #   List with theta.hat and the associated MCSE.
  
  # number of mcmc chains
  n_chains <- length(m100)
  
  # store results in this data strucure
  mtrees <- list(
    DE = list(
      m100 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m200 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m300 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m400 = matrix(data = NA_real_, nrow = n_chains, ncol = 7) 
    ),
    IE0 = list(
      m100 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m200 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m300 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m400 = matrix(data = NA_real_, nrow = n_chains, ncol = 7) 
    ),
    IE1 = list(
      m100 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m200 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m300 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      m400 = matrix(data = NA_real_, nrow = n_chains, ncol = 7) 
    )
  )
  
  for (k in 1:n_chains){
    # DE(g)
    mtrees$DE$m100[k, ] <- resCalc(m100[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    mtrees$DE$m200[k, ] <- resCalc(m200[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    mtrees$DE$m300[k, ] <- resCalc(m300[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    mtrees$DE$m400[k, ] <- resCalc(m400[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    # IE(z=0,g)
    mtrees$IE0$m100[k, ] <- resCalc(m100[[k]]$IE0[,-1], X.mat = X, 
                                    type = res.type, ce.type = "IE0")
    mtrees$IE0$m200[k, ] <- resCalc(m200[[k]]$IE0[,-1], X.mat = X, 
                                    type = res.type, ce.type = "IE0")
    mtrees$IE0$m300[k, ] <- resCalc(m300[[k]]$IE0[,-1], X.mat = X,
                                    type = res.type, ce.type = "IE0")
    mtrees$IE0$m400[k, ] <- resCalc(m400[[k]]$IE0[,-1], X.mat = X,
                                    type = res.type, ce.type = "IE0")
    # IE(z=1,g)
    mtrees$IE1$m100[k, ] <- resCalc(m100[[k]]$IE1[,-1], X.mat = X, 
                                    type = res.type, ce.type = "IE1")
    mtrees$IE1$m200[k, ] <- resCalc(m200[[k]]$IE1[,-1], X.mat = X, 
                                    type = res.type, ce.type = "IE1")
    mtrees$IE1$m300[k, ] <- resCalc(m300[[k]]$IE1[,-1], X.mat = X, 
                                    type = res.type, ce.type = "IE1")
    mtrees$IE1$m400[k, ] <- resCalc(m400[[k]]$IE1[,-1], X.mat = X, 
                                    type = res.type, ce.type = "IE1")
  }

  mtrees.mean <- list(
    DE = cbind(
      colMeans(mtrees$DE$m100), 
      colMeans(mtrees$DE$m200),  
      colMeans(mtrees$DE$m300),  
      colMeans(mtrees$DE$m400)
    ),
    IE0 = cbind(
      colMeans(mtrees$IE0$m100), 
      colMeans(mtrees$IE0$m200),  
      colMeans(mtrees$IE0$m300),  
      colMeans(mtrees$IE0$m400)
    ),
    IE1 = cbind(
      colMeans(mtrees$IE1$m100), 
      colMeans(mtrees$IE1$m200),  
      colMeans(mtrees$IE1$m300),  
      colMeans(mtrees$IE1$m400)
    )
  )
  
  colnames(mtrees.mean$DE) <- c("m=100", "m=200", "m=300", "m=400")  
  colnames(mtrees.mean$IE0) <- c("m=100", "m=200", "m=300", "m=400")  
  colnames(mtrees.mean$IE1) <- c("m=100", "m=200", "m=300", "m=400")  
  
  mcse <- mcseTreeCalc(mtrees = mtrees, mtrees.mean = mtrees.mean)

  return(list(
    mean = mtrees.mean,
    mcse = mcse
  ))
}

resCalc <- function(res, X.mat, type = "sd", ce.type = "DE"){
  # Calculate summary statistic using output analysis.
  # Input:
  #   res: matrix with analysis output
  #   X.mat: covariate matrix
  #   type: type of summary statistic ('sd', 'bias', 'rmse')
  #   ce.type: type of causal effect (DE, IE0, or IE1)
  # Output:
  #   matrix containing calculated summary statistics
  
  if (type == "sd"){
    apply(1000 * res, 2, sd)
  } else if (type == "bias"){
    te <- truePM2(
      X = as.matrix(X.mat[,-1]), 
      g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
    )
    if (ce.type == "DE"){
      apply(1000 * res, 2, mean) - te$DE[-1]
    } else if (ce.type == "IE0"){
      apply(1000 * res, 2, mean) - te$IE0[-1]
    } else if (ce.type == "IE1"){
      apply(1000 * res, 2, mean) - te$IE1[-1] 
    }
  } else if (type == "rmse"){
    te <- truePM2(
      X = as.matrix(X.mat[,-1]), 
      g_seq = c(0.25, seq(from = 0.3, to = 0.9, by = 0.1))
    )
    if (ce.type == "DE"){
      apply(t(apply(1000 * res, 1, FUN = function(v){(v - te$DE[-1])^2})), 
            2, FUN = function(x){sqrt(mean(x))})
    } else if (ce.type == "IE0"){
      apply(t(apply(1000 * res, 1, FUN = function(v){(v - te$IE0[-1])^2})), 
            2, FUN = function(x){sqrt(mean(x))})
    } else if (ce.type == "IE1"){
      apply(t(apply(1000 * res, 1, FUN = function(v){(v - te$IE1[-1])^2})), 
            2, FUN = function(x){sqrt(mean(x))}) 
    }
  }
}

mcseTreeCalc <- function(mtrees, mtrees.mean){
  # Calculate the Monte Carlo standard error for a given summary statistic.
  # Input:
  #   mtrees: results from sensitivity study on number of trees
  #   mtrees.mean: mean estimate of summary statistic using MC samples
  # Output:
  #   List with MCSE.
  
  mcse <- list(
    DE = cbind(
      apply(apply(mtrees$DE$m100, 1, FUN = function(est){(est - mtrees.mean$DE[,1])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$DE$m200, 1, FUN = function(est){(est - mtrees.mean$DE[,2])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$DE$m300, 1, FUN = function(est){(est - mtrees.mean$DE[,3])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$DE$m400, 1, FUN = function(est){(est - mtrees.mean$DE[,4])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      )
    ), 
    IE0 = cbind(
      apply(apply(mtrees$IE0$m100, 1, FUN = function(est){(est - mtrees.mean$IE0[,1])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$IE0$m200, 1, FUN = function(est){(est - mtrees.mean$IE0[,2])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$IE0$m300, 1, FUN = function(est){(est - mtrees.mean$IE0[,3])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$IE0$m400, 1, FUN = function(est){(est - mtrees.mean$IE0[,4])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      )
    ), 
    IE1 = cbind(
      apply(apply(mtrees$IE1$m100, 1, FUN = function(est){(est - mtrees.mean$IE1[,1])^2}), 1, 
            FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$IE1$m200, 1, FUN = function(est){(est - mtrees.mean$IE1[,2])^2}), 1, 
            FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$IE1$m300, 1, FUN = function(est){(est - mtrees.mean$IE1[,3])^2}), 1, 
            FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(mtrees$IE1$m400, 1, FUN = function(est){(est - mtrees.mean$IE1[,4])^2}), 1, 
            FUN = function(x){sqrt(sum(x)) / length(x)}
      )
    )
  )
  
  colnames(mcse$DE) <- c("m=100", "m=200", "m=300", "m=400")
  colnames(mcse$IE0) <- c("m=100", "m=200", "m=300", "m=400")
  colnames(mcse$IE1) <- c("m=100", "m=200", "m=300", "m=400")

  return(mcse)  
}

plotNTreesRes <- function(
  plot.name,
  mtree,
  res.type = "sd",
  ce.type = "DE"
){
  # Create a figure showing the estimated summary statistics from the 
  #   sensitivity analysis on number of trees in BART model.
  # Input:
  #   plot.name: figure name
  #   mtree: list with point estimates and monte carlo standard error
  #   res.type: result of interest ('sd', 'bias', 'rmse')
  #   ce.type: causal effect of interest ('DE', 'IE0', 'IE1')
  # Output:
  #   Number of trees plots in supplement (Figures 31 -- 33).
  
  if (res.type == "sd"){
    y.label <- "posterior s.d."
    main.label.p1 <- "Number of trees:  marginal posterior s.d.," 
  } else if (res.type == "bias"){
    y.label <- "bias"
    main.label.p1 <- "Number of trees:  estimated bias,"
  } else if (res.type == "rmse"){
    y.label <- "RMSE"
    main.label.p1 <- "Number of trees:  RMSE,"
  }
  
  if (ce.type == "DE"){
    main.label <- paste(main.label.p1, "DE(g)", sep = " ")
    effect.i <- 1
  } else if (ce.type == "IE0"){
    main.label <- paste(main.label.p1, "IE(0,g)", sep = " ")
    effect.i <- 2
  } else if (ce.type == "IE1"){
    main.label <- paste(main.label.p1, "IE(1,g)", sep = " ")
    effect.i <- 3
  }
  
  png(
    file = plot.name, height = 6, width = 8,
    units = "in", res = 300, bg = "white"
  )
  
  par(
    mfrow = c(2,4), 
    mar = c(5,5,1,1),
    oma = c(0,0,2,0)
  )
  
  # g = 0.3
  k <- 1
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.3",
    col = "blue",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )
  
  # g = 0.4
  k <- 2
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.4",
    col = "darkgreen",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )
  
  # g = 0.5
  k <- 3
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.5",
    col = "red",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )
  
  # iv. g = 0.6
  k <- 4
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.6",
    col = "darkorange",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )
  
  # v. g = 0.7
  k <- 5
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.7",
    col = "purple",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )

  # vi. g = 0.8
  k <- 6
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.8",
    col = "darkred",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )
  
  # vii. g = 0.9
  k <- 7
  plotTreeResults(
    theta.hat = mtree$mean[[effect.i]][k, ],
    theta.mcse = mtree$mcse[[effect.i]][k,],
    main.title = "g = 0.9",
    col = "gold1",
    xlim = c(75, 425),
    ylab = y.label,
    xlab = "number of trees"
  )
  
  title(main = main.label, outer = TRUE, line = 0.25, adj = 0.05, cex.main = 1.25)
  
  dev.off()
}

plotTreeResults <- function(theta.hat, theta.mcse, main.title, col, ...){
  # Individual plot showing point estimates with uncertainty bounds.
  # Input:
  #   theta.hat: point estimates
  #   theta.mcse: monte carlo standard error
  #   main.title: title of plot (e.g., 'g = 0.3')
  #   col: color of plot
  #   ...: other plotting arguments
  # Output:
  #   Plot with point estimates and uncertainty bounds for m = 100, ..., 400.
  
  ub = theta.hat + theta.mcse
  lb = theta.hat - theta.mcse
  
  x_vals <- c(100, 200, 300, 400)
  plot(
    x_vals, theta.hat, pch = 20, bty = "l", 
    ylim = c(min(lb), max(ub)), xaxt = "n", col = col, ...
  )
  arrows(x_vals, theta.hat, x_vals, ub, 
         angle = 90, length = 0.05, col = col)
  arrows(x_vals, theta.hat, x_vals, lb, 
         angle = 90, length = 0.05, col = col)
  
  abline(h = 0, col = "black", lty = 2)
  
  axis(1, at=x_vals, labels=c("100", "200", "300", "400"))
  title(main = main.title, line = -1, adj = 0.05, cex.main = 1)
}


mcseDepthCalc <- function(power, power.mean){
  # Calculate the Monte Carlo standard error for a given summary statistic.
  # Input:
  #   power: results from sensitivity study on tree depth
  #   power.mean: mean estimate of summary statistic using MC samples
  # Output:
  #   List with MCSE.
  
  mcse <- list(
    DE = cbind(
      apply(apply(power$DE$p05, 1, FUN = function(est){(est - power.mean$DE[,1])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p1, 1, FUN = function(est){(est - power.mean$DE[,2])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p15, 1, FUN = function(est){(est - power.mean$DE[,3])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p2, 1, FUN = function(est){(est - power.mean$DE[,4])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p25, 1, FUN = function(est){(est - power.mean$DE[,5])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p3, 1, FUN = function(est){(est - power.mean$DE[,6])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p4, 1, FUN = function(est){(est - power.mean$DE[,7])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$DE$p5, 1, FUN = function(est){(est - power.mean$DE[,8])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      )
    ), 
    IE0 = cbind(
      apply(apply(power$IE0$p05, 1, FUN = function(est){(est - power.mean$IE0[,1])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p1, 1, FUN = function(est){(est - power.mean$IE0[,2])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p15, 1, FUN = function(est){(est - power.mean$IE0[,3])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p2, 1, FUN = function(est){(est - power.mean$IE0[,4])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p25, 1, FUN = function(est){(est - power.mean$IE0[,5])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p3, 1, FUN = function(est){(est - power.mean$IE0[,6])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p4, 1, FUN = function(est){(est - power.mean$IE0[,7])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE0$p5, 1, FUN = function(est){(est - power.mean$IE0[,8])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      )
    ), 
    IE1 = cbind(
      apply(apply(power$IE1$p05, 1, FUN = function(est){(est - power.mean$IE1[,1])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p1, 1, FUN = function(est){(est - power.mean$IE1[,2])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p15, 1, FUN = function(est){(est - power.mean$IE1[,3])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p2, 1, FUN = function(est){(est - power.mean$IE1[,4])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p25, 1, FUN = function(est){(est - power.mean$IE1[,5])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p3, 1, FUN = function(est){(est - power.mean$IE1[,6])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p4, 1, FUN = function(est){(est - power.mean$IE1[,7])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      ),
      apply(apply(power$IE1$p5, 1, FUN = function(est){(est - power.mean$IE1[,8])^2}), 1, 
        FUN = function(x){sqrt(sum(x)) / length(x)}
      )
    )
  )
  
  colnames(mcse$DE) <- c("p=0.5", "p=1", "p=1.5", "p=2", "p=2.5", "p=3", "p=4", "p=5") 
  colnames(mcse$IE0) <- c("p=0.5", "p=1", "p=1.5", "p=2", "p=2.5", "p=3", "p=4", "p=5") 
  colnames(mcse$IE1) <- c("p=0.5", "p=1", "p=1.5", "p=2", "p=2.5", "p=3", "p=4", "p=5") 

  return(mcse)  
}

depthStatCalc <- function(p05, p1, p15, p2, p25, p3, p4, p5, X, res.type = "sd"){
  # Estimate results (marginal standard deviation, bias, and RMSE) from
  #   the sensitivity analysis (tree depth).
  # Input:
  #   p05: results from power = 0.5 analysis
  #   p1: results from power = 1 analysis
  #   p15: results from power = 1.5 analysis
  #   p2: results from power = 2 analysis
  #   p25: results from power = 2.5 analysis
  #   p3: results from power = 3 analysis
  #   p4: results from power = 4 analysis
  #   p5: results from power = 5 analysis
  #   X: covariate matrix
  #   res.type: type of result to calculate ("sd", "bias", "rmse")
  # Output:
  #   List with theta.hat and the associated MCSE.
  
  # number of mcmc chains
  n_chains <- length(p5)
  
  # store results in this data strucure
  power <- list(
    DE = list(
      p05 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p1 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p15 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p2 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p25 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p3 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p4 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p5 = matrix(data = NA_real_, nrow = n_chains, ncol = 7)
    ),
    IE0 = list(
      p05 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p1 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p15 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p2 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p25 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p3 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p4 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p5 = matrix(data = NA_real_, nrow = n_chains, ncol = 7)
    ),
    IE1 = list(
      p05 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p1 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p15 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p2 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p25 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p3 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p4 = matrix(data = NA_real_, nrow = n_chains, ncol = 7),
      p5 = matrix(data = NA_real_, nrow = n_chains, ncol = 7)
    )
  )
  
  for (k in 1:n_chains){
    # DE(g)
    power$DE$p05[k, ] <- resCalc(p05[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p1[k, ] <- resCalc(p1[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p15[k, ] <- resCalc(p15[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p2[k, ] <- resCalc(p2[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p25[k, ] <- resCalc(p25[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p3[k, ] <- resCalc(p3[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p4[k, ] <- resCalc(p4[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")
    power$DE$p5[k, ] <- resCalc(p5[[k]]$DE[,-1], X.mat = X, 
                                   type = res.type, ce.type = "DE")    
    # IE(z=0,g)
    power$IE0$p05[k, ] <- resCalc(p05[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p1[k, ] <- resCalc(p1[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p15[k, ] <- resCalc(p15[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p2[k, ] <- resCalc(p2[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p25[k, ] <- resCalc(p25[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p3[k, ] <- resCalc(p3[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p4[k, ] <- resCalc(p4[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0")
    power$IE0$p5[k, ] <- resCalc(p5[[k]]$IE0[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE0") 
    # IE(z=1,g)
    power$IE1$p05[k, ] <- resCalc(p05[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p1[k, ] <- resCalc(p1[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p15[k, ] <- resCalc(p15[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p2[k, ] <- resCalc(p2[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p25[k, ] <- resCalc(p25[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p3[k, ] <- resCalc(p3[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p4[k, ] <- resCalc(p4[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1")
    power$IE1$p5[k, ] <- resCalc(p5[[k]]$IE1[,-1], X.mat = X, 
                                   type = res.type, ce.type = "IE1") 
  }

  power.mean <- list(
    DE = cbind(
      colMeans(power$DE$p05), 
      colMeans(power$DE$p1),  
      colMeans(power$DE$p15),  
      colMeans(power$DE$p2),
      colMeans(power$DE$p25), 
      colMeans(power$DE$p3),  
      colMeans(power$DE$p4),  
      colMeans(power$DE$p5)
    ),
    IE0 = cbind(
      colMeans(power$IE0$p05), 
      colMeans(power$IE0$p1),  
      colMeans(power$IE0$p15),  
      colMeans(power$IE0$p2),
      colMeans(power$IE0$p25), 
      colMeans(power$IE0$p3),  
      colMeans(power$IE0$p4),  
      colMeans(power$IE0$p5)
    ),
    IE1 = cbind(
      colMeans(power$IE1$p05), 
      colMeans(power$IE1$p1),  
      colMeans(power$IE1$p15),  
      colMeans(power$IE1$p2),
      colMeans(power$IE1$p25), 
      colMeans(power$IE1$p3),  
      colMeans(power$IE1$p4),  
      colMeans(power$IE1$p5)
    )
  )
  
  colnames(power.mean$DE) <- c("p=0.5", "p=1", "p=1.5", "p=2", "p=2.5", "p=3", "p=4", "p=5")  
  colnames(power.mean$IE0) <- c("p=0.5", "p=1", "p=1.5", "p=2", "p=2.5", "p=3", "p=4", "p=5")  
  colnames(power.mean$IE1) <- c("p=0.5", "p=1", "p=1.5", "p=2", "p=2.5", "p=3", "p=4", "p=5")  
  
  mcse <- mcseDepthCalc(power = power, power.mean = power.mean)

  return(list(
    mean = power.mean,
    mcse = mcse
  ))
}

plotDepthRes <- function(
  plot.name,
  power,
  res.type = "sd",
  ce.type = "DE"
){
  # Create a figure showing the estimated summary statistics from the 
  #   sensitivity analysis on tree depth prior in BART model.
  # Input:
  #   plot.name: figure name
  #   power: list with point estimates and monte carlo standard error
  #   res.type: result of interest ('sd', 'bias', 'rmse')
  #   ce.type: causal effect of interest ('DE', 'IE0', 'IE1')
  # Output:
  #   Number of trees plots in supplement (Figures 31 -- 33).
  
  if (res.type == "sd"){
    y.label <- "posterior s.d."
    main.label.p1 <- "Power (controls tree depth):  marginal posterior s.d.," 
  } else if (res.type == "bias"){
    y.label <- "bias"
    main.label.p1 <- "Power (controls tree depth):  estimated bias,"
  } else if (res.type == "rmse"){
    y.label <- "RMSE"
    main.label.p1 <- "Power (controls tree depth):  RMSE,"
  }
  
  if (ce.type == "DE"){
    main.label <- paste(main.label.p1, "DE(g)", sep = " ")
    effect.i <- 1
  } else if (ce.type == "IE0"){
    main.label <- paste(main.label.p1, "IE(0,g)", sep = " ")
    effect.i <- 2
  } else if (ce.type == "IE1"){
    main.label <- paste(main.label.p1, "IE(1,g)", sep = " ")
    effect.i <- 3
  }
  
  png(
    file = plot.name, height = 6, width = 8,
    units = "in", res = 300, bg = "white"
  )
  
  par(
    mfrow = c(2,4), 
    mar = c(5,5,1,1),
    oma = c(0,0,2,0)
  )
  
  # g = 0.3
  k <- 1
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.3",
    col = "blue",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )
  
  # g = 0.4
  k <- 2
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.4",
    col = "darkgreen",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )
  
  # g = 0.5
  k <- 3
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.5",
    col = "red",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )
  
  # iv. g = 0.6
  k <- 4
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.6",
    col = "darkorange",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )
  
  # v. g = 0.7
  k <- 5
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.7",
    col = "purple",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )

  # vi. g = 0.8
  k <- 6
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.8",
    col = "darkred",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )
  
  # vii. g = 0.9
  k <- 7
  plotPowerResults(
    theta.hat = power$mean[[effect.i]][k, ],
    theta.mcse = power$mcse[[effect.i]][k,],
    main.title = "g = 0.9",
    col = "gold1",
    ylab = y.label,
    xlab = "beta (controls tree depth)"
  )
  
  title(main = main.label, outer = TRUE, line = 0.25, adj = 0.05, cex.main = 1.25)
  
  dev.off()
}

plotPowerResults <- function(theta.hat, theta.mcse, main.title, col, ...){
  # Individual plot showing point estimates with uncertainty bounds.
  # Input:
  #   theta.hat: point estimates
  #   theta.mcse: monte carlo standard error
  #   main.title: title of plot (e.g., 'g = 0.3')
  #   col: color of plot
  #   ...: other plotting arguments
  # Output:
  #   Plot with point estimates and uncertainty bounds for m = 100, ..., 400.
  
  ub = theta.hat + theta.mcse
  lb = theta.hat - theta.mcse
  
  x_vals <- c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5)
  plot(
    x_vals, theta.hat, pch = 20, bty = "l", 
    ylim = c(min(lb), max(ub)), xaxt = "n", col = col, ...
  )
  arrows(x_vals, theta.hat, x_vals, ub, 
         angle = 90, length = 0.05, col = col)
  arrows(x_vals, theta.hat, x_vals, lb, 
         angle = 90, length = 0.05, col = col)
  
  abline(h = 0, col = "black", lty = 2)
  
  axis(1, at=x_vals, labels=c("0.5", "1", "1.5", "2", "2.5", "3", "4", "5"))
  title(main = main.title, line = -1, adj = 0.95, cex.main = 1)
}

simSplit <- function(d, alpha, beta){
  # Simulate tree split using (alpha,beta) prior.
  # Input:
  #   d: depth of tree
  #   alpha: hyperparameter 1
  #   beta: hyperparameter 2
  # Output:
  #   New tree (either grown or not)

  prob_split <- alpha * (1 + d)^(-beta)
  split <- rbinom(n = 1, size = 1, prob = prob_split)
  
  if (split){
    # return new branches
    tree_ls <- list(
      simSplit(d + 1, alpha, beta), 
      simSplit(d + 1, alpha, beta)
    )
  } else {
    # return terminal branch
    tree_ls <- list(1)
  }
  
  return(tree_ls)
}

simTermNodes <- function(n_sims, alpha, beta){
  # Simulate the number of terminal nodes for a given set of hyperparameters.
  # Input:
  #   n_sims: number of simulations
  #   alpha: hyperparameter 1
  #   beta: hyperparameter 2
  # Output:
  #   Distribution of number of terminal nodes given (alpha, beta).

  term_nodes <-rep(0, n_sims)
  
  for (k in 1:n_sims){
    term_nodes[k] <- sum(unlist(simSplit(d = 0, alpha, beta)))
  }
  
  return(term_nodes)
}

simMaxDepth <- function(n_sims, alpha, beta){
  # Simulate the max depth of the tree for a given set of hyperparameters.
  # Input:
  #   n_sims: number of simulations
  #   alpha: hyperparameter 1
  #   beta: hyperparameter 2
  # Output:
  #   Distribution of tree depth for given (alpha, beta) hyperparameters.
  max_depth <-rep(0, n_sims)
  
  for (k in 1:n_sims){
    tree_s <- FromListSimple(simSplit(d = 0, alpha, beta))
    max_depth[k] <- tree_s$height
  }

  return(max_depth)
}

plotTermNodes <- function(plot.name, s1, s2, s3, s4, s5, s6, s7, s8){
  # Plot showing the distribution of # of terminal nodes for different
  #   values of (alpha, beta).
  # Input:
  #   plot.name: figure name
  #   s1 - s8: samples from distribution for given (alpha, beta) combo
  # Output:
  #   Figure showing distribution of # of terminal nodes (Figure 34 in Supp)
  
  png(
    file = plot.name, height = 6, width = 8,
    units = "in", res = 300, bg = "white"
  )
  
  par(
    mfrow = c(2,4), 
    mar = c(5,5,1,1),
    oma = c(0,0,2,0)
  )

  mean_prop_loc <- 0.7  

  # i. alpha = 0.95, beta = 0.5
  hist(s1, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=0.5")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim1), 1), sep = ""))

  # ii. alpha = 0.95, beta = 1
  hist(s2, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=1")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim2), 1), sep = ""))

  # iii. alpha = 0.95, beta = 1.5
  hist(s3, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=1.5")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim3), 1), sep = ""))

  # iv. alpha = 0.95, beta = 2
  hist(s4, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=2")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim4), 1), sep = ""))

  # v. alpha = 0.95, beta = 2.5
  hist(s5, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=2.5")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, labels = paste("mean = ", round(mean(tree_sim5), 1), sep = ""))

  # vi. alpha = 0.95, beta = 3
  hist(s6, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=3")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim6), 1), sep = ""))

  # vii. alpha = 0.95, beta = 4
  hist(s7, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=4")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim7), 1), sep = ""))

  # viii. alpha = 0.95, beta = 5
  hist(s8, 
     breaks = seq(from=0, to=75, by=1), 
     xlab = "Number of Terminal Nodes", 
     prob = TRUE,
     ylab = "Probability",
     xlim = c(0,30),
     main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=5")), line = -2, adj = 0.95)
  text(x = 20, y = par('usr')[4] * mean_prop_loc, 
     labels = paste("mean = ", round(mean(tree_sim8), 1), sep = ""))

  title(main = "Number of Terminal Nodes", outer = TRUE, line = 0, adj = 0.05)

  dev.off()
}

plotDepthDist <- function(plot.name, s1, s2, s3, s4, s5, s6, s7, s8){
  # Plot showing the distribution of tree depth for different
  #   values of (alpha, beta).
  # Input:
  #   plot.name: figure name
  #   s1 - s8: samples from distribution for given (alpha, beta) combo
  # Output:
  #   Figure showing tree depth distributions (Figure 35 in Supp)
  
  png(
    file = plot.name, height = 6, width = 8,
    units = "in", res = 300, bg = "white"
  )
  
  par(
    mfrow = c(2,4), 
    mar = c(5,5,1,1),
    oma = c(0,0,2,0)
  )
  
  mean_prop_loc <- 0.7  
  
  # i. alpha = 0.95, beta = 0.5
  hist(s1, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=0.5")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim1), 1), sep = ""))
  
  # ii. alpha = 0.95, beta = 1
  hist(s2, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=1")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim2), 1), sep = ""))
  
  # iii. alpha = 0.95, beta = 1.5
  hist(s3, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=1.5")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim3), 1), sep = ""))
  
  # iv. alpha = 0.95, beta = 2
  hist(s4, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=2")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim4), 1), sep = ""))
  
  # v. alpha = 0.95, beta = 2.5
  hist(s5, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=2.5")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim5), 1), sep = ""))
  
  # vi. alpha = 0.95, beta = 3
  hist(s6, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=3")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim6), 1), sep = ""))
  
  # vii. alpha = 0.95, beta = 4
  hist(s7, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=4")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim7), 1), sep = ""))
  
  # viii. alpha = 0.95, beta = 5
  hist(s8, 
       breaks = seq(from=0, to=75, by=1), 
       xlab = "Max Tree Depth", 
       prob = TRUE,
       ylab = "Probability",
       xlim = c(0,20),
       main = "")
  title(main = expression(paste(alpha, "=0.95, ", beta, "=5")), line = -2, adj = 0.95)
  text(x = 15, y = par('usr')[4] * mean_prop_loc, 
       labels = paste("mean = ", round(mean(depth_sim8), 1), sep = ""))
  
  title(main = "Max Tree Depth", outer = TRUE, line = 0, adj = 0.05)
  
  dev.off()
}
