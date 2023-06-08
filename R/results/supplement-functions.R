### supplement-functions.R
### Nathan Wikle


#########################################################################
### 1. Functions for HyADS comparison
#########################################################################

meanSurface <- function(theta, X_tx, X_mx, mats){
  # Estimate the expected annual sulfate concentration using a given sample
  #   of theta from the sulfate model posterior.
  # Input:
  #   theta: vector of sulfate model parameters
  #   X_tx: emissions vector for US power plants
  #   X_mx: emissions vector for MX power plants
  #   mats: list containing advection-diffusion FVM matrices
  # Output:
  #   Vector representing the estimated expected annual sulfate concentration.

  # grab parameters
  gamma <- theta[1]
  xi <- theta[2]
  beta <- theta[3]
  s2 <- theta[4]
  delta <- theta[5]
  alpha <- theta[6]
  b0 <- theta[7]

  D <- mats$D
  C <- mats$C

  # FVM approximation to advection-diffusion
  n_cells <- nrow(D)
  A.1 <- gamma * D + alpha * C + Matrix::Diagonal(n_cells, delta)
  A.2 <- (gamma / xi) * D + (alpha / xi) * C + Matrix::Diagonal(n_cells, 1)
  # calculate mean (no intercept)
  sources <- solve(A.2, X_tx * beta) + solve(A.2, X_mx * beta)
  mu <- solve(A.1, sources) 
  # return mean
  return(as.vector(mu))
}

breakpoint.creation <- function(ras.obj, n, min.val = NA, max.val = NA) {
  # This function creates a similar color scale for plotting across rasters.
  # Input:
  #   ras.obj: a raster object to be plotted
  #   n: number of color breaks
  #   min.val: minimum value of color scale (default = NA)
  #   max.val: maximum value of color scale (default = NA)
  # Output:
  #   Sequence with breakpoint values to be used when plotting raster.

  if (is.na(min.val)) {
    min.val <- min(values(ras.obj))
  }

  if (is.na(max.val)) {
    max.val <- max(values(ras.obj))
  }

  seq(from = min.val, to = max.val, length.out = n + 1)
}

plotMeanSO4Comp <- function(
  plot.name,
  ou.raster, 
  hyads.raster
){
  # Plot a comparison of the estimated annual average SO4 attributed to coal
  #   power plants using two models: the OU sulfate model and HyADS.
  # Input:
  #   plot.name: name of plot
  #   ou.raster: raster with estimates from OU model
  #   hyads.raster: raster with estimates from HyADS
  # Output:
  #   Supplementary Figure 1.
  
  png(
    file = plot.name,
    width = 12, height = 5, units = "in",
    res = 300, bg = "white" # "transparent",
  )

  par(mfrow = c(1,2), mar = c(3,1,3,8), oma = c(1,1,2,0))


  ### Predicted SO4 surface
  plot(ou.raster,
    legend = FALSE, axes = FALSE,
    breaks = breakpoint.creation(
      ou.raster,
      n = 255,
      min.val = 0,
      max.val = 1.25
    ),
    col = rev(terrain.colors(255)),
    main = "", cex.main = 1.5,
    box = TRUE
  )

  maps::map("state", add = TRUE)

  plot(ou.raster,
    legend.only = TRUE,
    legend.width = 0.75, legend.shrink = 0.55,
    breaks = breakpoint.creation(
      ou.raster,
      n = 255,
      min.val = 0,
      max.val = 1.25 
    ),
    col = rev(terrain.colors(255)),
    axis.args = list(
      at = c(0, 0.25, 0.5, 0.75, 1, 1.25),
      labels = c(0, 0.25, 0.5, 0.75, 1, 1.25),
      cex.axis = 0.75
    ),
    legend.args =
      list(
        text = expression(paste(S * O[4], "  (", mu * g / m^3, ")", sep = "")),
        side = 4, font = 2, line = 2.75, cex = 1
      )
  )
  
  title(
    main = expression(paste("A. ", S * O[4], " from Coal-Fired Power Plants (OU Model)", sep = "")),
    line = 1, adj = 0
  )


  ### HyADS PM2.5 surface
  plot(hyads.raster,
    legend = FALSE, axes = FALSE,
    breaks = breakpoint.creation(
      hyads.raster,
      n = 255,
      min.val = 0,
      max.val = 1.25
    ),
    col = rev(terrain.colors(255)),
    main = "", cex.main = 1.5,
    box = TRUE
  )

  maps::map("state", add = TRUE)

  plot(hyads.raster,
    legend.only = TRUE,
    legend.width = 0.75, legend.shrink = 0.55,
    breaks = breakpoint.creation(
      hyads.raster,
      n = 255,
      min.val = 0,
      max.val = 1.25 
    ),
    col = rev(terrain.colors(255)),
    axis.args = list(
      at = c(0, 0.25, 0.5, 0.75, 1, 1.25),
      labels = c(0, 0.25, 0.5, 0.75, 1, 1.25),
      cex.axis = 0.75
    ),
    legend.args =
      list(
        text = expression(paste(P * M[2.5], "  (", mu * g / m^3, ")", sep = "")),
        side = 4, font = 2, line = 2.75, cex = 1
      )
  )
  
  title(
    main = expression(paste("B. Coal ", P * M[2.5], " (HyADS)", sep = "")),
    line = 1, adj = 0
  )
  
  title(main = "Estimated Pollution Due to Coal-Fired Power Plant Emissions in 2016", 
        line = 0, adj = 0, outer = TRUE)

  dev.off()
}

plotPropSO4Comp <- function(
  plot.name,
  ou.raster, 
  hyads.raster
){
  # Plot a comparison of the estimated proportion of SO4 attributed to the
  #   Big Brown Power Plant using two models: the OU model and HyADS.
  # Input:
  #   plot.name: name of plot
  #   ou.raster: raster with estimates from OU model
  #   hyads.raster: raster with estimates from HyADS
  # Output:
  #   Supplementary Figure 2.
  png(
    file = plot.name,
    width = 12, height = 5, units = "in",
    res = 300, bg = "white" # "transparent",
  )

  par(mfrow = c(1,2), mar = c(3,1,3,8), oma = c(1,1,2,0))

  ### Predicted SO4 surface
  plot(ou.raster,
    legend = FALSE, axes = FALSE,
    breaks = breakpoint.creation(
      ou.raster,
      n = 255,
      min.val = 0,
      max.val = 0.45
    ),
    col = rev(heat.colors(255)),
    main = "", cex.main = 1.5,
    box = TRUE
  )

  maps::map("state", add = TRUE)

  plot(ou.raster,
    legend.only = TRUE,
    legend.width = 0.75, legend.shrink = 0.55,
    breaks = breakpoint.creation(
      ou.raster,
      n = 255,
      min.val = 0,
      max.val = 0.45 
    ),
    col = rev(heat.colors(255)),
    axis.args = list(
      at = c(0, 0.1, 0.2, 0.3, 0.4),
      labels = c(0, 0.1, 0.2, 0.3, 0.4),
      cex.axis = 0.75
    ),
    legend.args =
      list(
        text = "Proportion",
        side = 4, font = 2, line = 2.75, cex = 1
      )
  )
  
  title(
    main = expression(paste("A. Proportion of ", S * O[4], " Due to Fac. 3497 (OU Model)", sep = "")),
    line = 1, adj = 0
  )


  ### HyADS PM2.5 surface
  plot(hyads.raster,
    legend = FALSE, axes = FALSE,
    breaks = breakpoint.creation(
      hyads.raster,
      n = 255,
      min.val = 0,
      max.val = 0.45
    ),
    col = rev(heat.colors(255)),
    main = "", cex.main = 1.5,
    box = TRUE
  )

  maps::map("state", add = TRUE)

  plot(hyads.raster,
    legend.only = TRUE,
    legend.width = 0.75, legend.shrink = 0.55,
    breaks = breakpoint.creation(
      hyads.raster,
      n = 255,
      min.val = 0,
      max.val = 0.45 
    ),
    col = rev(heat.colors(255)),
    axis.args = list(
      at = c(0, 0.1, 0.2, 0.3, 0.4),
      labels = c(0, 0.1, 0.2, 0.3, 0.4),
      cex.axis = 0.75
    ),
    legend.args =
      list(
        text = "Proportion",
        side = 4, font = 2, line = 2.75, cex = 1
      )
  )
  
  title(
    main = expression(paste("B. Proportion of Coal ", P * M[2.5], " Due to Fac. 3497 (HyADS)", sep = "")),
    line = 1, adj = 0
  )
  
  title(main = "Estimated Proportion of Pollution Due to Emissions from Facility 3497", 
      line = 0, adj = 0, outer = TRUE)

  dev.off()
}


#########################################################################
### 2. Functions for comparison of covariate balance and overlap 
#########################################################################

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

upwindCovar <- function(theta, mats, X_mat, P_mat, s_vec, key, covars, version = 2){
  # Calculate the upwind parameters for a given sample of theta 
  #   (sulfate model parameters).
  # Input:
  #   theta:   Vector of SO4 model parameters.
  #   mats:    List containing advection and diffusion FVM matrices: C and D.
  #   X_mat:   Matrix of counterfactual emissions scenarios. 
  #   P_mat:   Projection matrix mapping raster elements to ZCTAs.
  #   s_vec:   Vector listing treatment status of key-associated power plants.
  #   key:     Vector of key-associated power plants.
  #   covars:  Matrix containing covariates used to in upwind calculation.
  #   version: The version of upwind treatments to calculate (default = 2)
  # Output:
  #   A list with the calculated upwind covariate values for each ZCTA.
  
  # sulfate model parameters
  gamma <- theta[1]
  xi <- theta[2]
  beta <- theta[3]
  s2 <- theta[4]
  delta <- theta[5]
  alpha <- theta[6]
  b0 <- theta[7]

  # advection-diffusion FVM matrices
  D <- mats$D
  C <- mats$C

  # calculate mu
  n_cells <- nrow(D)
  A.1 <- gamma * D + alpha * C + Diagonal(n_cells, delta)
  A.2 <- (gamma / xi) * D + (alpha / xi) * C + Diagonal(n_cells, 1)
  sources <- solve(A.2, X_mat * beta)
  mu <- solve(A.1, sources)

  # determine source-receptor matrix (T)
  T_mat <- Matrix::crossprod(P_mat, mu)

  # weighted degree
  T_sum <- rowSums(T_mat) - T_mat[cbind(1:nrow(T_mat), key)]

  # save results in a list
  res_list <- list()
  res_list[[1]] <- T_sum

  if (version == 2){
    # determine average upwind treatment effect
    T_prop_mat <- T_mat / T_sum
    T_prop_mat[cbind(1:nrow(T_mat), key)] <- 0
    heat_upwind <- T_prop_mat %*% as.matrix(covars)
    res_list[[2]] <- log(as.vector(heat_upwind))
  }

  return(res_list)
}

calc.balance <- function(cov.matrix, trt1, trt2, weights, d.type = 2){
  # Calculate the absolute mean difference in treatment.
  # Input:
  #   cov.matrix: covariate matrix
  #   trt1: vector indicating indices of treated units
  #   trt2: vector indicating indices of untreated units
  #   weights: propensity score weights
  #   d.type: type of standard deviation used (default = 2: pooled standard deviation)
  # Output:
  #   Matrix with estimated absolute mean differences for each covariate.
  diff <- balance.est(
    m1 = cov.matrix[trt1, ],
    wt1 = weights[trt1],
    m2 = cov.matrix[trt2, ],
    wt2 = weights[trt2],
    d.type = 2
  )
  
  # return mean differences
  return(diff)
}

balance.est <- function(m1, m2, wt1, wt2, d.type){
  # Estimate covariate balance (absolute mean difference).
  # Input:
  #   m1: covariate matrix for treated units
  #   wt1: propensity score weights for treated units
  #   m2: covariate matrix for untreated units
  #   wt2: propensity score weights for untreated units
  #   d.type: standard deviation used 
  #     (1 = sd of treated group, 2 = pooled standard deviation)
  # Output: 
  #   Matrix with estimated absolute mean difference between groups
  
  if (d.type == 1){
    # standard deviation of the treated group
    denom <- apply(m1, 2, sd)
  } else {
    # pooled sd (square root of the mean of the group variances)
    trt.var <- apply(m1, 2, var)
    notrt.var <- apply(m2, 2, var)
    
    denom <- sqrt(rowMeans(cbind(trt.var, notrt.var)))
  }
  
  diff.mat <- matrix(NA_real_, nrow = ncol(m1), ncol = 2)
  colnames(diff.mat) <- c("Diff.Un", "Diff.Adj")
  rownames(diff.mat) <- colnames(m1)
  
  diff.mat[,1] <- abs(colMeans(m1) - colMeans(m2)) / denom
  
  m1.wt <- apply(m1, 2, function(x) {
    wt1 * x / sum(wt1)
  })
  
  m2.wt <- apply(m2, 2, function(x) {
    wt2 * x / sum(wt2)
  })
  
  diff.mat[,2] <- abs(colSums(m1.wt) - colSums(m2.wt)) / denom
  
  return(diff.mat)
}

gCBPS.est <- function(g.new, cbps.output){
  # Returns the generalized propensity score from fitted CBPS output.
  # Input:
  #   g.new: 
  #   cbps.output: output from call to CBPS
  # Output:
  #   Generalized CBPS.
  
  # standardize g
  g.bar <- attr(cbps.output$Ttilde, "scaled:center")
  g.sd <- attr(cbps.output$Ttilde, "scaled:scale")
  g.star <- (g.new - g.bar) / g.sd
  
  # calculate propensity score
  Xbeta.tilde <- cbps.output$Xtilde %*% cbps.output$beta.tilde
  s2.tilde <- cbps.output$sigmasq.tilde
  exp(-(g.star - Xbeta.tilde)^2 / 2 / s2.tilde) / sqrt(2 * pi * s2.tilde)
}

plotZBalance <- function(plot.name, diff.mat){
  # Creates plot showing covariate balance between Z groups.
  # Input:
  #   plot.name: name of figure
  #   diff.mat: matrix containing absolute standardized mean differences
  # Output:
  #   Figure 3 from supplement.
  
  png(
    file = plot.name, 
    height = 10, width = 8, units = "in",
    res = 300, bg = "white"
  )
  
  par(
    mfrow = c(1, 1),
    mar = c(5, 7, 5, 4)
  )
  
  plot(
    abs(diff.mat[, 1]), 1:nrow(diff.mat),
    bty = "l",
    xlab = "Absolute Standardized Mean Differences",
    xlim = c(0, max(abs(diff.mat))),
    ylab = "",
    pch = 15, col = "black",
    main = "Balance Assessment on Z",
    yaxt = "n",
  )
  
  points(
    abs(diff.mat[, 2]), 1:nrow(diff.mat),
    pch = 16, col = "red"
  )
  
  points(
    abs(diff.mat[, 3]), 1:nrow(diff.mat),
    pch = 17, col = "blue"
  )
  
  axis(2,
       at = 1:nrow(diff.mat),
       labels = rownames(diff.mat), las = 2
  )
  
  abline(h = 1:nrow(diff.mat), col = "grey", lty = "dotted")
  
  legend(0.4, 20, bg = "white",
         legend = c("Weighted (CBPS)", "Weighted (BART)", "Unweighted"),
         col = c("blue", "red", "black"), pch = c(17, 16, 15),
         cex = 1, box.lty = 1, y.intersp = 2
  )
  
  abline(v = 0.05, col = "red")
  
  dev.off()
}

plotZOverlap <- function(plot.name, z, prob.z){
  # Plot the propensity score overlap for Z (Figure s4).
  # Input:
  #   plot.name: figure name
  #   z: vector of z values
  #   prob.z: estimated probability of treatment
  # Output:
  #   Figure 4 from supplement
  
  png(
    file = plot.name, 
    height = 6, width = 7, units = "in",
    res = 300, bg = "white"
  )

  par(mfrow = c(1,1), mar = c(4,4,2,2), oma = c(2,2,2,2))

  # histograms
  hist(prob.z[z == 1], col = alpha("red", 0.4), 
       main = "", xlab = "propensity score")
  hist(prob.z[z == 0], col = alpha("yellow", 0.4), add = TRUE)
  title(main = "Propensity Overlap: Direct Treatment", adj = 0.05, line = 1)

  # add labels
  text(x = 0.18, y = 175, "Z = 0", font = 2)
  # text(x = 0.18, y = 170, "Z = 0")
  text(x = 0.7, y = 175, "Z = 1", font = 2)

  dev.off()
}

plotGBalance <- function(plot.name, corr.mat){
  # Creates plot showing covariate balance on G.
  # Input:
  #   plot.name: name of figure
  #   corr.mat: matrix containing propensity-weighted correlations
  # Output:
  #   Figure 5 from supplement.
  
  png(
    file = plot.name, 
    height = 10, width = 8, units = "in",
    res = 300, bg = "white"
  )
  
  par(
    mfrow = c(1, 1),
    mar = c(5, 7, 5, 4)
  )
  
  plot(
    abs(corr.mat[, 1]), 1:nrow(corr.mat),
    bty = "l",
    xlab = "Correlation",
    xlim = c(0, max(abs(corr.mat))),
    ylab = "",
    pch = 15, col = "black",
    main = "Balance Assessment on G",
    yaxt = "n",
  )
  
  points(
    abs(corr.mat[, 2]), 1:nrow(corr.mat),
    pch = 16, col = "red"
  )
  
  points(
    abs(corr.mat[, 3]), 1:nrow(corr.mat),
    pch = 17, col = "blue"
  )
  
  axis(2,
       at = 1:nrow(corr.mat),
       labels = rownames(corr.mat), las = 2
  )
  
  abline(h = 1:nrow(corr.mat), col = "grey", lty = "dotted")
  
  legend(0.35, 10, bg = "white",
         legend = c("Weighted (CBPS)", "Weighted (BART)", "Unweighted"),
         col = c("blue", "red", "black"), pch = c(17, 16, 15),
         cex = 1, box.lty = 1, y.intersp = 2
  )
  
  abline(v = 0.1, col = "red")
  
  dev.off()
}

plotGOverlap <- function(
  plot.name, g.cutoff, g.propensity, g.quantiles, cols, trans
){
  # Plot the generalized propensity score overlap for G (Figure s6).
  # Input:
  #   plot.name: figure name
  #   g.cutoff: quantile values of g
  #   g.propensity: fitted CBPS model for g
  #   g.quantiles: list with indices dependent on quantile of each unit
  #   cols: colors used in plot
  #   trans: transparency of color
  # Output:
  #   Figure 6 from supplement
  
  png(
    file = plot.name, height = 8, width = 9,
    units = "in", res = 300, bg = "white"
  )
  
  par(mfrow = c(3,4),
      mar = c(2,4,3,2), 
      oma = c(3,1,1,1))
  
  # 1st quartile, f(g = 1st quartile | x)
  rhat.low <- gCBPS.est(g.new = g.cutoff[1], cbps.output = g.propensity) 

  hist(rhat.low[g.quantiles[[1]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[1], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "First 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)
  
  hist(rhat.low[g.quantiles[[2]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[2], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Second 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)
  
  hist(rhat.low[g.quantiles[[3]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[3], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Third 25%", adj = 0.1, line = 0, cex.main = 1, font.main = 1)
  
  hist(rhat.low[g.quantiles[[4]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[4], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Fourth 25%", adj = 0.1, line = 0, cex.main = 1, font.main = 1)
  
  mtext(bquote(paste('A. Overlap of f(g = Q'[1]*' | x)')), line = -2.3, 
        outer = TRUE, adj = 0)

  # 2nd quartile, f(g = median | x)
  rhat.mid <- gCBPS.est(g.new = g.cutoff[2], cbps.output = g.propensity) 

  # par(mfrow = c(2,2))
  hist(rhat.mid[g.quantiles[[1]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[1], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "First 25%", adj = 0.3, line = -3, cex.main = 1, font.main = 1)

  hist(rhat.mid[g.quantiles[[2]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[2], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Second 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)

  hist(rhat.mid[g.quantiles[[3]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[3], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Third 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)

  hist(rhat.mid[g.quantiles[[4]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[4], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Fourth 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)

  mtext(bquote(paste('B. Overlap of f(g = Q'[2]*' | x)')), line = -21.3, 
        outer = TRUE, adj = 0)
  

  ## C. 3rd quartile, f(g = 3rd quartile | x)
  rhat.high <- gCBPS.est(g.new = g.cutoff[3], cbps.output = g.propensity)  

  # par(mfrow = c(2,2))
  hist(rhat.high[g.quantiles[[1]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[1], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "First 25%", adj = 0.3, line = -3, cex.main = 1, font.main = 1)
  
  hist(rhat.high[g.quantiles[[2]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[2], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Second 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)
  
  hist(rhat.high[g.quantiles[[3]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[3], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Third 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)
  
  hist(rhat.high[g.quantiles[[4]]], main = "", xlab = "", xaxt='n', 
       col = alpha(cols[4], trans))
  axis(side = 1, at = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1))
  title(main = "Fourth 25%", adj = 0.1, line = -3, cex.main = 1, font.main = 1)
  
  mtext(bquote(paste('C. Overlap of f(g = Q'[3]*' | x)')), font = 2, 
        line = -40.4, outer = TRUE, adj = 0)
  
  
  mtext("Generalized Propensity Score: f(g | x)", side = 1, line = 1, 
        cex = 0.75, outer = TRUE)
  
  dev.off()
}


#########################################################################
### 3. Functions for plotting causal estimates stratified by covariates
#########################################################################

plotEffect <- function(
    res.nounc, res.unc, y.lims = NA, 
    c1 = "darkblue", c2 = "darkred", trans = 0.5,
    type = 1
){
  # Plots the estimated causal effect curves using given output from 
  #   models fit with and without uncertainty propagation.
  # Input:
  #   res.nounc: estimated effects from the model without uncertainty 
  #     propagation
  #   res.unc: estimated effects from teh model with uncertainty propagation
  #   y.lims: y-axis limits for the plot; NA = estimated from 'res' data
  #   c1: color for estimates with uncertainty propagation
  #   c2: color for estimates without uncertainty propagation
  #   trans: color transparency
  #   type: outcome used in the analysis (1 = medicare, 2 = asthma)
  # Output:
  #   Figure showing estimated causal effect curves with and without 
  #     uncertainty propagation.
  
  if (type == 1){
    ylab.title = "mortality rate per 1000"
  } else if (type == 2){
    ylab.title = "asthma ED rate per 1000"
  }
  
  # two step
  mu.nounc.hat <- 1000 * colMeans(res.nounc)
  mu.nounc.lb <- 1000 * apply(res.nounc, 2, quantile, probs = 0.025, na.rm = TRUE)
  mu.nounc.ub <- 1000 * apply(res.nounc, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  # cut feedback
  mu.unc.hat <- 1000 * colMeans(res.unc, na.rm = TRUE)
  mu.unc.lb <- 1000 * apply(res.unc, 2, quantile, probs = 0.025, na.rm = TRUE)
  mu.unc.ub <- 1000 * apply(res.unc, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  g.seq <- seq(from = 0.25, to = 0.9, by = 0.01)
  
  if (is.na(y.lims[1])){
    y.lims <- c(min(mu.nounc.lb, mu.unc.lb), max(mu.nounc.ub, mu.unc.ub))
  } 
  
  # plot with boundaries
  plot(g.seq, mu.nounc.hat,
       type = "l", ylab = ylab.title, xlab = "g",
       bty = "l", col = "gray90", ylim = y.lims,
       cex.axis = 1.25, cex.lab = 1.5
  )
  
  # cut feedback (wider intervals...)
  polygon(
    c(g.seq, rev(g.seq)), c(mu.unc.lb, rev(mu.unc.ub)),
    col = alpha(c1, trans), lty = 0
  )
  
  # plugin (no uncertainty propagation)
  polygon(
    c(g.seq, rev(g.seq)), c(mu.nounc.lb, rev(mu.nounc.ub)),
    col = alpha(c2, trans), lty = 0
  )
  
  # plug-in lines
  lines(g.seq, mu.nounc.hat, col = c2, lwd = 3)

  # cut feedback lines
  lines(g.seq, mu.unc.hat, col = c1, lty = 2, lwd = 3)
}

plotStratEffects <- function(
  plot.name,
  bart.plugin,
  bart.cut,
  y.lims,
  cols,
  outcome.type = 1,
  ce.type = "DE"
){
  # Creates a plot showing the estimated causal effects stratified by 
  #   degree and distance from key-assocaited facility.
  # Input:
  #   plot.name: figure name
  #   bart.plugin: results from BART model with plugin inference
  #   bart.cut: results from BART model with uncertainty propagation
  #   y.lims: y-axis limits
  #   cols: vector of colors
  #   outcome.type: outcome used in the analysis (1 = medicare, 2 = asthma)
  #   ce.type: type of causal effect to plot (DE, IE0, or IE1)
  # Output:
  #   Figure 7 or 8 from supplement
  
  if (ce.type == "DE"){
    res.k <- 3
  } else if (ce.type == "IE0"){
    res.k <- 4
  } else if (ce.type == "IE1"){
    res.k <- 5
  }

  png(
    file = plot.name, height = 8, width = 20,
    units = "in", res = 300, bg = "white"
  )
  
  par(mfrow = c(2,4), mar = c(6,4.3,3,2), oma = c(1,6,6,1))
  
  # Degree
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$ldeg,
    res.unc = bart.cut[[res.k]]$ldeg,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "Low", line = -1, adj = 0.1, cex.main = 1.5)
  
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$mdeg,
    res.unc = bart.cut[[res.k]]$mdeg,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "Medium", line = -1, adj = 0.1, cex.main = 1.5)
  
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$hdeg,
    res.unc = bart.cut[[res.k]]$hdeg,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "High", line = -1, adj = 0.1, cex.main = 1.5)
  
  
  plot.new()
  legend(x = 0.1, y = 0.8, legend = c("plug-in", "cut feedback"), col = cols[4:3], 
         lty = c(1, 2), lwd = 4, bty = "n", cex = 2, y.intersp = 1.5)
  
  
  mtext("Degree", line = -2, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  # Distance
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$d1,
    res.unc = bart.cut[[res.k]]$d1,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "[0, 50) km", line = -1, adj = 0.1, cex.main = 1.5)
  
  
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$d2,
    res.unc = bart.cut[[res.k]]$d2,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "[50, 100) km", line = -1, adj = 0.1, cex.main = 1.5)
  
  
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$d3,
    res.unc = bart.cut[[res.k]]$d3,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "[100, 200) km", line = -1, adj = 0.1, cex.main = 1.5)
  
  plotEffect(
    res.nounc = bart.plugin[[res.k]]$d4,
    res.unc = bart.cut[[res.k]]$d4,
    y.lims = y.lims, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = outcome.type
  )
  title(main = "200+ km", line = -1, adj = 0.1, cex.main = 1.5)
  
  mtext("Distance to Key-Associated Facility", line = -28.5, adj = 0, 
        outer = TRUE, font = 2, cex = 1.5)
  
  if (ce.type == "DE"){
     mtext("A.  Direct Effects:  DE(g)", line = 1, adj = -0.02, 
        outer = TRUE, font = 2, cex = 1.7)
  } else if (ce.type == "IE0"){
    mtext("B.  Indirect Effects:  IE(0,g)", line = 1, adj = -0.02, 
          outer = TRUE, font = 2, cex = 1.7)
  } else if (ce.type == "IE1"){
    mtext("C.  Indirect Effects:  IE(1,g)", line = 1, adj = -0.02, 
          outer = TRUE, font = 2, cex = 1.7)
  }
  
  dev.off()
}


#########################################################################
### 4. Functions for Major Findings from Simulation Study 
#########################################################################

plotPM1Results <- function(plot.name, p.plugin, p.cut, b.plugin, b.cut, p.var,b.var){
  # Plots the coverage rate and proportion of variance results from simulation
  #   study PM1.
  # Input:
  #   plot.name: figure name
  #   p.plugin: coverage results for Poisson plugin model
  #   p.cut: coverage results for Poisson cut model
  #   b.plugin: coverage results for log-linear BART plugin model
  #   b.cut: coverage results for log-linear BART cut model
  #   p.var: variance results for Poisson cut model
  #   b.var: variance results for log-linear BART cut model
  # Output:
  #   Figure 10 from the supplement.

  png(
    file = plot.name, height = 6, width = 12,
    units = "in", res = 300, bg = "white"
  )
  
  par(
    mfrow = c(2,4),
    mar = c(4,4,4,2),
    oma = c(1,1,1,1)
  )

  g.vals <- c(0.3,0.4,0.5,0.6,0.7,0.8)

  # DE(g)
  k = 1
  plot(
    g.vals, p.plugin[,k], type = "l", lty = 2, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "coverage rate", xlab = "g"
  )

  lines(g.vals, p.cut[,k], lty = 1, lwd = 2, col = "darkred")
  lines(g.vals, b.cut[,k], col = "darkblue", lwd = 2)
  title(main = "DE(g)", adj = 0, line = 0.5, cex.main = 1)

  # IE(0,g)
  k = 2
  plot(
    g.vals, p.plugin[,k], type = "l", lty = 2, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "coverage rate", xlab = "g"
  )
  lines(g.vals, p.cut[,k], lty = 1, lwd = 2, col = "darkred")
  lines(g.vals, b.cut[,k], col = "darkblue", lwd = 2)
  title(main = "IE(0,g)", adj = 0, line = 0.5, cex.main = 1)

  # IE(1,g)
  k = 3
  plot(
    g.vals, p.plugin[,k], type = "l", lty = 2, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "coverage rate", xlab = "g"
  )
  lines(g.vals, p.cut[,k], lty = 1, lwd = 2, col = "darkred")
  lines(g.vals, b.cut[,k], col = "darkblue", lwd = 2)
  title(main = "IE(1,g)", adj = 0, line = 0.5, cex.main = 1)
  
  # legend
  par(mar = c(0,0,1,1) + 0.1)
  plot.new()
  legend(x = 0, y = 0.8,
    legend = c("Poisson plug-in", "Poisson cut", "log-linear BART cut*"),
    bty = "n", col = c("darkred", "darkred", "darkblue"),
    lty = c(2,1,1),
    lwd = c(2,2,2),
    cex = 1.25
  )

  text(x = 0.1, y = 0.4,  " * Coverage rates for log-linear BART ", adj = 0)
  text(x = 0.1, y = 0.34, "     were at or near 1 for both plug-in", adj = 0)
  text(x = 0.1, y = 0.28,  "     and cut methods.", adj = 0)
  
  mtext("A.  95% Credible Interval Coverage Rates", font = 1, line = -2,
        outer = TRUE, adj = 0)

  ### proportion of variance
  par(mar = c(4,4,4,2))
  g.vals <- c(0.3,0.4,0.5,0.6,0.7,0.8)
  
  # DE(g)
  k = 1
  plot(
    g.vals, p.var[,k], type = "l", lty = 3, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "proportion of variance", xlab = "g"
  )
  lines(g.vals, b.var[,k], lty = 1, lwd = 2, col = "darkblue")
  title(main = "DE(g)", adj = 0, line = 0.5, cex.main = 1)
  
  # IE(0,g)
  k = 2
  plot(
    g.vals, p.var[,k], type = "l", lty = 3, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "proportion of variance", xlab = "g"
  )
  lines(g.vals, b.var[,k], lty = 1, lwd = 2, col = "darkblue")
  title(main = "IE(0,g)", adj = 0, line = 0.5, cex.main = 1)
  
  # IE(1,g)
  k = 3
  plot(
    g.vals, p.var[,k], type = "l", lty = 3, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "proportion of variance", xlab = "g"
  )
  lines(g.vals, b.var[,k], lty = 1, lwd = 2, col = "darkblue")
  title(main = "IE(1,g)", adj = 0, line = 0.5, cex.main = 1)
  
  par(mar = c(0,0,1,1) + 0.1)
  plot.new()
  legend(x = 0, y = 0.8, 
         legend = c("Poisson", "log-linear BART"),
         bty = "n", col = c("darkred", "darkblue"), lty = c(3,1), 
         lwd = c(2,2), cex = 1.25
  )
  
  mtext("B.  Proportion of Variance Attributed to Interference", 
        font = 1, line = -23.7, outer = TRUE, adj = 0)
  
  dev.off()
}

plotPM3Results <- function(plot.name, p.plugin, p.cut, b.plugin, b.cut, p.var, b.var){
  # Plots the coverage rate and proportion of variance results from simulation
  #   study PM3.
  # Input:
  #   plot.name: figure name
  #   p.plugin: coverage results for Poisson plugin model
  #   p.cut: coverage results for Poisson cut model
  #   b.plugin: coverage results for log-linear BART plugin model
  #   b.cut: coverage results for log-linear BART cut model
  #   p.var: variance results for Poisson cut model
  #   b.var: variance results for log-linear BART cut model
  # Output:
  #   Figure 11 from the supplement.
  
  png(
    file = plot.name, height = 6, width = 12,
    units = "in", res = 300, bg = "white"
  )
  
  par(
    mfrow = c(2,4),
    mar = c(4,4,4,2),
    oma = c(1,1,1,1)
  )
  
  g.vals <- c(0.3,0.4,0.5,0.6,0.7,0.8)
  
  # DE(g)
  k = 1
  plot(
    g.vals, p.plugin[,k], type = "l", lty = 2, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "coverage rate", xlab = "g"
  )
  lines(g.vals, p.cut[,k], lty = 1, lwd = 2, col = "darkred")
  lines(g.vals, b.plugin[,k], lty = 2, col = "darkblue", lwd = 2)
  lines(g.vals, b.cut[,k], col = "darkblue", lwd = 2)
  title(main = "DE(g)", adj = 0, line = 0.5, cex.main = 1)
  
  # IE(0,g)
  k = 2
  plot(
    g.vals, p.plugin[,k], type = "l", lty = 2, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "coverage rate", xlab = "g"
  )
  lines(g.vals, p.cut[,k], lty = 1, lwd = 2, col = "darkred")
  lines(g.vals, b.plugin[,k], lty = 2, col = "darkblue", lwd = 2)
  lines(g.vals, b.cut[,k], col = "darkblue", lwd = 2)
  title(main = "IE(0,g)", adj = 0, line = 0.5, cex.main = 1)
  
  # IE(1,g)
  k = 3
  plot(
    g.vals, p.plugin[,k], type = "l", lty = 2, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "coverage rate", xlab = "g"
  )
  lines(g.vals, p.cut[,k], lty = 1, lwd = 2, col = "darkred")
  lines(g.vals, b.plugin[,k], lty = 2, col = "darkblue", lwd = 2)
  lines(g.vals, b.cut[,k], col = "darkblue", lwd = 2)
  title(main = "IE(1,g)", adj = 0, line = 0.5, cex.main = 1)
  
  # legend
  par(mar = c(0,0,1,1) + 0.1)
  plot.new()
  legend(x = 0, y = 0.8, 
    legend = c(
      "Poisson plug-in", "Poisson cut", 
      "log-linear BART plug-in", "log-linear BART cut"
    ),
    bty = "n", col = c("darkred", "darkred", "darkblue", "darkblue"), 
    lty = c(2,1,2,1), 
    lwd = c(2,2,2,2), 
    cex = 1.25
  )
  
  mtext("A.  95% Credible Interval Coverage Rates", font = 1, line = -2, 
        outer = TRUE, adj = 0)
  
  
  ### proportion of variance
  
  par(mar = c(4,4,4,2))
  
  g.vals <- c(0.3,0.4,0.5,0.6,0.7,0.8)
  
  # DE(g)
  k = 1
  plot(
    g.vals, p.var[,k], type = "l", lty = 3, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "proportion of variance", xlab = "g"
  )
  lines(g.vals, b.var[,k], lty = 1, lwd = 2, col = "darkblue")
  title(main = "DE(g)", adj = 0, line = 0.5, cex.main = 1)
  
  # IE(0,g)
  k = 2
  plot(
    g.vals, p.var[,k], type = "l", lty = 3, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "proportion of variance", xlab = "g"
  )
  lines(g.vals, b.var[,k], lty = 1, lwd = 2, col = "darkblue")
  title(main = "IE(0,g)", adj = 0, line = 0.5, cex.main = 1)
  
  # IE(1,g)
  k = 3
  plot(
    g.vals, p.var[,k], type = "l", lty = 3, lwd = 2, bty = "l", col = "darkred",
    ylim = c(0,1), ylab = "proportion of variance", xlab = "g"
  )
  lines(g.vals, b.var[,k], lty = 1, lwd = 2, col = "darkblue")
  title(main = "IE(1,g)", adj = 0, line = 0.5, cex.main = 1)
  
  # legend
  par(mar = c(0,0,1,1) + 0.1)
  plot.new()
  legend(x = 0, y = 0.8, 
    legend = c("Poisson", "log-linear BART"),
    bty = "n", col = c("darkred", "darkblue"), lty = c(3,1), 
    lwd = c(2,2), cex = 1.25
  )
  
  mtext("B.  Proportion of Variance Attributed to Interference", 
        font = 1, line = -23.7, outer = TRUE, adj = 0)
  
  dev.off()
}

deltaCoverage <- function(results){
  # Combines coverage and variance statistics from each simulation study
  #   into a single data frame.
  # Input:
  #   results: results from the simulation studies
  # Output:
  #   Data frame containing coverage and variance pairs.
  
  # coverage results
  coverage <- results$coverage
  # variance results
  variance <- results$variance
  
  # change in coverage
  delta <- c(
    c(coverage[[1]][[1]]$cut[-1,] - coverage[[1]][[1]]$plugin[-1,]),
    c(coverage[[1]][[2]]$cut[-1,] - coverage[[1]][[2]]$plugin[-1,]),
    c(coverage[[2]][[1]]$cut[-1,] - coverage[[2]][[1]]$plugin[-1,]),
    c(coverage[[2]][[2]]$cut[-1,] - coverage[[2]][[2]]$plugin[-1,]),
    c(coverage[[3]][[1]]$cut[-1,] - coverage[[3]][[1]]$plugin[-1,]),
    c(coverage[[3]][[2]]$cut[-1,] - coverage[[3]][[2]]$plugin[-1,]),
    c(coverage[[4]][[1]]$cut[-1,] - coverage[[4]][[1]]$plugin[-1,]),
    c(coverage[[4]][[2]]$cut[-1,] - coverage[[4]][[2]]$plugin[-1,]),
    c(coverage[[5]][[1]]$cut[-1,] - coverage[[5]][[1]]$plugin[-1,]),
    c(coverage[[5]][[2]]$cut[-1,] - coverage[[5]][[2]]$plugin[-1,]),
    c(coverage[[6]][[1]]$cut[-1,] - coverage[[6]][[1]]$plugin[-1,]),
    c(coverage[[6]][[2]]$cut[-1,] - coverage[[6]][[2]]$plugin[-1,])
  )
  
  # % variance due to interference
  inter.var <- 100 * c(
    c(variance[[1]][[1]]$DE[-1,2] / variance[[1]][[1]]$DE[-1,3]),
    c(variance[[1]][[1]]$IE0[-1,2] / variance[[1]][[1]]$IE0[-1,3]),
    c(variance[[1]][[1]]$IE1[-1,2] / variance[[1]][[1]]$IE1[-1,3]),
    c(variance[[1]][[2]]$DE[-1,2] / variance[[1]][[2]]$DE[-1,3]),
    c(variance[[1]][[2]]$IE0[-1,2] / variance[[1]][[2]]$IE0[-1,3]),
    c(variance[[1]][[2]]$IE1[-1,2] / variance[[1]][[2]]$IE1[-1,3]),
    c(variance[[2]][[1]]$DE[-1,2] / variance[[2]][[1]]$DE[-1,3]),
    c(variance[[2]][[1]]$IE0[-1,2] / variance[[2]][[1]]$IE0[-1,3]),
    c(variance[[2]][[1]]$IE1[-1,2] / variance[[2]][[1]]$IE1[-1,3]),
    c(variance[[2]][[2]]$DE[-1,2] / variance[[2]][[2]]$DE[-1,3]),
    c(variance[[2]][[2]]$IE0[-1,2] / variance[[2]][[2]]$IE0[-1,3]),
    c(variance[[2]][[2]]$IE1[-1,2] / variance[[2]][[2]]$IE1[-1,3]),    
    c(variance[[3]][[1]]$DE[-1,2] / variance[[3]][[1]]$DE[-1,3]),
    c(variance[[3]][[1]]$IE0[-1,2] / variance[[3]][[1]]$IE0[-1,3]),
    c(variance[[3]][[1]]$IE1[-1,2] / variance[[3]][[1]]$IE1[-1,3]),
    c(variance[[3]][[2]]$DE[-1,2] / variance[[3]][[2]]$DE[-1,3]),
    c(variance[[3]][[2]]$IE0[-1,2] / variance[[3]][[2]]$IE0[-1,3]),
    c(variance[[3]][[2]]$IE1[-1,2] / variance[[3]][[2]]$IE1[-1,3]),
    c(variance[[4]][[1]]$DE[-1,2] / variance[[4]][[1]]$DE[-1,3]),
    c(variance[[4]][[1]]$IE0[-1,2] / variance[[4]][[1]]$IE0[-1,3]),
    c(variance[[4]][[1]]$IE1[-1,2] / variance[[4]][[1]]$IE1[-1,3]),
    c(variance[[4]][[2]]$DE[-1,2] / variance[[4]][[2]]$DE[-1,3]),
    c(variance[[4]][[2]]$IE0[-1,2] / variance[[4]][[2]]$IE0[-1,3]),
    c(variance[[4]][[2]]$IE1[-1,2] / variance[[4]][[2]]$IE1[-1,3]),
    c(variance[[5]][[1]]$DE[-1,2] / variance[[5]][[1]]$DE[-1,3]),
    c(variance[[5]][[1]]$IE0[-1,2] / variance[[5]][[1]]$IE0[-1,3]),
    c(variance[[5]][[1]]$IE1[-1,2] / variance[[5]][[1]]$IE1[-1,3]),
    c(variance[[5]][[2]]$DE[-1,2] / variance[[5]][[2]]$DE[-1,3]),
    c(variance[[5]][[2]]$IE0[-1,2] / variance[[5]][[2]]$IE0[-1,3]),
    c(variance[[5]][[2]]$IE1[-1,2] / variance[[5]][[2]]$IE1[-1,3]),    
    c(variance[[6]][[1]]$DE[-1,2] / variance[[6]][[1]]$DE[-1,3]),
    c(variance[[6]][[1]]$IE0[-1,2] / variance[[6]][[1]]$IE0[-1,3]),
    c(variance[[6]][[1]]$IE1[-1,2] / variance[[6]][[1]]$IE1[-1,3]),
    c(variance[[6]][[2]]$DE[-1,2] / variance[[6]][[2]]$DE[-1,3]),
    c(variance[[6]][[2]]$IE0[-1,2] / variance[[6]][[2]]$IE0[-1,3]),
    c(variance[[6]][[2]]$IE1[-1,2] / variance[[6]][[2]]$IE1[-1,3])
  )
  
  data.frame(
    effect = rep(c(rep("DE",7), rep("IE0",7), rep("IE1",7)), times = 12),
    perc.inter = inter.var,
    delta = delta,
    estimator = rep(c(rep("parametric", 21), rep("BART", 21)), times = 6)
  )
}

plotDeltaComps <- function(plot.name, delta.comps){
  # Plot comparing % variance attributed to interference and
  #   the estimated change in coverage across all simulation studies.
  # Input:
  #   plot.name: figure name
  #   delta.comps: data frame generated by 'deltaCoverage'
  # Output:
  #   Figure 12 from supplement

  png(
    file = plot.name, width = 6, height = 5,
    units = "in", res = 300, bg = "white"
  )
  
  par(mfrow = c(1,1))
  plot(delta.comps[,2], delta.comps[,3], pch = 20, cex = 0.4,
     bty = "l", xlab = "% variance attributed to interference", 
     ylab = "change in coverage")

  # add linear regression
  delta.corr <- lm(delta ~ perc.inter, data = delta.comps)
  abline(delta.corr, col = "red", lty = 3, lwd = 2)

  # add cubic smoothing spline
  delta.ss <- npreg::ss(delta.comps[,2], delta.comps[,3], m = 2, nknots = 5)
  lines(delta.ss$x, delta.ss$y, lty = 2, col = "darkblue", lwd = 2)

  # add legend
  legend("topleft", inset = 0.05, legend=c("OLS Regression", "Cubic Smoothing Spline"),
       col=c("red", "darkblue"), lty=3:2, lwd = 2, cex=0.8,
       box.lty=0)
  
  dev.off()
}
