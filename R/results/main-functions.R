
### main-functions.R
### Nathan Wikle

### Functions used to generate the results in the main text.

##############################################################
### 1. Plotting Functions
##############################################################

### Plot power plant facilities with scrubber status

scrubberPlot <- function(plot.name, coordinates, scrubbers){
  # Creates Figure 1a, plotting power plant facilities and their
  #   corresponding scrubber status.
  # Input:
  #   plot.name: the name of the figure (should end with .png)
  #   coordinates: sf object containing facility coordinates
  #   scrubbers: vector with facilities' scrubbers status
  # Output:
  #   PNG showing power plant facilities and scrubber status (Fig 1a).
  
  png(
    file = plot.name,
    height = 5.5, width = 7, units = "in", res = 300, bg = "white"
  )

  par(mar = c(0, 2, 3, 2))

  plot(
    coordinates,
    main = "",
    ylim = c(26.25, 40), pch = 20,
    col = "black"
  )
  map("state", add = TRUE)
  map(database = "state", regions = "texas", col = "blue", add = TRUE)
  plot(
    coordinates[which(scrubbers == 1), ],
    add = TRUE, pch = 20, cex = 1.1,
    col = "red"
  )

  legend("topleft",
    title = "Power Plants",
    legend = c("scrubbed", "unscrubbed"),
    bty = "o", pch = c(20, 20),
    pt.cex = c(1.1, 1.1),
    col = c("red", "black"),
    box.col = "black", bg = "white"
  )
  
  dev.off()

}

### Plot asthma rates

asthmaPlot <- function(plot.name, shp_file, tx_shp, type = "counts", geo = "zcta"){
  # Creates Figure 1b, plotting the pediatric asthma rates across TX ZCTAs.
  # Input:
  #   plot.name: the name of the figure (should end with .png)
  #   shp_file: sf object with ashthma rates for each ZCTA
  #   tx_shp: sf shape file of Texas
  #   type: type of output (e.g., "counts" or "rates")
  #   geo: geographic spatial resolution
  # Output:
  #   PNG with pediatric asthma rates in TX ZCTAs (Fig 1b).
  
  # interesection with texas
  plot_sf <- st_intersection(shp_file, tx_shp)

  geo_color = NA
  if (geo == "county"){
    geo_color = "white"
  }

  if (type == "counts") {
    # group by count totals
    plot_sf$asthma_groups <- base::cut(
      plot_sf$asthma_total,
      breaks = c(-Inf, 5, 25, 50, 100, 200, 500, 1000, Inf),
      labels = c("0-5", "6-25", "26-50", "51-100", "101-200", "201-500", "501-1000", "1000+")
    )

    asthma_map <- ggplot(data = plot_sf) +
      geom_sf(data = tx_shp$geometry, fill = "grey90") + 
      geom_sf(
        aes(fill = asthma_groups),
        lwd = 0,
        colour = geo_color
      ) +
      geom_sf(data = tx_shp$geometry, fill = NA) +
      scale_fill_manual(
        values = brewer.pal(8, name = "Blues"),
        name = "Total Hosp."
      ) + theme_void() +
      theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.05, .95),
        # change title placement
        plot.title = element_text(hjust = 0.05, vjust = -3)
      ) +
      ggtitle("Asthma hospitalizations (both inpatient and outpatient) in 2016.")
  } else if (type == "rates"){

    # create breaks based on rates
    # plot_sf$rate_breaks <- base::cut(
    #   plot_sf$asthma_rate,
    #   breaks = c(-Inf, 2, 4, 6, 8, 10, 25, 100, Inf),
    #   labels = c("[0, 2]", "(2, 4]", "(4, 6]", "(6, 8]", "(8, 10]", "(10, 25]", "(25, 100]", "100+")
    # )

    plot_sf$rate_breaks <- base::cut(
      plot_sf$asthma_rate,
      breaks = c(-Inf, 5, 10, 20, 100, Inf),
      labels = c("[0, 5]", "(5, 10]", "(10, 20]", "(20, 100]", "100+")
    )

    # # create breaks based on rates
    # plot_sf$rate_breaks <- base::cut(
    #   plot_sf$asthma_rates,
    #   breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, Inf),
    #   labels = c("[0, 0.2]", "(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", "(0.8, 1]", "(1, 1.2]", "(1.2, 1.4]")
    # )

    # create plot
    asthma_map <- ggplot(data = plot_sf) +
      geom_sf(data = tx_shp$geometry, fill = "grey80") +
      geom_sf(
        aes(fill = rate_breaks),
        lwd = 0,
        colour = geo_color
      ) +
      geom_sf(data = tx_shp$geometry, fill = NA) +
      scale_fill_manual(
        values = brewer.pal(5, name = "Reds"),
        name = "ED visits (per 1000)",
        na.value = "grey80", 
        na.translate = FALSE
      ) + theme_void() +
      theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.02, .93),
        # change title placement
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 18, face = "bold")
      ) # +
      ggtitle("Rate of asthma ED visits per 1000 people in 2016") 
  } else {
    stop("Error: 'type' must = either 'counts' or 'rates'")
  }

  # save ggplot
  ggsave(
    filename = plot.name,
    plot = asthma_map,
    width = 7,
    height = 5.5,
    units = "in",
    dpi = "print",
  )
}

### Sulfate plots

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

sulfateDataPlot <- function(
  plot.name, tx
){
  # Creates Figure 2a, plotting the 2016 annual sulfate data.
  # Input:
  #   plot.name: the name of the figure (should end with .png)
  #   tx: list containing SO4 raster, wind, and emissions data
  # Output:
  #   PNG of 2016 annual average SO4 concentrations (Fig 2a).
  
  png(
    file = plot.name,
    width = 7, height = 5.75, units = "in",
    res = 300, bg = "white" # "transparent",
  )

  uwind = tx$wind.small$wind.2016$uwind
  vwind = tx$wind.small$wind.2016$vwind 

  # coordinate system
  xy <- coordinates(uwind)

  # u and v wind values
  u.w.vals <- values(uwind)
  v.w.vals <- values(vwind)

  # create a raster for plotting
  n.row <- dim(uwind)[1]
  n.col <- dim(uwind)[2]

  row.spacing <- (1:(floor(n.row / 4)) * 4) - 1
  arrow.pts <- sapply(row.spacing, function(x) {
    (n.col * x) + (1:floor(n.col / 4) * 4)
  })

  par(
    xpd = FALSE,
    mar = c(2.5, 1, 2.5, 1),
    oma = c(0, 0, 0, 0)
  )
  
  em.mult <- 0.00004

  ### Actual SO4 surface
  plot(tx$so4,
    legend = FALSE, axes = FALSE,
    breaks = breakpoint.creation(
      tx$so4,
      n = 255,
      min.val = 0,
      max.val = 3
    ),
    col = rev(terrain.colors(255)),
    main = "", cex.main = 1.5,
    box = TRUE
  )

  map("state", add = TRUE)

  arrows(
    x0 = xy[arrow.pts, 1], y0 = xy[arrow.pts, 2],
    x1 = xy[arrow.pts, 1] + 0.3 * u.w.vals[arrow.pts],
    y1 = xy[arrow.pts, 2] + 0.3 * v.w.vals[arrow.pts],
    length = 0.03, col = "grey40"
  )

  plot(tx$so4,
    legend.only = TRUE,
    legend.width = 0.75, legend.shrink = 0.55,
    breaks = breakpoint.creation(
      tx$so4,
      n = 255,
      min.val = 0,
      max.val = 3
    ),
    col = rev(terrain.colors(255)),
    axis.args = list(
      at = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      cex.axis = 0.75
    ),
    legend.args =
      list(
        text = expression(paste(S * O[4], "  (", mu * g / m^3, ")", sep = "")),
        side = 4, font = 2, line = 2.75, cex = 1
      )
  )

  points(
    tx$em$em.usa$fac$Fac.Longitude,
    tx$em$em.usa$fac$Fac.Latitude,
    col = "black",
    pch = 16, cex = em.mult * tx$em$em.usa$fac$totSO2emissions * 12
  )

  points(
    tx$em$em.mx$longitude[164],
    tx$em$em.mx$latitude[164],
    col = "darkred",
    pch = 18, cex = 3 * tx$em$em.mx$total.SO2[164] / 200000
  )

  points(
    tx$em$em.mx$longitude[-c(164, 167)],
    tx$em$em.mx$latitude[-c(164, 167)],
    col = "black",
    pch = 16, cex = em.mult * tx$em$em.mx$total.SO2[-c(164, 167)]
  )

  # add legend depicting emissions size
  par(xpd = TRUE)
  legend("topleft",
    title = expression(paste("2016 ",
      S * O[2], " Emissions",
      sep = ""
    )),
    legend = c("25k tons", "50k tons", "200k tons"),
    bty = "o", pch = c(16, 16, 18),
    pt.cex = c(25000 * em.mult, 50000 * em.mult, 3),
    col = c("black", "black", "darkred"),
    box.col = "black", bg = "white"
  )

  dev.off()
}

expectedSO4Plot <- function(plot.name, mats, so4.data, so4.res){
  # Creates a figure showing the expected annual sulfate in 2016
  #   attributed to observed coal-fired power plants, using the 
  #   estimated posterior mean of the sulfate model parameters.
  # Input:
  #   plot.name: the figure name
  #   mats: advection-diffusion FVM matrices
  #   so4.data: sulfate data used in the analysis
  #   so4.res: mcmc samples from the analysis
  # Output:
  #   Figure showing the estimated annual sulfate on the same scale as the
  #     observed annual sulfate in 2016 (Figure 2b).

  # posterior mean of the sulfate model parameters
  burnin <- 25000
  theta_mu <- colMeans(so4.res$samples[-c(1:burnin), ])

  # convert to expected sulfate given observed emissions
  mu_hat <- meanSurface(
    theta = theta_mu,
    X_tx = so4.data$X$X.usa,
    X_mx = so4.data$X$X.mexico,
    mats = list(D = mats$D, C = mats$C)
  )
  
  # convert mu_hat to a raster
  mean_raster <- so4.data$so4
  values(mean_raster) <- mu_hat

  png(
    file = plot.name,
    width = 7, height = 5.75, 
    units = "in", res = 300, bg = "white"
  )

  par(
    xpd = FALSE,
    mar = c(2.5, 1, 2.5, 1),
    oma = c(0, 0, 0, 0)
  )

  ### predicted SO4 surface
  plot(
    mean_raster,
    legend = FALSE, axes = FALSE,
    breaks = breakpoint.creation(
      so4.data$so4,
      n = 255,
      min.val = 0,
      max.val = 3.0
    ),
    col = rev(terrain.colors(255)),
    main = "", cex.main = 1.5,
    box = TRUE
  )
  map("state", add = TRUE)

  plot(
    so4.data$so4,
    legend.only = TRUE,
    legend.width = 0.75, legend.shrink = 0.55,
    breaks = breakpoint.creation(
      so4.data$so4,
      n = 255,
      min.val = 0,
      max.val = 3.0
    ),
    col = rev(terrain.colors(255)),
    axis.args = list(
      at = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      cex.axis = 0.75
    ),
    legend.args = list(
      text = expression(paste(S * O[4], "  (", mu * g / m^3, ")", sep = "")),
      side = 4, font = 2, line = 2.75, cex = 1
    )
  )
  
  em.mult <- 0.00004

  points(
    so4.data$em$em.usa$fac$Fac.Longitude,
    so4.data$em$em.usa$fac$Fac.Latitude,
    col = "black",
    pch = 16, cex = em.mult * so4.data$em$em.usa$fac$totSO2emissions * 12
  )

  points(
    so4.data$em$em.mx$longitude[164],
    so4.data$em$em.mx$latitude[164],
    col = "darkred",
    pch = 18, cex = 3 * so4.data$em$em.mx$total.SO2[164] / 200000
  )

  points(
    so4.data$em$em.mx$longitude[-c(164, 167)],
    so4.data$em$em.mx$latitude[-c(164, 167)],
    col = "black",
    pch = 16, cex = em.mult * so4.data$em$em.mx$total.SO2[-c(164, 167)]
  )

  # add legend depicting emissions size
  par(xpd = TRUE)
  legend("topleft",
    title = expression(paste("2016 ",
      S * O[2], " Emissions",
      sep = ""
    )),
    legend = c("25k tons", "50k tons", "200k tons"),
    bty = "o", pch = c(16, 16, 18),
    pt.cex = c(25000 * em.mult, 50000 * em.mult, 2.25),
    col = c("black", "black", "darkred"),
    box.col = "black", bg = "white"
  )

  dev.off()
}

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

### Plot Treatment Levels

plotTreatment <- function(
    plot_name, 
    shp_file,
    tx_shp,
    trt = "Z", 
    title = "",
    fac_scrub = NULL,
    fac_noscrub = NULL
){
  # Creates Figure 3, plotting the treatment values for TX ZCTAs.
  # Input:
  #   plot.name: the name of the figure (should end with .png)
  #   shp_file: sf object with treatment levels for each ZCTA
  #   tx_shp: sf shape file of Texas
  #   trt: type of treatment to plot ('Z', 'G', or 'G sd')
  #   title: main figure title
  #   fac_scrub: add scrubbed facilities to the plot (if NULL, no facilities added)
  #   fac_noscrub: add unscrubbed facilities to the plot (if NULL, no facilities added)
  # Output:
  #   PNG with treatment levels across TX ZCTAs (Fig 3).
  
  # interesection with texas
  plot_sf <- st_intersection(shp_file, tx_shp)
  geo_color = NA
  
  if (trt == "Z"){
    # z plot
    trt_map <- ggplot(data = plot_sf) +
      geom_sf(data = tx_shp$geometry, fill = "grey85") +
      geom_sf(
        aes(fill = factor(1 - Z)), 
        lwd = 0,
        colour = geo_color
      ) + 
      geom_sf(data = tx_shp$geometry, fill = NA) +
      scale_fill_manual(
        values = c("#079aab", "#BC544B"), 
        labels = c("Scrubbed", "Unscrubbed"),
        name = "Direct Treatment (Z)"
      )
    
    if (!is.null(fac_scrub)){
      # add facilities to plot
      trt_map <- trt_map +
        geom_sf(
          data = fac_scrub,
          size = 3, shape = 20,
          col = "#e3980d"
        ) +
        geom_sf(
          data = fac_noscrub, 
          size = 3, shape = 20, 
          col = "black"
        )
    }
    
    trt_map <- trt_map +  theme_void() +
      theme(   
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.62, .97),
        # change title placement
        plot.title = element_text(hjust = 0.05, vjust = -2),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm')
      ) +
      ggtitle(title)
    
  } else if (trt == "G"){
    # g plot
    trt_map <- ggplot(data = plot_sf) +
      geom_sf(data = tx_shp$geometry, fill = "grey85") +
      geom_sf(
        aes(fill = G), 
        lwd = 0,
        colour = geo_color
      ) + 
      geom_sf(data = tx_shp$geometry, fill = NA) +
      scale_fill_carto_c(
        name = "Mean Upwind Treatment (G)",
        type = "diverging", palette = "Earth",
        limits = c(0,1),
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.vjust = 0.5, 
          title.hjust = 0
        )
      ) + theme_void() + 
      # scale_fill_continuous(type = "viridis",
      #   name = "Indirect Treatment (G)"
      # ) + theme_void() +
      theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.62, .97),
        # legend.position = c(0.03, .9),
        # change title placement
        plot.title = element_text(hjust = 0.05, vjust = -3),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm')
      ) +
      ggtitle(title)
  } else if (trt == "G sd"){
    
    # g plot
    trt_map <- ggplot(data = plot_sf) +
      geom_sf(data = tx_shp$geometry, fill = "grey85") +
      geom_sf(
        aes(fill = G_sd), 
        lwd = 0,
        colour = geo_color
      ) + 
      geom_sf(data = tx_shp$geometry, fill = NA) +
      scale_fill_carto_c(
        name = "Standard deviation of G",
        type = "divergent", palette = "Geyser",
        limits = c(0,0.2), 
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.vjust = 0.5, 
          title.hjust = 0
        )
      ) + theme_void() + 
      theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.62, .97),
        # change title placement
        plot.title = element_text(hjust = 0.05, vjust = -3),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm')
      ) +
      ggtitle(title)
    
  } else if (trt == "degree"){
    
    # g plot
    trt_map <- ggplot(data = plot_sf) +
      geom_sf(data = tx_shp$geometry, fill = "grey85") +
      geom_sf(
        aes(fill = degree), 
        lwd = 0,
        colour = geo_color
      ) + 
      geom_sf(data = tx_shp$geometry, fill = NA) +
      scale_fill_carto_c(
        name = "Degree",
        type = "sequential", palette = "SunsetDark",
        limits = c(0, 1)
      ) + theme_void() + 
      theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.03, .9),
        # change title placement
        plot.title = element_text(hjust = 0.05, vjust = -3)
      ) +
      ggtitle(title)
    
    
  } else {
    # incorrect trt value
    stop("This value of 'trt' is not accepted. Please specify 'Z', 'G', or 'G sd'.")
  }
  
  # save ggplot
  ggsave(
    filename = plot_name,
    plot = trt_map,
    width = 7,
    height = 7,
    units = "in",
    dpi = "screen",
  )
  
  # return(trt_map)
}

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

varComponents <- function(res, n_samps){
  # Estimates the total, between, and within variance for a given estimate.
  # Input:
  #   res: casual effect output from 'cut feedback' method
  #   n_samps: number of multiple imputation samples
  # Output:
  #   The between, within, and total variance of the estimate.
  
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
  # Calculate the proportion of variance attributed to uncertainty 
  #   in the interference structure.
  # Input:
  #   results: list of cut feedback results
  #   loglin: whether the output is from a loglinear (e.g., count) model
  # Output:
  #   List with between, within, and total variance for DE, IE0, and IE1.
  
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
    var.DE <- varComponents(res = mult * results[[s]]$bart.cut$DE, n_samps = 100)
    var.IE0 <- varComponents(res = mult * results[[s]]$bart.cut$IE0, n_samps = 100)
    var.IE1 <- varComponents(res = mult * results[[s]]$bart.cut$IE1, n_samps = 100)
    
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

plotVar <- function(prop.var, c1, c2){
  # Plot the proportion of variance attributed to uncertainty in the 
  #   interference structure as a function of g.
  # Input:
  #   prop.var: proportion of variance
  #   c1: color 1
  #   c2: color 2
  # Output:
  #   Plot with proportion of variance for parameteric and nonparametric 
  #     estimators.
  
  plot(prop.var[,1], prop.var[,2], type = "l", bty = "l",
       ylim = c(min(prop.var[,2:3]), max(prop.var[,2:3])),
       xlab = "g", ylab = "proportion of variance", col = c1, lwd = 2, 
       cex.axis = 1.25, cex.lab = 1.5)
  lines(prop.var[,1], prop.var[,3], type = "l", col = c2, lwd = 2, lty = 2)
  
}

plotAsthmaEffects <- function(
  plot.name = here::here("output", "main", "FIG5-asthma_results.png"), 
  pois.plugin, 
  pois.cut,
  bart.plugin, 
  bart.cut,
  prop.var.DE,
  prop.var.IE0,
  prop.var.IE1,
  y.lims.de, 
  y.lims.ie0,
  y.lims.ie1,
  cols
){
  # Creates a plot of the estimated causal effects from the asthma 
  #   analysis (Fig 5).
  # Input:
  #   plot.name: figure name
  #   pois.plugin: results from Poisson plugin inference
  #   pois.cut: results from Poisson cut inference 
  #   bart.plugin: results from log-linear BART plugin inference
  #   bart.cut: results from log-linear BART cut inference 
  #   prop.var.DE: proportion of variance ('DE')
  #   prop.var.IE0: proportion of variance ('IE0')
  #   prop.var.IE1: proportion of variance ('IE1')
  #   y.lims.de: y-axis limits ('DE')
  #   y.lims.ie0: y-axis limits ('IE0')
  #   y.lims.ie1: y-axis limits ('IE1')
  #   cols: vector with plot colors
  # Output:
  #   Figure 5 from main text.
  
  png(
    file = plot.name,
    height = 12, width = 15, units = "in", res = 300, bg = "white"
  )
  
  # par(mfrow = c(2,3), mar = c(6,4.1,5,9), oma = c(1,1,1,1))
  par(mfrow = c(3,3), mar = c(6,4.3,5,9), oma = c(1,1,1,1))
  
  ### DE(g)
  
  plotEffect(
    res.nounc = pois.plugin$DE$all,
    res.unc = pois.cut$DE$all,
    y.lims = y.lims.de, 
    c1 = cols[1], 
    c2 = cols[2],
    trans = 0.6,
    type = 2
  )
  
  title(main = "Poisson", line = 1, adj = 0, cex.main = 1.5)
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = -2 + de.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = -2, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -2 - de.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 3.2 + de.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 3.2, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 3.2 - de.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  
  plotEffect(
    res.nounc = bart.plugin$DE$all,
    res.unc = bart.cut$DE$all,
    y.lims = y.lims.de, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = 2
  )
  title(main = "log-linear BART", line = 1, adj = 0, cex.main = 1.5)
  
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 2.0 + de.ratio, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 2.0, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 0, 
        col = cols[4], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 0 - de.ratio, 
        col = cols[4], font = 2, outer = FALSE)
  
  
  plotVar(prop.var = prop.var.DE, c1 = cols[2], c2 = cols[3])
  title(main = "Proportion of Variance Attributed to Interference", 
        line = 1, adj = 0, cex.main = 1.55)
  
  mtext("Poisson", side = 4, las = 1, adj = 0, line = 0, at = 0.78, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("log-linear BART", side = 4, las = 1, adj = 0, line = 0, at = 0.43, 
        col = cols[3], font = 2, outer = FALSE)
  
  
  mtext("A.  Direct Effects:  DE(g)", line = -2, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  ### IE(0,g)
  plotEffect(
    res.nounc = pois.plugin$IE0$all,
    res.unc = pois.cut$IE0$all,
    y.lims = y.lims.ie0, 
    c1 = cols[1], 
    c2 = cols[2],
    trans = 0.6,
    type = 2
  )
  
  title(main = "Poisson", line = 1, adj = 0, cex.main = 1.5)
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 12.5 + ie0.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 12.5,
        col = cols[1], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 12.5 - ie0.ratio,
        col = cols[1], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 1.5 + ie0.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 1.5, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 1.5 - ie0.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  
  plotEffect(
    res.nounc = bart.plugin$IE0$all,
    res.unc = bart.cut$IE0$all,
    y.lims = y.lims.ie0, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = 2
  )
  title(main = "log-linear BART", line = 1, adj = 0, cex.main = 1.5)
  
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 2 + ie0.ratio, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 2, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = -0.6, 
        col = cols[4], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -0.6 - ie0.ratio, 
        col = cols[4], font = 2, outer = FALSE)
  
  
  plotVar(prop.var = prop.var.IE0[-1,], c1 = cols[2], c2 = cols[3])
  title(main = "Proportion of Variance Attributed to Interference", 
        line = 1, adj = 0, cex.main = 1.55)
  
  mtext("Poisson", side = 4, las = 1, adj = 0, line = 0, at = 0.95, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("log-linear BART", side = 4, las = 1, adj = 0, line = 0, at = 0.45, 
        col = cols[3], font = 2, outer = FALSE)
  
  
  # mtext("B.  Upwind Effects:  IE(0,g)", line = -31.1, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  mtext("B.  Upwind Effects:  IE(0,g)", line = -31.4, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  ### IE(1,g)
  plotEffect(
    res.nounc = pois.plugin$IE1$all,
    res.unc = pois.cut$IE1$all,
    y.lims = y.lims.ie1, 
    c1 = cols[1], 
    c2 = cols[2],
    trans = 0.6, 
    type = 2
  )
  
  title(main = "Poisson", line = 1, adj = 0, cex.main = 1.5)
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 12.5 + ie1.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 12.5, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 12.5 - ie1.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 5 + ie1.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 5, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 5 - ie1.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  
  plotEffect(
    res.nounc = bart.plugin$IE1$all,
    res.unc = bart.cut$IE1$all,
    y.lims = y.lims.ie1, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = 2
  )
  title(main = "log-linear BART", line = 1, adj = 0, cex.main = 1.5)
  
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 2.7 + ie1.ratio, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 2.7, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = -1, 
        col = cols[4], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -1 - ie1.ratio, 
        col = cols[4], font = 2, outer = FALSE)
  
  plotVar(prop.var = prop.var.IE1[-1,], c1 = cols[2], c2 = cols[3])
  title(main = "Proportion of Variance Attributed to Interference", 
        line = 1, adj = 0, cex.main = 1.55)
  
  mtext("Poisson", side = 4, las = 1, adj = 0, line = 0, at = 0.95, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("log-linear BART", side = 4, las = 1, adj = 0, line = 0, at = 0.44, 
        col = cols[3], font = 2, outer = FALSE)
  
  mtext("C.  Upwind Effects:  IE(1,g)", line = -61, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  dev.off()
}  


plotMedicareEffects <- function(
  plot.name = here::here("output", "main", "FIG5-asthma_results.png"), 
  pois.plugin, 
  pois.cut,
  bart.plugin, 
  bart.cut,
  prop.var.DE,
  prop.var.IE0,
  prop.var.IE1,
  y.lims.de, 
  y.lims.ie0,
  y.lims.ie1,
  cols
){
  # Creates a plot of the estimated causal effects from the medicare 
  #   analysis (Fig 4).
  # Input:
  #   plot.name: figure name
  #   pois.plugin: results from Poisson plugin inference
  #   pois.cut: results from Poisson cut inference 
  #   bart.plugin: results from log-linear BART plugin inference
  #   bart.cut: results from log-linear BART cut inference 
  #   prop.var.DE: proportion of variance ('DE')
  #   prop.var.IE0: proportion of variance ('IE0')
  #   prop.var.IE1: proportion of variance ('IE1')
  #   y.lims.de: y-axis limits ('DE')
  #   y.lims.ie0: y-axis limits ('IE0')
  #   y.lims.ie1: y-axis limits ('IE1')
  #   cols: vector with plot colors
  # Output:
  #   Figure 4 from main text.
  
  png(
    file = plot.name,
    height = 12, width = 15, units = "in", res = 300, bg = "white"
  )
  
  # par(mfrow = c(2,3), mar = c(6,4.1,5,9), oma = c(1,1,1,1))
  par(mfrow = c(3,3), mar = c(6,4.3,5,9), oma = c(1,1,1,1))
  
  ### DE(g)
  
  plotEffect(
    res.nounc = pois.plugin$DE$all,
    res.unc = pois.cut$DE$all,
    y.lims = y.lims.de, 
    c1 = cols[1], 
    c2 = cols[2],
    trans = 0.6,
    type = 1
  )
  
  title(main = "Poisson", line = 1, adj = 0, cex.main = 1.5)
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 3.7 + de.ratio * 2, 
      col = cols[1], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 3.7 + de.ratio, 
      col = cols[1], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 3.7, 
      col = cols[1], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 0, 
      col = cols[2], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 0 - de.ratio, 
      col = cols[2], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 0 - 2 * de.ratio, 
      col = cols[2], font = 2, outer = FALSE)
  
  plotEffect(
    res.nounc = bart.plugin$DE$all,
    res.unc = bart.cut$DE$all,
    y.lims = y.lims.de, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = 1
  )
  title(main = "log-linear BART", line = 1, adj = 0, cex.main = 1.5)
  
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 7 + de.ratio / 2, 
      col = cols[3], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 7 - de.ratio / 2, 
      col = cols[3], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 0 + de.ratio / 2, 
      col = cols[4], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 0 - de.ratio / 2, 
      col = cols[4], font = 2, outer = FALSE)
  
  
  plotVar(prop.var = prop.var.DE, c1 = cols[2], c2 = cols[3])
  title(main = "Proportion of Variance Attributed to Interference", 
        line = 1, adj = 0, cex.main = 1.55)
  
  mtext("Poisson", side = 4, las = 1, adj = 0, line = 0, at = 0.26, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("log-linear BART", side = 4, las = 1, adj = 0, line = 0, at = 0.175, 
        col = cols[3], font = 2, outer = FALSE)
  
  
  mtext("A.  Direct Effects:  DE(g)", line = -2, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  ### IE(0,g)
  plotEffect(
    res.nounc = pois.plugin$IE0$all,
    res.unc = pois.cut$IE0$all,
    y.lims = y.lims.ie0, 
    c1 = cols[1], 
    c2 = cols[2],
    trans = 0.6,
    type = 1
  )
  
  title(main = "Poisson", line = 1, adj = 0, cex.main = 1.5)
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = -21 + ie0.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = -21, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -21 - ie0.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 0 + 2 * ie0.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 0 + ie0.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 0, 
        col = cols[2], font = 2, outer = FALSE)
  
  plotEffect(
    res.nounc = bart.plugin$IE0$all,
    res.unc = bart.cut$IE0$all,
    y.lims = y.lims.ie0, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = 1
  )
  title(main = "log-linear BART", line = 1, adj = 0, cex.main = 1.5)
  
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 20 + ie0.ratio / 2, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 20 - ie0.ratio / 2, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = -16 + ie0.ratio / 2, 
        col = cols[4], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -16 - ie0.ratio / 2, 
        col = cols[4], font = 2, outer = FALSE)
  
  
  plotVar(prop.var = prop.var.IE0[-1,], c1 = cols[2], c2 = cols[3])
  title(main = "Proportion of Variance Attributed to Interference", 
        line = 1, adj = 0, cex.main = 1.55)
  
  mtext("Poisson", side = 4, las = 1, adj = 0, line = 0, at = 0.47, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("log-linear BART", side = 4, las = 1, adj = 0, line = 0, at = 0.25, 
        col = cols[3], font = 2, outer = FALSE)
  
  # mtext("B.  Upwind Effects:  IE(0,g)", line = -31.1, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  mtext("B.  Upwind Effects:  IE(0,g)", line = -31.4, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  ### IE(1,g)
  plotEffect(
    res.nounc = pois.plugin$IE1$all,
    res.unc = pois.cut$IE1$all,
    y.lims = y.lims.ie1, 
    c1 = cols[1], 
    c2 = cols[2],
    trans = 0.6,
    type = 1
  )
  
  title(main = "Poisson", line = 1, adj = 0, cex.main = 1.5)
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = -22 + ie1.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = -22, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -22 - ie1.ratio, 
        col = cols[1], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = 0 + 2 * ie1.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("interference", side = 4, las = 1, adj = 0, line = 0, at = 0 + ie1.ratio, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 0, 
        col = cols[2], font = 2, outer = FALSE)
  
  plotEffect(
    res.nounc = bart.plugin$IE1$all,
    res.unc = bart.cut$IE1$all,
    y.lims = y.lims.ie1, 
    c1 = cols[3],
    c2 = cols[4],
    trans = 0.6,
    type = 1
  )
  title(main = "log-linear BART", line = 1, adj = 0, cex.main = 1.5)
  
  mtext("with", side = 4, las = 1, adj = 0, line = 0, at = 20 + ie1.ratio / 2, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = 20 - ie1.ratio / 2, 
        col = cols[3], font = 2, outer = FALSE)
  mtext("without", side = 4, las = 1, adj = 0, line = 0, at = -16 + ie1.ratio / 2, 
        col = cols[4], font = 2, outer = FALSE)
  mtext("uncertainty", side = 4, las = 1, adj = 0, line = 0, at = -16 - ie1.ratio / 2, 
        col = cols[4], font = 2, outer = FALSE)
  
  plotVar(prop.var = prop.var.IE1[-1,], c1 = cols[2], c2 = cols[3])
  title(main = "Proportion of Variance Attributed to Interference", 
        line = 1, adj = 0, cex.main = 1.55)
  
  mtext("Poisson", side = 4, las = 1, adj = 0, line = 0, at = 0.67, 
        col = cols[2], font = 2, outer = FALSE)
  mtext("log-linear BART", side = 4, las = 1, adj = 0, line = 0, at = 0.24, 
        col = cols[3], font = 2, outer = FALSE)
  
  mtext("C.  Upwind Effects:  IE(1,g)", line = -61, adj = 0, outer = TRUE, font = 2, cex = 1.5)
  
  dev.off()
}  
