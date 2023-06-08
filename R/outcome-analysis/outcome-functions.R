### outcome-functions.R
### Nathan Wikle

########################################################################
### 1. Data cleaning and set-up
########################################################################

### A. Convert Lat-Long emissions to raster locations

createX <- function(r, emissions){
  # Create the X design matrix for the pollution data.
  # Input:
  #   r: a single raster object
  #   emissions: a data frame with power plant emissions data.
  # Output:
  #   A design vector, X, representing the total amount of emissions in 
  #     each grid of the raster.
  
  # number of cells in rectangle
  n.cells <- length(values(r))
  
  # create an SO2 variable using the emissions data
  em.vals <- rep(0, n.cells)
  em.inds <- raster::cellFromXY(r, cbind(emissions$Fac.Longitude, emissions$Fac.Latitude))
  
  # make sure unique inds are
  unique.inds <- unique(em.inds)
  
  for(curr.ind in unique.inds){
    spec.points <- em.inds == curr.ind
    em.vals[curr.ind] <- sum(emissions$totSO2emissions[spec.points])
  }
  
  # combine all values into a single vector
  X <- em.vals
  
  # return X
  return(X)
}

### B. Create Finite Volume Method structures for given raster

FVMInds <- function(dims, bound.cond){
  # Create list denoting index for interior and boundary of raster with given
  #   spatial extent. This is used to create the FVM approximation of advection
  #   due to wind.
  # Input:
  #   dims: dimension of raster
  #   bound.cond: type of boundary condition (insulated or periodic).
  # Output:
  #   List with index of interior and boundary points.
  
  # rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]
  
  #######################
  ### interior points ###
  #######################
  interior <- matrix(NA_real_, nrow = (n.c - 2) * (n.r - 2), ncol = 5)
  colnames(interior) <- c("P", "N", "S", "E", "W")
  
  if (bound.cond == "insulated"){
    # create a diffusion matrix with insulated boundaries
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){
        
        # index in original scheme
        orig.ind <- (i - 1) * n.c + j
        
        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1
        
        cur.ind <- cur.ind + 1
      }
    }
    
    ### west boundary points
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(west.bound) <- c("P", "N", "S", "E")
    
    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### east boundary
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(east.bound) <- c("P", "N", "S", "W")
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      
      # index in original scheme
      orig.ind <- i * n.c
      
      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### north boundary
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(north.bound) <- c("P", "S", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j
      
      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + n.c
      north.bound[cur.ind, 3] <- orig.ind + 1
      north.bound[cur.ind, 4] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### south boundary
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(south.bound) <- c("P", "N", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)
      
      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind + 1
      south.bound[cur.ind, 4] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ### corners
    nw <- c(1, NA_real_, 1 + n.c, 2, NA_real_)
    ne <- c(n.c, NA_real_, 2 * n.c, NA_real_, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, NA_real_, (n.r - 1) * n.c + 2, NA_real_)
    se <- c(n.r * n.c, (n.r - 1) * n.c, NA_real_, NA_real_, n.r * n.c - 1)
    
    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")
    
  } else if (bound.cond == "periodic"){
    
    # create a diffusion matrix with periodic boundaries
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){
        
        # index in original scheme
        orig.ind <- (i - 1) * n.c + j
        
        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1
        
        cur.ind <- cur.ind + 1
      }
    }
    
    ############################
    ### west boundary points ###
    ############################
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(west.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1
      west.bound[cur.ind, 5] <- orig.ind + (n.c - 1)
      
      cur.ind <- cur.ind + 1
    }
    
    #####################
    ### east boundary ###
    #####################
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(east.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (i in 2:(n.r-1)){
      
      # index in original scheme
      orig.ind <- i * n.c
      
      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - (n.c - 1)
      east.bound[cur.ind, 5] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ######################
    ### north boundary ###
    ######################
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(north.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j
      
      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + (n.r - 1) * n.c
      north.bound[cur.ind, 3] <- orig.ind + n.c
      north.bound[cur.ind, 4] <- orig.ind + 1
      north.bound[cur.ind, 5] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ######################
    ### south boundary ###
    ######################
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(south.bound) <- c("P", "N", "S", "E", "W")
    
    cur.ind <- 1
    for (j in 2:(n.c-1)){
      
      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)
      
      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind - (n.r - 1) * n.c
      south.bound[cur.ind, 4] <- orig.ind + 1
      south.bound[cur.ind, 5] <- orig.ind - 1
      
      cur.ind <- cur.ind + 1
    }
    
    ###############
    ### corners ###
    ###############
    nw <- c(1, 1 + (n.r - 1) * n.c, 1 + n.c, 2, n.c)
    ne <- c(n.c, n.c * n.r, 2 * n.c, 1, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, 1, (n.r - 1) * n.c + 2, n.r * n.c)
    se <- c(n.r * n.c, (n.r - 1) * n.c, n.c, 1 + (n.r - 1) * n.c, n.r * n.c - 1)
    
    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")
  }
  
  # combine all indice matrices into a single list
  indices <- list(interior = interior,
                  north = north.bound,
                  south = south.bound,
                  east = east.bound,
                  west = west.bound,
                  corners = corners)
  
  # return indices
  return(indices)
}

diffFVM <- function(dims, inds, bound.cond){
  # Create finite volume matrix approximation to homogeneous diffusion.
  # Input:
  #   dims: dimension of raster
  #   inds: indices of interior/boundary points, created with FVMInds
  #   bound.cond: boundary conditions (insulated, periodic, Neumann)
  # Output:
  #   Sparse matrix used for discrete approximation of homogenous diffusion.

  # rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]

  # total number of cells
  n.cells <- n.r * n.c

  if (bound.cond == "insulated"){

    # i, j, and x values for interior points
    n.interior <- nrow(inds$interior)
    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(rep(4, n.interior), rep(-1, n.interior * 4)))

    # i, j, and x values for west boundary points
    n.west <- nrow(inds$west)
    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west),
                         c(rep(3, n.west), rep(-1, n.west * 3)))

    # i, j, and x values for east boundary points
    n.east <- nrow(inds$east)
    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east),
                         c(rep(3, n.east), rep(-1, n.east * 3)))

    # i, j, and x values for north boundary points
    n.north <- nrow(inds$north)
    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north),
                          c(rep(3, n.north), rep(-1, n.north * 3)))

    # i, j, and x values for south boundary points
    n.south <- nrow(inds$south)
    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south),
                          c(rep(3, n.south), rep(-1, n.south * 3)))

    # i, j, and x values for corner points

    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(rep(2, 4), rep(-1, 8)))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

  } else if (bound.cond == "periodic"){

    # i, j, and x values for interior points
    n.interior <- nrow(inds$interior)
    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(rep(4, n.interior), rep(-1, n.interior * 4)))


    # i, j, and x values for west boundary points
    n.west <- nrow(inds$west)
    sparse.west <- cbind(rep(inds$west[,1], 5),
                         as.vector(inds$west),
                         c(rep(4, n.west), rep(-1, n.west * 4)))

    # i, j, and x values for east boundary points
    n.east <- nrow(inds$east)
    sparse.east <- cbind(rep(inds$east[,1], 5),
                         as.vector(inds$east),
                         c(rep(4, n.east), rep(-1, n.east * 4)))

    # i, j, and x values for north boundary points
    n.north <- nrow(inds$north)
    sparse.north <- cbind(rep(inds$north[,1], 5),
                          as.vector(inds$north),
                          c(rep(4, n.north), rep(-1, n.north * 4)))

    # i, j, and x values for south boundary points
    n.south <- nrow(inds$south)
    sparse.south <- cbind(rep(inds$south[,1], 5),
                          as.vector(inds$south),
                          c(rep(4, n.south), rep(-1, n.south * 4)))

    # i, j, and x values for corner points
    sparse.corner <- cbind(rep(inds$corners[,1], 5),
                           as.vector(inds$corners),
                           c(rep(4, 4), rep(-1, 4 * 4)))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

  } else if (bound.cond == "Neumann"){

    # i, j, and x values for interior points
    n.interior <- nrow(inds$interior)
    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(rep(4, n.interior), rep(-1, n.interior * 4)))

    # i, j, and x values for west boundary points
    n.west <- nrow(inds$west)
    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west),
                         c(rep(4, n.west), rep(-1, n.west * 2), rep(-2, n.west)))

    # i, j, and x values for east boundary points
    n.east <- nrow(inds$east)
    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east),
                         c(rep(4, n.east), rep(-1, n.east * 2), rep(-2, n.east)))

    # i, j, and x values for north boundary points
    n.north <- nrow(inds$north)
    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north),
                          c(rep(4, n.north), rep(-2, n.north), rep(-1, n.north * 2)))

    # i, j, and x values for south boundary points
    n.south <- nrow(inds$south)
    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south),
                          c(rep(4, n.south), rep(-2, n.south), rep(-1, n.south * 2)))

    # i, j, and x values for corner points
    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(rep(4, 4), rep(-2, 8)))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)
  }

  # quickly fill the sparse matrix
  D.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
                        dims = c(n.cells, n.cells))

  # return the sparse diffusion matrix
  return(D.mat)
}

FVMInds <- function(dims, bound.cond){
  # Create list denoting index for interior and boundary of raster with given
  #   spatial extent. This is used to create the FVM approximation of advection
  #   due to wind.
  # Input:
  #   dims: dimension of raster
  #   bound.cond: type of boundary condition (insulated or periodic).
  # Output:
  #   List with index of interior and boundary points.

  # rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]

  #######################
  ### interior points ###
  #######################
  interior <- matrix(NA_real_, nrow = (n.c - 2) * (n.r - 2), ncol = 5)
  colnames(interior) <- c("P", "N", "S", "E", "W")

  if (bound.cond == "insulated"){
    # create a diffusion matrix with insulated boundaries

    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){

        # index in original scheme
        orig.ind <- (i - 1) * n.c + j

        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1

        cur.ind <- cur.ind + 1
      }
    }

    ### west boundary points
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(west.bound) <- c("P", "N", "S", "E")

    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1

      cur.ind <- cur.ind + 1
    }

    ### east boundary
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 4)
    colnames(east.bound) <- c("P", "N", "S", "W")

    cur.ind <- 1
    for (i in 2:(n.r-1)){

      # index in original scheme
      orig.ind <- i * n.c

      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ### north boundary
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(north.bound) <- c("P", "S", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j

      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + n.c
      north.bound[cur.ind, 3] <- orig.ind + 1
      north.bound[cur.ind, 4] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ### south boundary
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 4)
    colnames(south.bound) <- c("P", "N", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)

      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind + 1
      south.bound[cur.ind, 4] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ### corners
    nw <- c(1, NA_real_, 1 + n.c, 2, NA_real_)
    ne <- c(n.c, NA_real_, 2 * n.c, NA_real_, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, NA_real_, (n.r - 1) * n.c + 2, NA_real_)
    se <- c(n.r * n.c, (n.r - 1) * n.c, NA_real_, NA_real_, n.r * n.c - 1)

    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")

  } else if (bound.cond == "periodic"){

    # create a diffusion matrix with periodic boundaries

    cur.ind <- 1
    for (i in 2:(n.r-1)){
      for (j in 2:(n.c - 1)){

        # index in original scheme
        orig.ind <- (i - 1) * n.c + j

        # indices (with neighbors) for interior
        interior[cur.ind, 1] <- orig.ind
        interior[cur.ind, 2] <- orig.ind - n.c
        interior[cur.ind, 3] <- orig.ind + n.c
        interior[cur.ind, 4] <- orig.ind + 1
        interior[cur.ind, 5] <- orig.ind - 1

        cur.ind <- cur.ind + 1
      }
    }

    ############################
    ### west boundary points ###
    ############################
    west.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(west.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (i in 2:(n.r - 1)){
      # index in original scheme
      orig.ind <- (i - 1) * n.c + 1
      # indices (with neighbors) for west boundary
      west.bound[cur.ind, 1] <- orig.ind
      west.bound[cur.ind, 2] <- orig.ind - n.c
      west.bound[cur.ind, 3] <- orig.ind + n.c
      west.bound[cur.ind, 4] <- orig.ind + 1
      west.bound[cur.ind, 5] <- orig.ind + (n.c - 1)

      cur.ind <- cur.ind + 1
    }

    #####################
    ### east boundary ###
    #####################
    east.bound <- matrix(NA_real_, nrow = (n.r - 2), ncol = 5)
    colnames(east.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (i in 2:(n.r-1)){

      # index in original scheme
      orig.ind <- i * n.c

      # indices (with neighbors) for interior
      east.bound[cur.ind, 1] <- orig.ind
      east.bound[cur.ind, 2] <- orig.ind - n.c
      east.bound[cur.ind, 3] <- orig.ind + n.c
      east.bound[cur.ind, 4] <- orig.ind - (n.c - 1)
      east.bound[cur.ind, 5] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ######################
    ### north boundary ###
    ######################
    north.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(north.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j

      # indices (with neighbors) for interior
      north.bound[cur.ind, 1] <- orig.ind
      north.bound[cur.ind, 2] <- orig.ind + (n.r - 1) * n.c
      north.bound[cur.ind, 3] <- orig.ind + n.c
      north.bound[cur.ind, 4] <- orig.ind + 1
      north.bound[cur.ind, 5] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ######################
    ### south boundary ###
    ######################
    south.bound <- matrix(NA_real_, nrow = (n.c - 2), ncol = 5)
    colnames(south.bound) <- c("P", "N", "S", "E", "W")

    cur.ind <- 1
    for (j in 2:(n.c-1)){

      # index in original scheme
      orig.ind <- j + n.c * (n.r - 1)

      # indices (with neighbors) for interior
      south.bound[cur.ind, 1] <- orig.ind
      south.bound[cur.ind, 2] <- orig.ind - n.c
      south.bound[cur.ind, 3] <- orig.ind - (n.r - 1) * n.c
      south.bound[cur.ind, 4] <- orig.ind + 1
      south.bound[cur.ind, 5] <- orig.ind - 1

      cur.ind <- cur.ind + 1
    }

    ###############
    ### corners ###
    ###############
    nw <- c(1, 1 + (n.r - 1) * n.c, 1 + n.c, 2, n.c)
    ne <- c(n.c, n.c * n.r, 2 * n.c, 1, n.c - 1)
    sw <- c((n.r - 1) * n.c + 1, (n.r - 2) * n.c + 1, 1, (n.r - 1) * n.c + 2, n.r * n.c)
    se <- c(n.r * n.c, (n.r - 1) * n.c, n.c, 1 + (n.r - 1) * n.c, n.r * n.c - 1)

    corners <- rbind(nw, ne, sw, se)
    colnames(corners) <- c("P", "N", "S", "E", "W")
  }

  # combine all indice matrices into a single list
  indices <- list(interior = interior,
                  north = north.bound,
                  south = south.bound,
                  east = east.bound,
                  west = west.bound,
                  corners = corners)

  # return indices
  return(indices)
}

interpWind <- function(r.old, r.new, direction, with.bound = TRUE){
  # Linear interpolation of wind velocity across surface boundaries.
  # Input:
  #   r.old: wind raster with larger spatial extent.
  #   r.new: wind raster with smaller spatial extent.
  #   direction: direction of wind (either "u" or "v").
  #   with.bound: boolean determining if velocities should be saved on
  #     the boundary.
  # Output:
  #   Matrix with interpolated matrix velocities.

  # grab xy coordinates from rasters
  xy.old <- coordinates(r.old)
  xy.new <- coordinates(r.new)

  # dimensions of rasters (n.rows by n.cols)
  dim.old <- dim(r.old)[1:2]
  dim.new <- dim(r.new)[1:2]

  # cell indices for trimmed case
  cell.inds <- cellFromXY(r.old, xy.new)

  # indices for the top boundary
  top.bound <- (cell.inds[1] - dim.old[2]) + (-1):(dim.new[2])

  # cell indices with boundaries
  all.bounds <- top.bound

  for (k in 1:(dim.new[1] + 1)){
    next.bound <- top.bound + (k * dim.old[2])
    all.bounds <- c(all.bounds, next.bound)
  }

  # wind values for all points (including boundaries)
  orig.vals <- values(r.old)[all.bounds]
  # matrix of face value velocities
  face.vals <- matrix(NA_real_, nrow = length(all.bounds), ncol = 2)

  # new dimensions
  n.rows <- dim.new[1] + 2
  n.cols <- dim.new[2] + 2

  if (direction == "u"){
    ### west to east

    # name columns appropriately
    colnames(face.vals) <- c("west", "east")

    for(i in 1:n.rows){
      for(j in 1:n.cols){

        # current indice
        cur.ind <- (i - 1) * n.cols + j

        if (j == 1){
          # no west face velocity
          west <- 0
          east <- sum(orig.vals[cur.ind + 0:1]) / 2
        } else if (j == n.cols){
          # no east face velocity
          west <-  sum(orig.vals[cur.ind + -1:0]) / 2
          east <- 0
        } else {
          # interior point
          west <- sum(orig.vals[cur.ind + -1:0]) / 2
          east <- sum(orig.vals[cur.ind + 0:1]) / 2
        }

        # add values to boundary matrix
        face.vals[cur.ind,] <- c(west, east)
      }
    }
  } else if (direction == "v"){
    # south to north
    # name columns appropriately
    colnames(face.vals) <- c("south", "north")

    for(i in 1:n.rows){
      for(j in 1:n.cols){

        # current indice
        cur.ind <- (i - 1) * n.cols + j

        if (i == 1){
          # no north face velocity
          north <- 0
          south <- (orig.vals[cur.ind] + orig.vals[cur.ind + n.cols]) / 2
        } else if (i == n.rows){
          # no south face velocity
          north <-  (orig.vals[cur.ind] + orig.vals[cur.ind - n.cols]) / 2
          south <- 0
        } else {
          # interior point
          north <-  (orig.vals[cur.ind] + orig.vals[cur.ind - n.cols]) / 2
          south <- (orig.vals[cur.ind] + orig.vals[cur.ind + n.cols]) / 2
        }

        # add values to boundary matrix
        face.vals[cur.ind,] <- c(south, north)
      }
    }
  }

  if (with.bound){
    # return face.vals as is
    return(face.vals)

  } else {
    # remove boundary points from face.vals
    bound.points <- c(1:n.cols,
                      (1:(n.rows-2) * n.cols + 1),
                      (2:(n.rows-1) * n.cols),
                      ((n.rows-1) * n.cols + 1:n.cols))

    return(face.vals[-bound.points,])
  }
}

windConvFVM <- function(dims, velocity.faces, inds, bc, bound = FALSE){
  # Create an advection matrix using the Finite Volume Method (FVM).
  # Input:
  #   dims: dimension of the raster
  #   velocity.faces: interpolated velocities on grid edges (from interpWind)
  #   inds: matrix indices of raster grid cells (from FVMinds).
  #   bc: boundary condition
  #   bound: not used?
  # Output:
  #   Sparse FVM matrix representation of advection due to wind.

  ## rows and column totals
  n.r <- dims[1]
  n.c <- dims[2]

  # total number of cells
  n.cells <- n.r * n.c

  if (bc == "insulated"){

    ### i, j, and x values for interior points
    n.interior <- nrow(inds$interior)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$interior[,1], 2] -
      velocity.faces[inds$interior[,1], 1]
    F.diff.ns <- velocity.faces[inds$interior[,1], 4] -
      velocity.faces[inds$interior[,1], 3]
    F.diff.comb <- F.diff.ew + F.diff.ns

    a.w <- velocity.faces[inds$interior[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$interior[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$interior[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$interior[,1], 4]; a.n[a.n < 0] <- 0

    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(F.diff.comb + a.n + a.s + a.e + a.w, -a.n, -a.s, -a.e, -a.w))


    ### i, j, and x values for west boundary points
    n.west <- nrow(inds$west)

    F.diff.ns <- velocity.faces[inds$west[,1], 4] -
      velocity.faces[inds$west[,1], 3]
    F.diff.comb <- F.diff.ns + velocity.faces[inds$west[,1], 2]

    a.e <- -velocity.faces[inds$west[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$west[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- velocity.faces[inds$west[,1], 4]; a.n[a.n > 0] <- 0

    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west[,-5]),
                         c(F.diff.comb + a.n + a.s + a.e, -a.n, -a.s, -a.e))

    ### i, j, and x values for east boundary points
    n.east <- nrow(inds$east)

    # coefficient calculations
    F.diff.ns <- velocity.faces[inds$east[,1], 4] -
      velocity.faces[inds$east[,1], 3]
    F.diff.comb <- F.diff.ns - velocity.faces[inds$east[,1], 1]

    a.w <- velocity.faces[inds$east[,1], 1]; a.w[a.w < 0] <- 0
    a.s <- velocity.faces[inds$east[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$east[,1], 4]; a.n[a.n < 0] <- 0

    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east[,-4]),
                         c(F.diff.comb + a.n + a.s + a.w, -a.n, -a.s, -a.w))

    ### i, j, and x values for north boundary points
    n.north <- nrow(inds$north)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$north[,1], 2] -
      velocity.faces[inds$north[,1], 1]
    F.diff.comb <- F.diff.ew - velocity.faces[inds$north[,1], 3]

    a.w <- velocity.faces[inds$north[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$north[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$north[,1], 3]; a.s[a.s < 0] <- 0

    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north[,-2]),
                          c(F.diff.comb + a.s + a.e + a.w, -a.s, -a.e, -a.w))

    ### i, j, and x values for south boundary points
    n.south <- nrow(inds$south)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$south[,1], 2] -
      velocity.faces[inds$south[,1], 1]
    F.diff.comb <- F.diff.ew + velocity.faces[inds$south[,1], 4]

    a.w <- velocity.faces[inds$south[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$south[,1], 2]; a.e[a.e < 0] <- 0
    a.n <- -velocity.faces[inds$south[,1], 4]; a.n[a.n < 0] <- 0

    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south[,-3]),
                          c(F.diff.comb + a.n + a.e + a.w, -a.n, -a.e, -a.w))

    ### i, j, and x values for corner points

    # northwest values
    nw.ae <- max(c(0, -velocity.faces[inds$corners[1,1], 2]))
    nw.as <- max(c(0, velocity.faces[inds$corners[1,1], 3]))
    nw.ap <- nw.ae + nw.as + velocity.faces[inds$corners[1,1], 2] -
      velocity.faces[inds$corners[1,1], 3]

    # northeast values
    ne.aw <- max(c(0, velocity.faces[inds$corners[2,1], 1]))
    ne.as <- max(c(0, velocity.faces[inds$corners[2,1], 3]))
    ne.ap <- ne.aw + ne.as - velocity.faces[inds$corners[2,1], 1] -
      velocity.faces[inds$corners[2,1], 3]

    # southwest values
    sw.ae <- max(c(0, -velocity.faces[inds$corners[3,1], 2]))
    sw.an <- max(c(0, -velocity.faces[inds$corners[3,1], 4]))
    sw.ap <- sw.ae + sw.an + velocity.faces[inds$corners[3,1], 2] +
      velocity.faces[inds$corners[3,1], 4]

    # southeast values
    se.aw <- max(c(0, velocity.faces[inds$corners[4,1], 1]))
    se.an <- max(c(0, -velocity.faces[inds$corners[4,1], 4]))
    se.ap <- se.aw + se.an - velocity.faces[inds$corners[4,1], 1] +
      velocity.faces[inds$corners[4,1], 4]

    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(nw.ap, ne.ap, sw.ap, se.ap, -sw.an, -se.an, -nw.as,
                             -ne.as, -nw.ae, -sw.ae, -ne.aw, -se.aw))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

    # quickly fill the sparse matrix
    A.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
                          dims = c(n.cells, n.cells))

  } else if (bc == "endless world"){

    ### i, j, and x values for interior points
    n.interior <- nrow(inds$interior)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$interior[,1], 2] -
      velocity.faces[inds$interior[,1], 1]
    F.diff.ns <- velocity.faces[inds$interior[,1], 4] -
      velocity.faces[inds$interior[,1], 3]
    F.diff.comb <- F.diff.ew + F.diff.ns

    a.w <- velocity.faces[inds$interior[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$interior[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$interior[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$interior[,1], 4]; a.n[a.n < 0] <- 0

    sparse.interior <- cbind(rep(inds$interior[,1], 5),
                             as.vector(inds$interior),
                             c(F.diff.comb + a.n + a.s + a.e + a.w, -a.n, -a.s, -a.e, -a.w))


    ### i, j, and x values for west boundary points
    n.west <- nrow(inds$west)

    # coefficient values
    F.diff.ns <- velocity.faces[inds$west[,1], 4] -
      velocity.faces[inds$west[,1], 3]
    F.diff.comb <- F.diff.ns + velocity.faces[inds$west[,1], 2]

    a.e <- -velocity.faces[inds$west[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$west[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- velocity.faces[inds$west[,1], 4]; a.n[a.n > 0] <- 0

    west.faces <- velocity.faces[inds$west[,1], 1]
    diag.vals <- F.diff.comb + a.n + a.s + a.e
    diag.vals[west.faces < 0] <- (diag.vals - west.faces)[west.faces < 0]

    sparse.west <- cbind(rep(inds$west[,1], 4),
                         as.vector(inds$west[,-5]),
                         c(diag.vals, -a.n, -a.s, -a.e))

    ### i, j, and x values for east boundary points
    n.east <- nrow(inds$east)

    # coefficient calculations
    F.diff.ns <- velocity.faces[inds$east[,1], 4] -
      velocity.faces[inds$east[,1], 3]
    F.diff.comb <- F.diff.ns - velocity.faces[inds$east[,1], 1]

    a.w <- velocity.faces[inds$east[,1], 1]; a.w[a.w < 0] <- 0
    a.s <- velocity.faces[inds$east[,1], 3]; a.s[a.s < 0] <- 0
    a.n <- -velocity.faces[inds$east[,1], 4]; a.n[a.n < 0] <- 0

    # if east face value is positive, change diagonal value to include it leaving system
    east.faces <- velocity.faces[inds$east[,1],2]
    diag.vals <- F.diff.comb + a.n + a.s + a.w
    diag.vals[east.faces > 0] <- (diag.vals + east.faces)[east.faces > 0]

    sparse.east <- cbind(rep(inds$east[,1], 4),
                         as.vector(inds$east[,-4]),
                         c(diag.vals, -a.n, -a.s, -a.w))

    ### i, j, and x values for north boundary points
    n.north <- nrow(inds$north)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$north[,1], 2] -
      velocity.faces[inds$north[,1], 1]
    F.diff.comb <- F.diff.ew - velocity.faces[inds$north[,1], 3]

    a.w <- velocity.faces[inds$north[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$north[,1], 2]; a.e[a.e < 0] <- 0
    a.s <- velocity.faces[inds$north[,1], 3]; a.s[a.s < 0] <- 0

    north.faces <- velocity.faces[inds$north[,1], 4]
    diag.vals <- F.diff.comb + a.s + a.e + a.w
    diag.vals[north.faces > 0] <- (diag.vals + north.faces)[north.faces > 0]

    sparse.north <- cbind(rep(inds$north[,1], 4),
                          as.vector(inds$north[,-2]),
                          c(diag.vals, -a.s, -a.e, -a.w))

    ### i, j, and x values for south boundary points
    n.south <- nrow(inds$south)

    # coefficient calculations
    F.diff.ew <- velocity.faces[inds$south[,1], 2] -
      velocity.faces[inds$south[,1], 1]
    F.diff.comb <- F.diff.ew + velocity.faces[inds$south[,1], 4]

    a.w <- velocity.faces[inds$south[,1], 1]; a.w[a.w < 0] <- 0
    a.e <- -velocity.faces[inds$south[,1], 2]; a.e[a.e < 0] <- 0
    a.n <- -velocity.faces[inds$south[,1], 4]; a.n[a.n < 0] <- 0

    south.faces <- velocity.faces[inds$south[,1], 3]
    diag.vals <- F.diff.comb + a.n + a.e + a.w
    diag.vals[south.faces < 0] <- (diag.vals - south.faces)[south.faces < 0]

    sparse.south <- cbind(rep(inds$south[,1], 4),
                          as.vector(inds$south[,-3]),
                          c(F.diff.comb + a.n + a.e + a.w, -a.n, -a.e, -a.w))

    ### i, j, and x values for corner points

    ## northwest values
    nw.ae <- max(c(0, -velocity.faces[inds$corners[1,1], 2]))
    nw.as <- max(c(0, velocity.faces[inds$corners[1,1], 3]))
    nw.ap <- nw.ae + nw.as + velocity.faces[inds$corners[1,1], 2] -
      velocity.faces[inds$corners[1,1], 3]

    # north face
    nw.nf <- velocity.faces[inds$corners[1,1], 4]
    # west face
    nw.wf <- velocity.faces[inds$corners[1,1], 1]

    if (nw.nf > 0){
      # add to diagonal
      nw.ap <- nw.ap + nw.nf
    }

    if (nw.wf < 0){
      # subtract from diagonal
      nw.ap <- nw.ap - nw.wf
    }

    ## northeast values
    ne.aw <- max(c(0, velocity.faces[inds$corners[2,1], 1]))
    ne.as <- max(c(0, velocity.faces[inds$corners[2,1], 3]))
    ne.ap <- ne.aw + ne.as - velocity.faces[inds$corners[2,1], 1] -
      velocity.faces[inds$corners[2,1], 3]

    # north face
    ne.nf <- velocity.faces[inds$corners[2,1], 4]
    # east face
    ne.ef <- velocity.faces[inds$corners[2,1], 2]

    if (ne.nf > 0){
      ne.ap <- ne.ap + ne.nf
    }

    if (ne.ef > 0){
      ne.ap <- ne.ap + ne.ef
    }


    ## southwest values
    sw.ae <- max(c(0, -velocity.faces[inds$corners[3,1], 2]))
    sw.an <- max(c(0, -velocity.faces[inds$corners[3,1], 4]))
    sw.ap <- sw.ae + sw.an + velocity.faces[inds$corners[3,1], 2] +
      velocity.faces[inds$corners[3,1], 4]

    # south face
    sw.sf <- velocity.faces[inds$corners[3,1],3]
    # west face
    sw.wf <- velocity.faces[inds$corners[3,1],1]

    if (sw.sf < 0){
      sw.ap <- sw.ap - sw.sf
    }

    if (sw.wf < 0){
      sw.ap <- sw.ap - sw.wf
    }

    ## southeast values
    se.aw <- max(c(0, velocity.faces[inds$corners[4,1], 1]))
    se.an <- max(c(0, -velocity.faces[inds$corners[4,1], 4]))
    se.ap <- se.aw + se.an - velocity.faces[inds$corners[4,1], 1] +
      velocity.faces[inds$corners[4,1], 4]

    # south face
    se.sf <- velocity.faces[inds$corners[4,1], 3]
    # east face
    se.ef <- velocity.faces[inds$corners[4,1], 2]

    if (se.sf < 0){
      se.ap <- se.ap - se.sf
    }

    if (se.ef > 0){
      se.ap <- se.ap + se.ef
    }

    # combine all values together
    sparse.corner <- cbind(c(inds$corners[,1], inds$corners[3:4,1], inds$corners[1:2,1],
                             inds$corners[c(1,3),1], inds$corners[c(2,4),1]),
                           c(inds$corners[,1], inds$corners[3:4,2], inds$corners[1:2,3],
                             inds$corners[c(1,3),4], inds$corners[c(2,4),5]),
                           c(nw.ap, ne.ap, sw.ap, se.ap, -sw.an, -se.an, -nw.as,
                             -ne.as, -nw.ae, -sw.ae, -ne.aw, -se.aw))

    # combine i, j, and x matrices into single matrix
    sparse.ijx <- rbind(sparse.interior,
                        sparse.west,
                        sparse.east,
                        sparse.north,
                        sparse.south,
                        sparse.corner)

    # quickly fill the sparse matrix
    A.mat <- sparseMatrix(sparse.ijx[,1], sparse.ijx[,2], x = sparse.ijx[,3],
                          dims = c(n.cells, n.cells))
  }

  return(A.mat)
}

createC <- function(big, small, dims, inds, bound.cond = "endless world"){
  # Create matrix fo discretized advection process, due to wind.
  # Input:
  #   big: raster with wind data, large spatial extent
  #   small: raster with wind data, small spatial extent
  #   dims: dimension of the raster
  #   inds: list of indices, created with FVMInds
  #   bound.cond: boundary condition (periodic or endless world)
  # Output:
  #   Sparse matrix representing advection due to wind, as specified by
  #     wind velocities found in 'small'.

  big.uwind <- big$uwind
  big.vwind <- big$vwind

  small.uwind <- small$uwind
  small.vwind <- small$vwind

  u.faces <- interpWind(big.uwind, small.uwind, "u", with.bound = FALSE)
  v.faces <- interpWind(big.vwind, small.vwind, "v", with.bound = FALSE)
  velocity <- cbind(u.faces, v.faces)

  # wind vmatrix
  C <- windConvFVM(dims, velocity, inds, bc = bound.cond, bound = FALSE)
  return(C)
}


########################################################################
### B. Outcome Analysis Functions
########################################################################

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

tuningMatFull <- function() {
  # Creates an initial covariance matrix (32x32) for the M-H proposal distribution.
  
  matrix(
    data = c(
      4.29e-01, -1.18e-03,  1.39e-03, -2.15e-05, -3.74e-03, 
      -4.48e-03, -7.09e-03, -3.16e-03, -4.15e-03, -1.07e-03,  
      3.35e-04, -2.14e-03, -7.90e-04,  6.08e-04, -2.06e-04, 
      -2.59e-02,  8.83e-04, -5.48e-04,  4.32e-07,  9.05e-07, 
      -6.25e-05, -3.89e-03, -4.35e-05, -3.34e-03, -1.09e-04,  
      7.04e-03, -2.66e-07, -2.25e-02,  1.83e-02, -2.41e-04,  
      4.70e-03,  8.70e-04,      # row one
      -1.18e-03,  3.80e-03, -8.60e-04,  7.65e-07,  2.32e-04,  
      3.22e-04,  3.75e-04, -4.36e-04, -5.76e-04,  8.32e-06,  
      7.03e-04,  3.86e-04,  1.18e-04, -6.26e-05, -1.11e-05,  
      3.83e-04,  1.20e-05,  1.73e-05,  3.05e-08, -1.01e-08,
      -5.15e-06,  1.90e-04, -3.84e-06,  3.90e-06, -2.05e-05,  
      5.91e-05, -3.67e-07, -6.80e-05, -2.95e-04,  6.31e-05, 
      3.75e-05, -1.26e-04,      # row two
      1.39e-03, -8.60e-04,  5.65e-03,  3.44e-05,  3.91e-04, 
      9.64e-05,  3.99e-05,  3.18e-04, -2.96e-04, -2.48e-04,  
      1.76e-03, -4.06e-05,  2.73e-04, -8.82e-06,  9.60e-06, 
      -1.09e-03, -5.85e-05,  7.40e-05, -7.80e-09,  3.24e-07,
      2.13e-06,  1.35e-04, -3.18e-06, -5.87e-05,  8.49e-05,  
      1.33e-04,  1.37e-07, -2.00e-04, -2.14e-04,  2.29e-05,  
      6.21e-05, -9.96e-05,      # row three
      -2.15e-05,  7.65e-07,  3.44e-05,  5.49e-07,  2.69e-07,  
      1.02e-06,  3.96e-06,  4.57e-06,  1.60e-07, -2.79e-07,  
      1.90e-05, -2.93e-06,  1.75e-06, -3.46e-07,  3.00e-07, 
      -6.81e-06, -7.85e-07,  7.34e-07,  5.15e-10,  2.97e-09,
      2.50e-08,  1.39e-07, -1.32e-08,  4.49e-07, -8.04e-07,  
      2.94e-06,  4.84e-09, -1.78e-06,  8.07e-06, -6.72e-07, 
      -3.14e-06, -6.32e-07,     # row four
      -3.74e-03,  2.32e-04,  3.91e-04,  2.69e-07,  4.72e-03,
      4.46e-03,  4.19e-03,  1.39e-04,  9.25e-05,  1.00e-04,
      1.27e-04, -2.06e-04,  2.52e-04,  6.06e-05,  1.68e-05,
      -2.38e-03,  1.15e-06,  1.75e-05,  3.94e-08,  4.14e-08,
      -5.55e-06,  1.73e-04,  5.79e-06, -1.61e-05,  1.35e-05,
      1.75e-05, -3.22e-07, -1.25e-04,  8.28e-06,  9.87e-05,
      1.99e-04, -2.87e-04,      # row five
      -4.48e-03,  3.22e-04,  9.64e-05,  1.02e-06,  4.46e-03,
      4.88e-03,  4.14e-03,  1.63e-04,  1.31e-04,  2.09e-04,
      1.68e-04, -1.93e-04,  1.99e-04, -3.68e-05,  2.14e-05,
      -1.99e-03, -9.74e-06,  2.19e-05, -2.40e-08,  6.46e-08,
      -2.54e-06,  1.29e-04,  5.74e-06, -2.31e-05,  9.85e-06,
      -3.11e-05, -1.65e-07, -1.44e-04, -1.13e-04, -1.43e-05,
      9.71e-05, -7.96e-05,      # row six
      -7.09e-03,  3.75e-04,  3.99e-05,  3.96e-06,  4.19e-03,
      4.14e-03,  4.73e-03,  1.01e-03,  2.52e-04,  1.10e-04,
      3.44e-04,  2.17e-04,  1.61e-04, -4.39e-05,  1.20e-05,
      -1.29e-03, -1.48e-05, -1.73e-05,  9.82e-08, -1.66e-07,
      5.76e-06,  3.98e-04,  5.95e-06,  5.26e-05, -3.19e-05,
      2.27e-05, -1.41e-07, -2.95e-06,  1.20e-04,  9.93e-05,
      8.31e-05, -2.88e-04,     # row seven
      -3.16e-03, -4.36e-04,  3.18e-04,  4.57e-06,  1.39e-04,
      1.63e-04,  1.01e-03,  3.05e-03,  7.92e-04, -1.46e-04,
      -6.11e-05, -9.71e-04, -1.67e-04, -5.92e-05, -1.02e-05,
      5.89e-04, -6.74e-07, -8.80e-06,  1.30e-07, -5.17e-08,
      -2.60e-06,  3.41e-04, -8.28e-07,  5.93e-05,  5.07e-06,
      1.57e-04, -3.93e-08,  1.09e-04,  1.18e-04,  1.65e-04,
      3.94e-06, -3.04e-04,    # row eight
      -4.15e-03, -5.76e-04, -2.96e-04,  1.60e-07,  9.25e-05,
      1.31e-04,  2.52e-04,  7.92e-04,  2.93e-03,  5.49e-04,
      -6.18e-04, -1.80e-05, -1.46e-05, -7.80e-06, -1.53e-05,
      7.10e-04, -9.07e-06, -2.10e-05,  4.88e-08,  2.37e-07,
      -3.39e-06, -2.41e-05,  3.76e-06,  2.76e-05, -3.02e-05,
      3.36e-05,  2.88e-07, -1.62e-04,  2.27e-04, -1.46e-04,
      -1.05e-04,  2.03e-04,    # row nine
      -1.07e-03,  8.32e-06, -2.48e-04, -2.79e-07,  1.00e-04,
      2.09e-04,  1.10e-04, -1.46e-04,  5.49e-04,  2.99e-04,
      -2.08e-04, -1.73e-04,  1.72e-04, -1.95e-05, -3.19e-06,
      4.34e-04, -1.01e-05,  9.49e-06, -1.74e-08,  3.73e-08,
      1.30e-06, -7.13e-05,  9.15e-07, -1.27e-05, -1.85e-05,
      -8.13e-06,  1.51e-07, -1.31e-04,  2.39e-04, -7.07e-05,
      -3.66e-05,  1.16e-04,    # row ten
      3.35e-04,  7.03e-04,  1.76e-03,  1.90e-05,  1.27e-04,
      1.68e-04,  3.44e-04, -6.11e-05, -6.18e-04, -2.08e-04,
      3.20e-03, -8.96e-05, -7.18e-04, -7.08e-05,  4.18e-06,
      -7.38e-04, -3.22e-05,  3.61e-05, -1.94e-08,  8.30e-08,
      2.54e-06,  3.74e-04,  2.75e-06, -2.97e-06,  4.76e-06,
      1.13e-04, -2.81e-07, -4.60e-05, -3.95e-04,  4.05e-05,
      1.23e-04, -1.15e-04,    # row eleven
      -2.14e-03,  3.86e-04, -4.06e-05, -2.93e-06, -2.06e-04,
      -1.93e-04,  2.17e-04, -9.71e-04, -1.80e-05, -1.73e-04,
      -8.96e-05,  3.51e-03,  1.06e-04, -5.91e-05,  8.53e-06,
      8.76e-04,  3.98e-05, -5.86e-05, -1.44e-08, -4.13e-07,
      1.22e-05,  1.67e-04,  1.81e-06,  6.87e-05, -3.20e-05,
      -1.87e-04,  3.39e-07,  1.19e-04, -2.46e-04, -1.15e-04,
      -1.63e-04, 1.87e-04,     # row twelve
      -7.90e-04,  1.18e-04,  2.73e-04,  1.75e-06,  2.52e-04,
      1.99e-04,  1.61e-04, -1.67e-04, -1.46e-05,  1.72e-04,
      -7.18e-04,  1.06e-04,  8.82e-04, -2.88e-05, -1.13e-05,
      1.62e-04,  3.99e-06,  1.07e-05,  7.43e-09, -1.07e-07,
      3.07e-06, -1.44e-04,  7.93e-07, -1.98e-05, -3.20e-07,
      -2.30e-05,  3.18e-07, -1.00e-04,  4.35e-04, -3.73e-05,
      -1.28e-04,  5.25e-05,   # row thirteen
      6.08e-04, -6.26e-05, -8.82e-06, -3.46e-07,  6.06e-05,
      -3.68e-05, -4.39e-05, -5.92e-05, -7.80e-06, -1.95e-05,
      -7.08e-05, -5.91e-05, -2.88e-05,  2.53e-04, -3.30e-05,
      3.59e-05,  5.56e-06, -5.51e-06,  1.68e-09, -5.22e-08,
      1.24e-06, -3.62e-05, -2.31e-07, -2.86e-05,  3.04e-05,
      -3.10e-05, -4.36e-08,  2.53e-07,  1.24e-04,  7.17e-05,
      9.75e-05, -9.30e-05,    # row fourteen
      -2.06e-04, -1.11e-05,  9.60e-06,  3.00e-07,  1.68e-05,
      2.14e-05,  1.20e-05, -1.02e-05, -1.53e-05, -3.19e-06,
      4.18e-06,  8.53e-06, -1.13e-05, -3.30e-05,  1.02e-05,
      8.20e-05, -7.16e-06,  5.85e-06, -3.96e-09,  4.09e-08,
      -3.07e-08, -2.69e-05,  2.10e-07,  3.46e-06, -4.55e-06,
      2.02e-06,  5.56e-08,  3.36e-06,  3.37e-05, -1.09e-05,
      -1.83e-05,  1.42e-05,   # row fifteen
      -2.59e-02,  3.83e-04, -1.09e-03, -6.81e-06, -2.38e-03,
      -1.99e-03, -1.29e-03,  5.89e-04,  7.10e-04,  4.34e-04,
      -7.38e-04,  8.76e-04,  1.62e-04,  3.59e-05,  8.20e-05,
      2.39e-02, -9.66e-05,  3.33e-05, -5.15e-07,  1.68e-06,
      -6.09e-07,  4.41e-04, -1.65e-06, -1.19e-04,  5.72e-04,
      5.06e-04,  1.96e-06,  8.27e-04, -2.20e-03, -4.51e-04,
      -6.09e-04,  9.07e-04,    # row sixteen
      8.83e-04,  1.20e-05, -5.85e-05, -7.85e-07,  1.15e-06,
      -9.74e-06, -1.48e-05, -6.74e-07, -9.07e-06, -1.01e-05,
      -3.22e-05,  3.98e-05,  3.99e-06,  5.56e-06, -7.16e-06,
      -9.66e-05,  1.10e-04, -7.32e-05,  3.82e-08, -5.58e-07,
      -4.39e-06,  1.03e-05, -4.74e-07,  7.52e-06, -9.29e-06,
      -1.18e-04, -4.92e-07,  6.32e-05, -5.33e-05,  2.26e-05,
      -3.48e-05,  5.13e-06,    # row seventeen
      -5.48e-04,  1.73e-05,  7.40e-05,  7.34e-07,  1.75e-05,
      2.19e-05, -1.73e-05, -8.80e-06, -2.10e-05,  9.49e-06,
      3.61e-05, -5.86e-05,  1.07e-05, -5.51e-06,  5.85e-06,
      3.33e-05, -7.32e-05,  7.05e-05, -1.39e-08,  3.12e-07,
      2.02e-06,  4.45e-05,  2.01e-07,  1.55e-05, -2.55e-05,
      5.95e-05,  2.74e-07, -7.95e-05,  1.34e-04, -6.65e-05,
      4.83e-06,  5.04e-05,    # row eighteen
      4.32e-07,  3.05e-08, -7.80e-09,  5.15e-10,  3.94e-08,
      -2.40e-08,  9.82e-08,  1.30e-07,  4.88e-08, -1.74e-08,
      -1.94e-08, -1.44e-08,  7.43e-09,  1.68e-09, -3.96e-09,
      -5.15e-07,  3.82e-08, -1.39e-08,  4.07e-10, -2.40e-11,
      -1.99e-08, -1.12e-07,  3.43e-10,  7.56e-08, -1.02e-07,
      1.46e-07,  1.64e-10,  4.20e-09,  7.98e-07,  5.86e-08,
      -1.84e-07, -1.43e-07,   # row nineteen
      9.05e-07, -1.01e-08,  3.24e-07,  2.97e-09,  4.14e-08,
      6.46e-08, -1.66e-07, -5.17e-08,  2.37e-07,  3.73e-08,
      8.30e-08, -4.13e-07, -1.07e-07, -5.22e-08,  4.09e-08,
      1.68e-06, -5.58e-07,  3.12e-07, -2.40e-11,  5.09e-09,
      -3.92e-08, -1.58e-07,  3.06e-09, -1.65e-07,  1.96e-07,
      1.23e-06,  3.54e-09, -4.80e-07,  1.90e-06, -5.90e-07,
      -7.96e-07,  5.54e-07,   # row twenty
      -6.25e-05, -5.15e-06,  2.13e-06,  2.50e-08, -5.55e-06,
      -2.54e-06,  5.76e-06, -2.60e-06, -3.39e-06,  1.30e-06,
      2.54e-06,  1.22e-05,  3.07e-06,  1.24e-06, -3.07e-08,
      -6.09e-07, -4.39e-06,  2.02e-06, -1.99e-08, -3.92e-08,
      2.88e-06, -5.08e-06, -2.05e-08, -4.35e-06,  5.08e-06,
      -1.63e-05,  1.51e-08, -2.05e-06, -6.43e-05,  4.51e-07,
      2.48e-05,  5.45e-06,    # row twenty one
      -3.89e-03,  1.90e-04,  1.35e-04,  1.39e-07,  1.73e-04,
      1.29e-04,  3.98e-04,  3.41e-04, -2.41e-05, -7.13e-05,
      3.74e-04,  1.67e-04, -1.44e-04, -3.62e-05, -2.69e-05,
      4.41e-04,  1.03e-05,  4.45e-05, -1.12e-07, -1.58e-07,
      -5.08e-06,  2.10e-03,  2.48e-06,  1.45e-04, -1.11e-04,
      8.22e-05, -1.88e-06,  8.10e-05, -2.69e-03, -1.80e-05,
      8.11e-05, -3.42e-05,    # row twenty two
      -4.35e-05, -3.84e-06, -3.18e-06, -1.32e-08,  5.79e-06,
      5.74e-06,  5.95e-06, -8.28e-07,  3.76e-06,  9.15e-07,
      2.75e-06,  1.81e-06,  7.93e-07, -2.31e-07,  2.10e-07,
      -1.65e-06, -4.74e-07,  2.01e-07,  3.43e-10,  3.06e-09,
      -2.05e-08,  2.48e-06,  5.91e-07,  1.94e-07,  2.51e-08,
      -5.52e-08,  5.45e-09,  5.31e-07,  5.08e-06, -4.11e-06,
      -4.66e-06,  5.00e-06,   # row twenty three
      -3.34e-03,  3.90e-06, -5.87e-05,  4.49e-07, -1.61e-05,
      -2.31e-05,  5.26e-05,  5.93e-05,  2.76e-05, -1.27e-05,
      -2.97e-06,  6.87e-05, -1.98e-05, -2.86e-05,  3.46e-06,
      -1.19e-04,  7.52e-06,  1.55e-05,  7.56e-08, -1.65e-07,
      -4.35e-06,  1.45e-04,  1.94e-07,  1.97e-04, -1.73e-04,
      -4.93e-05, -1.09e-07,  1.03e-04, -1.38e-05, -5.76e-05,
      -9.46e-05,  4.21e-05,   # row twenty four
      -1.09e-04, -2.05e-05,  8.49e-05, -8.04e-07,  1.35e-05,
      9.85e-06, -3.19e-05,  5.07e-06, -3.02e-05, -1.85e-05,
      4.76e-06, -3.20e-05, -3.20e-07,  3.04e-05, -4.55e-06,
      5.72e-04, -9.29e-06, -2.55e-05, -1.02e-07,  1.96e-07,
      5.08e-06, -1.11e-04,  2.51e-08, -1.73e-04,  2.35e-04,
      5.52e-05,  2.28e-07,  9.36e-05, -2.89e-04,  9.53e-05,
      5.45e-05, -8.57e-05,    # row twenty five
      7.04e-03,  5.91e-05,  1.33e-04,  2.94e-06,  1.75e-05,
      -3.11e-05,  2.27e-05,  1.57e-04,  3.36e-05, -8.13e-06,
      1.13e-04, -1.87e-04, -2.30e-05, -3.10e-05,  2.02e-06,
      5.06e-04, -1.18e-04,  5.95e-05,  1.46e-07,  1.23e-06,
      -1.63e-05,  8.22e-05, -5.52e-08, -4.93e-05,  5.52e-05,
      2.38e-03,  2.13e-07, -6.00e-04, -4.37e-04,  4.45e-04,
      1.84e-04, -8.05e-04,    # row twenty six
      -2.66e-07, -3.67e-07,  1.37e-07,  4.84e-09, -3.22e-07,
      -1.65e-07, -1.41e-07, -3.93e-08,  2.88e-07,  1.51e-07,
      -2.81e-07,  3.39e-07,  3.18e-07, -4.36e-08,  5.56e-08,
      1.96e-06, -4.92e-07,  2.74e-07,  1.64e-10,  3.54e-09,
      1.51e-08, -1.88e-06,  5.45e-09, -1.09e-07,  2.28e-07,
      2.13e-07,  9.61e-09, -5.85e-07,  6.97e-06, -1.98e-06,
      -2.26e-06,  2.57e-06,   # row twenty seven
      -2.25e-02, -6.80e-05, -2.00e-04, -1.78e-06, -1.25e-04,
      -1.44e-04, -2.95e-06,  1.09e-04, -1.62e-04, -1.31e-04,
      -4.60e-05,  1.19e-04, -1.00e-04,  2.53e-07,  3.36e-06,
      8.27e-04,  6.32e-05, -7.95e-05,  4.20e-09, -4.80e-07,
      -2.05e-06,  8.10e-05,  5.31e-07,  1.03e-04,  9.36e-05,
      -6.00e-04, -5.85e-07,  1.57e-03, -1.26e-03,  1.53e-04,
      -2.68e-04, -1.79e-04,   # row twenty eight
      1.83e-02, -2.95e-04, -2.14e-04,  8.07e-06,  8.28e-06,
      -1.13e-04,  1.20e-04,  1.18e-04,  2.27e-04,  2.39e-04,
      -3.95e-04, -2.46e-04,  4.35e-04,  1.24e-04,  3.37e-05,
      -2.20e-03, -5.33e-05,  1.34e-04,  7.98e-07,  1.90e-06,
      -6.43e-05, -2.69e-03,  5.08e-06, -1.38e-05, -2.89e-04,
      -4.37e-04,  6.97e-06, -1.26e-03,  1.66e-02, -2.31e-03,
      -3.05e-03,  2.72e-03,   # row twenty nine
      -2.41e-04,  6.31e-05,  2.29e-05, -6.72e-07,  9.87e-05,
      -1.43e-05,  9.93e-05,  1.65e-04, -1.46e-04, -7.07e-05,
      4.05e-05, -1.15e-04, -3.73e-05,  7.17e-05, -1.09e-05,
      -4.51e-04,  2.26e-05, -6.65e-05,  5.86e-08, -5.90e-07,
      4.51e-07, -1.80e-05, -4.11e-06, -5.76e-05,  9.53e-05,
      4.45e-04, -1.98e-06,  1.53e-04, -2.31e-03,  2.74e-03,
      1.99e-03, -3.69e-03,   # row thirty
      4.70e-03,  3.75e-05,  6.21e-05, -3.14e-06,  1.99e-04,
      9.71e-05,  8.31e-05,  3.94e-06, -1.05e-04, -3.66e-05,
      1.23e-04, -1.63e-04, -1.28e-04,  9.75e-05, -1.83e-05,
      -6.09e-04, -3.48e-05,  4.83e-06, -1.84e-07, -7.96e-07,
      2.48e-05,  8.11e-05, -4.66e-06, -9.46e-05,  5.45e-05,
      1.84e-04, -2.26e-06, -2.68e-04, -3.05e-03,  1.99e-03,
      2.60e-03, -2.66e-03,   # row thirty one
      8.70e-04, -1.26e-04, -9.96e-05, -6.32e-07, -2.87e-04,
      -7.96e-05, -2.88e-04, -3.04e-04,  2.03e-04,  1.16e-04,
      -1.15e-04,  1.87e-04,  5.25e-05, -9.30e-05,  1.42e-05,
      9.07e-04,  5.13e-06,  5.04e-05, -1.43e-07,  5.54e-07,
      5.45e-06, -3.42e-05,  5.00e-06,  4.21e-05, -8.57e-05,
      -8.05e-04,  2.57e-06, -1.79e-04,  2.72e-03, -3.69e-03,
      -2.66e-03,  5.11e-03   # row thirty two
    ),
    nrow = 32, ncol = 32, byrow = TRUE
  )
}

IndirectEffectCalc <- function(mu){
  # Estimates IE(z,g) from mu(z,g).
  # Input:
  #   mu: posterior estimates of mu(z,g).
  # Output:
  #   Posterior samples of IE(z,g).
  
  t(apply(mu, 1, FUN = function(vec){vec - vec[1]})) 
}

poisEstimates <- function(
    fit, n_post, X, degree_idx, dist_idx, up_heat_idx, heat_idx,
    degree_cuts = c(0.16, 0.4),
    dist_cuts = c(50, 100, 200),
    heat_cut = 15.5
){
  # Estimates the direct and indirect causal effects using output from a 
  #   Poisson regression outcome model. This includes estimates stratified
  #   by levels of degree, distance to key-associated facilities, and 
  #   key-associated heat input.
  # Input:
  #   fit: Poisson regression output from 'pois_reg_cpp'.
  #   n_post: Number of samples from posterior to include in casual effect estimate.
  #   X: Design matrix.
  #   degree_idx: Weighted degree (covariate matrix column) index.
  #   dist_idx: Distance to key-associated facility index.
  #   up_heat_idx: Upwind heat input index
  #   heat_idx: Key-associated heat input index
  #   degree_cuts: Cut thresholds for degree.
  #   dist_cuts: Cut thresholds for distance.
  #   heat_cut: Cut thresholds for (log-) heat index.
  # Output: 
  #   List containing posterior estimates of mu(z,g), DE(g), and IE(z,g).
  
  # create G sequence
  g_seq <- seq(from = 0.25, to = 0.9, by = 0.01)
  n_g <- length(g_seq)
  
  # stratify by degree
  low_deg <- which(X[, degree_idx] < degree_cuts[1])
  med_deg <- which(X[, degree_idx] >= degree_cuts[1] & X[, degree_idx] < degree_cuts[2])
  high_deg <- which(X[, degree_idx] >= degree_cuts[2])
  
  # stratify by degree and upwind heat input
  ld_lh <- which((X[, degree_idx] < degree_cuts[1]) & (X[, up_heat_idx] < heat_cut))
  ld_hh <- which((X[, degree_idx] < degree_cuts[1]) & (X[, up_heat_idx] >= heat_cut))
  md_lh <- which((X[, degree_idx] >= degree_cuts[1] & X[, degree_idx] < degree_cuts[2]) & (X[, up_heat_idx] < heat_cut))
  md_hh <- which((X[, degree_idx] >= degree_cuts[1] & X[, degree_idx] < degree_cuts[2]) & (X[, up_heat_idx] >= heat_cut))
  hd_lh <- which((X[, degree_idx] >= degree_cuts[2]) & (X[, up_heat_idx] < heat_cut))
  hd_hh <- which((X[, degree_idx] >= degree_cuts[2]) & (X[, up_heat_idx] >= heat_cut))
  
  # stratify by distance to key associated facility
  dist1 <- which(X[, dist_idx] < dist_cuts[1])
  dist2 <- which(X[, dist_idx] >= dist_cuts[1] & X[, dist_idx] < dist_cuts[2])
  dist3 <- which(X[, dist_idx] >= dist_cuts[2] & X[, dist_idx] < dist_cuts[3])
  dist4 <- which(X[, dist_idx] >= dist_cuts[3])
  
  # stratify by distance and heat input
  d1_lh <- which((X[, dist_idx] < dist_cuts[1]) & (X[, heat_idx] < heat_cut))
  d1_hh <- which((X[, dist_idx] < dist_cuts[1]) & (X[, heat_idx] >= heat_cut))
  d2_lh <- which((X[, dist_idx] >= dist_cuts[1] & X[, dist_idx] < dist_cuts[2]) & (X[, heat_idx] < heat_cut))
  d2_hh <- which((X[, dist_idx] >= dist_cuts[1] & X[, dist_idx] < dist_cuts[2]) & (X[, heat_idx] >= heat_cut))
  d3_lh <- which((X[, dist_idx] >= dist_cuts[2] & X[, dist_idx] < dist_cuts[3]) & (X[, heat_idx] < heat_cut))
  d3_hh <- which((X[, dist_idx] >= dist_cuts[2] & X[, dist_idx] < dist_cuts[3]) & (X[, heat_idx] >= heat_cut))
  d4_lh <- which((X[, dist_idx] >= dist_cuts[3]) & (X[, heat_idx] < heat_cut))
  d4_hh <- which((X[, dist_idx] >= dist_cuts[3]) & (X[, heat_idx] >= heat_cut))
  
  # mu0 list
  mu0 <- list(
    all = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ldeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    mdeg = matrix(NA_real_, nrow = n_post, ncol = n_g), 
    hdeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_hh = matrix(NA_real_, nrow = n_post, ncol = n_g)
  )
  
  # m1 list
  mu1 <- list(
    all = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ldeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    mdeg = matrix(NA_real_, nrow = n_post, ncol = n_g), 
    hdeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_hh = matrix(NA_real_, nrow = n_post, ncol = n_g)
  )
  
  ### estimate mu(z,g) for different strata
  
  # subsample parameter values
  mcmc.samp <- sample(nrow(fit[[1]]), n_post)
  betas <- fit[[1]][mcmc.samp,]
  
  # z = 0
  z <- 0
  for (k in 1:n_g){
    
    # set up X matrix
    g_val <- g_seq[k]
    X_m <- cbind(rep(1, nrow(X)), X, z, g_val, z * g_val)
    colnames(X_m)[(ncol(X_m) - 2):ncol(X_m)] <- c("Z", "G", "Z.G")
    
    # predict rates
    f_hat <- exp(X_m %*% t(betas))
    
    # mu(z,g)
    mu0$all[, k] <- colMeans(f_hat)
    if (length(low_deg) > 1) mu0$ldeg[, k] <- colMeans(f_hat[low_deg, ])
    if (length(med_deg) > 1) mu0$mdeg[, k] <- colMeans(f_hat[med_deg, ])
    if (length(high_deg) > 1) mu0$hdeg[, k] <- colMeans(f_hat[high_deg, ])
    if (length(ld_lh) > 1) mu0$ld_lh[, k] <- colMeans(f_hat[ld_lh, ])
    if (length(ld_hh) > 1) mu0$ld_hh[, k] <- colMeans(f_hat[ld_hh, ])
    if (length(md_lh) > 1) mu0$md_lh[, k] <- colMeans(f_hat[md_lh, ])
    if (length(md_hh) > 1) mu0$md_hh[, k] <- colMeans(f_hat[md_hh, ])
    if (length(hd_lh) > 1) mu0$hd_lh[, k] <- colMeans(f_hat[hd_lh, ])
    if (length(hd_hh) > 1) mu0$hd_hh[, k] <- colMeans(f_hat[hd_hh, ])
    if (length(dist1) > 1) mu0$d1[, k] <- colMeans(f_hat[dist1, ])
    if (length(dist2) > 1) mu0$d2[, k] <- colMeans(f_hat[dist2, ])
    if (length(dist3) > 1) mu0$d3[, k] <- colMeans(f_hat[dist3, ])
    if (length(dist4) > 1) mu0$d4[, k] <- colMeans(f_hat[dist4, ])
    if (length(d1_lh) > 1) mu0$d1_lh[, k] <- colMeans(f_hat[d1_lh, ])
    if (length(d1_hh) > 1) mu0$d1_hh[, k] <- colMeans(f_hat[d1_hh, ])
    if (length(d2_lh) > 1) mu0$d2_lh[, k] <- colMeans(f_hat[d2_lh, ])
    if (length(d2_hh) > 1) mu0$d2_hh[, k] <- colMeans(f_hat[d2_hh, ])
    if (length(d3_lh) > 1) mu0$d3_lh[, k] <- colMeans(f_hat[d3_lh, ])
    if (length(d3_hh) > 1) mu0$d3_hh[, k] <- colMeans(f_hat[d3_hh, ])
    if (length(d4_lh) > 1) mu0$d4_lh[, k] <- colMeans(f_hat[d4_lh, ])
    if (length(d4_hh) > 1) mu0$d4_hh[, k] <- colMeans(f_hat[d4_hh, ])
  }
  
  # z = 1 
  z <- 1
  for (k in 1:n_g){
    
    # set up X matrix
    g_val <- g_seq[k]
    X_m <- cbind(rep(1, nrow(X)), X, z, g_val, z * g_val)
    colnames(X_m)[(ncol(X_m) - 2):ncol(X_m)] <- c("Z", "G", "Z.G")
    
    # predict rates
    f_hat <- exp(X_m %*% t(betas))
    
    # mu(z,g)
    mu1$all[, k] <- colMeans(f_hat)
    if (length(low_deg) > 1) mu1$ldeg[, k] <- colMeans(f_hat[low_deg, ])
    if (length(med_deg) > 1) mu1$mdeg[, k] <- colMeans(f_hat[med_deg, ])
    if (length(high_deg) > 1) mu1$hdeg[, k] <- colMeans(f_hat[high_deg, ])
    if (length(ld_lh) > 1) mu1$ld_lh[, k] <- colMeans(f_hat[ld_lh, ])
    if (length(ld_hh) > 1) mu1$ld_hh[, k] <- colMeans(f_hat[ld_hh, ])
    if (length(md_lh) > 1) mu1$md_lh[, k] <- colMeans(f_hat[md_lh, ])
    if (length(md_hh) > 1) mu1$md_hh[, k] <- colMeans(f_hat[md_hh, ])
    if (length(hd_lh) > 1) mu1$hd_lh[, k] <- colMeans(f_hat[hd_lh, ])
    if (length(hd_hh) > 1) mu1$hd_hh[, k] <- colMeans(f_hat[hd_hh, ])
    if (length(dist1) > 1) mu1$d1[, k] <- colMeans(f_hat[dist1, ])
    if (length(dist2) > 1) mu1$d2[, k] <- colMeans(f_hat[dist2, ])
    if (length(dist3) > 1) mu1$d3[, k] <- colMeans(f_hat[dist3, ])
    if (length(dist4) > 1) mu1$d4[, k] <- colMeans(f_hat[dist4, ])
    if (length(d1_lh) > 1) mu1$d1_lh[, k] <- colMeans(f_hat[d1_lh, ])
    if (length(d1_hh) > 1) mu1$d1_hh[, k] <- colMeans(f_hat[d1_hh, ])
    if (length(d2_lh) > 1) mu1$d2_lh[, k] <- colMeans(f_hat[d2_lh, ])
    if (length(d2_hh) > 1) mu1$d2_hh[, k] <- colMeans(f_hat[d2_hh, ])
    if (length(d3_lh) > 1) mu1$d3_lh[, k] <- colMeans(f_hat[d3_lh, ])
    if (length(d3_hh) > 1) mu1$d3_hh[, k] <- colMeans(f_hat[d3_hh, ])
    if (length(d4_lh) > 1) mu1$d4_lh[, k] <- colMeans(f_hat[d4_lh, ])
    if (length(d4_hh) > 1) mu1$d4_hh[, k] <- colMeans(f_hat[d4_hh, ])
  }
  
  # estimate DE(g) = mu(1, g) - mu(0, g)
  DE <- list(
    all = mu1$all - mu0$all,
    ldeg = mu1$ldeg - mu0$ldeg,
    mdeg = mu1$mdeg - mu0$mdeg,
    hdeg = mu1$hdeg - mu0$hdeg,
    ld_lh = mu1$ld_lh - mu0$ld_lh,
    ld_hh = mu1$ld_hh - mu0$ld_hh,
    md_lh = mu1$md_lh - mu0$md_lh,
    md_hh = mu1$md_hh - mu0$md_hh,
    hd_lh = mu1$hd_lh - mu0$hd_lh,
    hd_hh = mu1$hd_hh - mu0$hd_hh,
    d1 = mu1$d1 - mu0$d1,
    d2 = mu1$d2 - mu0$d2,
    d3 = mu1$d3 - mu0$d3,
    d4 = mu1$d4 - mu0$d4,
    d1_lh = mu1$d1_lh - mu0$d1_lh,
    d1_hh = mu1$d1_hh - mu0$d1_hh,
    d2_lh = mu1$d2_lh - mu0$d2_lh,
    d2_hh = mu1$d2_hh - mu0$d2_hh,
    d3_lh = mu1$d3_lh - mu0$d3_lh,
    d3_hh = mu1$d3_hh - mu0$d3_hh,
    d4_lh = mu1$d4_lh - mu0$d4_lh,
    d4_hh = mu1$d4_hh - mu0$d4_hh
  )
  
  # IE(0, g) = mu(0, g) - mu(0, g* = 0.25)
  IE0 <- list(
    all = IndirectEffectCalc(mu0$all),
    ldeg = IndirectEffectCalc(mu0$ldeg), 
    mdeg = IndirectEffectCalc(mu0$mdeg),
    hdeg = IndirectEffectCalc(mu0$hdeg),
    ld_lh = IndirectEffectCalc(mu0$ld_lh),
    ld_hh = IndirectEffectCalc(mu0$ld_hh),
    md_lh = IndirectEffectCalc(mu0$md_lh),
    md_hh = IndirectEffectCalc(mu0$md_hh),
    hd_lh = IndirectEffectCalc(mu0$hd_lh),
    hd_hh = IndirectEffectCalc(mu0$hd_hh),
    d1 = IndirectEffectCalc(mu0$d1),
    d2 = IndirectEffectCalc(mu0$d2),
    d3 = IndirectEffectCalc(mu0$d3),
    d4 = IndirectEffectCalc(mu0$d4),
    d1_lh = IndirectEffectCalc(mu0$d1_lh),
    d1_hh = IndirectEffectCalc(mu0$d1_hh),
    d2_lh = IndirectEffectCalc(mu0$d2_lh),
    d2_hh = IndirectEffectCalc(mu0$d2_hh),
    d3_lh = IndirectEffectCalc(mu0$d3_lh),
    d3_hh = IndirectEffectCalc(mu0$d3_hh),
    d4_lh = IndirectEffectCalc(mu0$d4_lh),
    d4_hh = IndirectEffectCalc(mu0$d4_hh)
  )
  
  # IE(1, g) = mu(1, g) - mu(1, g* = 0.25)
  IE1 <- list(
    all = IndirectEffectCalc(mu1$all),
    ldeg = IndirectEffectCalc(mu1$ldeg), 
    mdeg = IndirectEffectCalc(mu1$mdeg),
    hdeg = IndirectEffectCalc(mu1$hdeg),
    ld_lh = IndirectEffectCalc(mu1$ld_lh),
    ld_hh = IndirectEffectCalc(mu1$ld_hh),
    md_lh = IndirectEffectCalc(mu1$md_lh),
    md_hh = IndirectEffectCalc(mu1$md_hh),
    hd_lh = IndirectEffectCalc(mu1$hd_lh),
    hd_hh = IndirectEffectCalc(mu1$hd_hh),
    d1 = IndirectEffectCalc(mu1$d1),
    d2 = IndirectEffectCalc(mu1$d2),
    d3 = IndirectEffectCalc(mu1$d3),
    d4 = IndirectEffectCalc(mu1$d4),
    d1_lh = IndirectEffectCalc(mu1$d1_lh),
    d1_hh = IndirectEffectCalc(mu1$d1_hh),
    d2_lh = IndirectEffectCalc(mu1$d2_lh),
    d2_hh = IndirectEffectCalc(mu1$d2_hh),
    d3_lh = IndirectEffectCalc(mu1$d3_lh),
    d3_hh = IndirectEffectCalc(mu1$d3_hh),
    d4_lh = IndirectEffectCalc(mu1$d4_lh),
    d4_hh = IndirectEffectCalc(mu1$d4_hh)
  )
  
  # save results
  results <- list(
    mu0 = mu0, 
    mu1 = mu1, 
    DE = DE, 
    IE0 = IE0, 
    IE1 = IE1
  )
  
  return(results)
}

guncPOIS <- function(sim.k, n_burn = 250000, n_mcmc = 750000, n_post = 250, n_thin = 10){
  # Estimates direct and indirect effects using a Poisson regression outcome
  #   model, with uncertainty propagation from the sulfate module to the 
  #   effect estimates.
  # Input:
  #   sim.k: Multiple imputation sample index.
  #   n_burn: Number of burnin samples (default = 250000).
  #   n_mcmc: Number of MCMC samples to retain after burnin (default = 750000).
  #   n_post: Number of posterior samples to use in effect estimates (default = 250).
  #   n_thin: Thinning rate in MCMC sampler (default = 10).
  # Output:
  #   List containing posterior estimates of mu(z,g), DE(g), and IE(z,g) for
  #     the given theta sample.
  
  ### set seed
  set.seed(95 + sim.k)
  
  ### sample theta
  theta_k <- theta.samples[sim.k,]
  
  ### initalize all other variables 
  outcome = out.df
  fac_df = fac.df
  proj_mat = proj.mat
  em_mat = em.mat
  adv_mats = adv.mats
  scrubber_vec = scrubber.vec
  key_assoc = key.assoc
  
  ### finalize data
  
  # upwind covariates:
  #   degree and log(upwind heat input)
  up_covars <- upwindCovar(
    theta = theta_k,
    mats = adv_mats,
    X_mat = em_mat,
    P_mat = proj_mat,
    s_vec = scrubber_vec,
    key = key_assoc,
    covars = fac_df[, "totHeatInput"],
    version = 2
  )
  
  outcome$degree <- up_covars[[1]]
  outcome$log_up_heat <- up_covars[[2]]
  
  # G_i estimate
  G_k <- gCalc(
    theta = theta_k, 
    mats = adv_mats, 
    X_mat = em_mat, 
    P_mat = proj_mat, 
    s_vec = scrubber_vec, 
    key = key_assoc, 
    version = 3
  )
  
  outcome$G <- G_k
  
  # remove (denom == 0)
  out_final <- outcome[outcome$Ped_16 > 0, ]
  
  ### Format X for BART
  
  # covariates
  X_mat <- as.matrix(out_final %>%
    st_drop_geometry() %>%
    mutate(Z = scrubbed, Z.G = scrubbed * G) %>%
    dplyr::select(c(
      "female", "ped_prop", "median_age", "white_prop",
      "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income",
      "move_rate", "ins_prop", "renter_housing", "urban_prop",
      "log_density", "smokerate", "tmin", "tmax", "prcp", "vp", "rel_humid", "bc", 
      "log_heat", "log_optime", "pct_capacity", "dist_to_key", "log_up_heat", "degree",
      "Z", "G", "Z.G"
  )))
  
  # response
  y <- out_final %>% dplyr::pull(asthma)
  
  # offset
  offset_vec <- out_final %>% dplyr::pull(Ped_16)

  # GLM estimate
  glm.mean <- stats::glm(
    y ~ . + offset(log(denom)) - denom,
    family = stats::poisson, 
    data = data.frame(y, denom = offset_vec, X_mat)
  )
  
  # Poisson regression 
  pois_res <- pois_reg_cpp(
    y_vector = y,
    X_matrix = X_mat,
    offset = log(offset_vec),
    theta_0 = glm.mean$coefficients,
    n_mcmc = n_mcmc,
    thin = n_thin,
    burnin = n_burn,
    n_adapt = 5000,
    Sigma_0 = tuningMatFull()[-23,-23], 
    keep_burnin = FALSE,
    prior_var = 1
  )
  
  poisEstimates(
    fit = pois_res,
    n_post = n_post,
    X = X_mat[, -c((ncol(X_mat) - 2):ncol(X_mat))],
    degree_idx = which(colnames(X_mat) == "degree"),
    dist_idx = which(colnames(X_mat) == "dist_to_key"),
    up_heat_idx = which(colnames(X_mat) == "log_up_heat"),
    heat_idx = which(colnames(X_mat) == "log_heat")
  )
}

bartEstimates <- function(
    fit, n_post, r_mean, X, degree_idx, dist_idx, up_heat_idx, heat_idx,
    degree_cuts = c(0.16, 0.4),
    dist_cuts = c(50, 100, 200),
    heat_cut = 15.5
){
  # Estimates the direct and indirect causal effects using output from a 
  #   log-linear BART regression outcome model. This includes effects stratified
  #   by levels of degree, distance to key-associated facilities, and 
  #   key-associated heat input.
  # Input:
  #   fit: Log-linear BART regression output from ''.
  #   n_post: Number of samples from posterior to include in casual effect estimate.
  #   X: Design matrix.
  #   degree_idx: Weighted degree (covariate matrix column) index.
  #   dist_idx: Distance to key-associated facility index.
  #   up_heat_idx: Upwind heat input index
  #   heat_idx: Key-associated heat input index
  #   degree_cuts: Cut thresholds for degree.
  #   dist_cuts: Cut thresholds for distance.
  #   heat_cut: Cut thresholds for (log-) heat index.
  # Output: 
  #   List containing posterior estimates of mu(z,g), DE(g), and IE(z,g).
  
  # create G sequence
  g_seq <- seq(from = 0.25, to = 0.9, by = 0.01)
  n_g <- length(g_seq)
  
  # stratify by degree
  low_deg <- which(X[, degree_idx] < degree_cuts[1])
  med_deg <- which(X[, degree_idx] >= degree_cuts[1] & X[, degree_idx] < degree_cuts[2])
  high_deg <- which(X[, degree_idx] >= degree_cuts[2])
  
  # stratify by degree and upwind heat input
  ld_lh <- which((X[, degree_idx] < degree_cuts[1]) & (X[, up_heat_idx] < heat_cut))
  ld_hh <- which((X[, degree_idx] < degree_cuts[1]) & (X[, up_heat_idx] >= heat_cut))
  md_lh <- which((X[, degree_idx] >= degree_cuts[1] & X[, degree_idx] < degree_cuts[2]) & (X[, up_heat_idx] < heat_cut))
  md_hh <- which((X[, degree_idx] >= degree_cuts[1] & X[, degree_idx] < degree_cuts[2]) & (X[, up_heat_idx] >= heat_cut))
  hd_lh <- which((X[, degree_idx] >= degree_cuts[2]) & (X[, up_heat_idx] < heat_cut))
  hd_hh <- which((X[, degree_idx] >= degree_cuts[2]) & (X[, up_heat_idx] >= heat_cut))
  
  # stratify by distance to key associated facility
  dist1 <- which(X[, dist_idx] < dist_cuts[1])
  dist2 <- which(X[, dist_idx] >= dist_cuts[1] & X[, dist_idx] < dist_cuts[2])
  dist3 <- which(X[, dist_idx] >= dist_cuts[2] & X[, dist_idx] < dist_cuts[3])
  dist4 <- which(X[, dist_idx] >= dist_cuts[3])
  
  # stratify by distance and heat input
  d1_lh <- which((X[, dist_idx] < dist_cuts[1]) & (X[, heat_idx] < heat_cut))
  d1_hh <- which((X[, dist_idx] < dist_cuts[1]) & (X[, heat_idx] >= heat_cut))
  d2_lh <- which((X[, dist_idx] >= dist_cuts[1] & X[, dist_idx] < dist_cuts[2]) & (X[, heat_idx] < heat_cut))
  d2_hh <- which((X[, dist_idx] >= dist_cuts[1] & X[, dist_idx] < dist_cuts[2]) & (X[, heat_idx] >= heat_cut))
  d3_lh <- which((X[, dist_idx] >= dist_cuts[2] & X[, dist_idx] < dist_cuts[3]) & (X[, heat_idx] < heat_cut))
  d3_hh <- which((X[, dist_idx] >= dist_cuts[2] & X[, dist_idx] < dist_cuts[3]) & (X[, heat_idx] >= heat_cut))
  d4_lh <- which((X[, dist_idx] >= dist_cuts[3]) & (X[, heat_idx] < heat_cut))
  d4_hh <- which((X[, dist_idx] >= dist_cuts[3]) & (X[, heat_idx] >= heat_cut))
  
  # mu0 list
  mu0 <- list(
    all = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ldeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    mdeg = matrix(NA_real_, nrow = n_post, ncol = n_g), 
    hdeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_hh = matrix(NA_real_, nrow = n_post, ncol = n_g)
  )
  
  # m1 list
  mu1 <- list(
    all = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ldeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    mdeg = matrix(NA_real_, nrow = n_post, ncol = n_g), 
    hdeg = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    ld_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    md_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    hd_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4 = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d1_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d2_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d3_hh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_lh = matrix(NA_real_, nrow = n_post, ncol = n_g),
    d4_hh = matrix(NA_real_, nrow = n_post, ncol = n_g)
  )
  
  ### estimate mu(z,g) for different strata
  
  # z = 0
  z <- 0
  
  for (k in 1:n_g){
    
    g_val <- g_seq[k]
    X_m <- cbind(X, z, g_val)
    
    colnames(X_m)[(ncol(X_m) - 1):ncol(X_m)] <- c("Z", "G")
    out <- fit$tree_fit$tree_samples$coefs(as.matrix(t(X_m)))
    f_hat <- exp(out[1, , ] + log(r_mean))
    
    # mu(z,g)
    mu0$all[, k] <- colMeans(f_hat)
    if (length(low_deg) > 1) mu0$ldeg[, k] <- colMeans(f_hat[low_deg, ])
    if (length(med_deg) > 1) mu0$mdeg[, k] <- colMeans(f_hat[med_deg, ])
    if (length(high_deg) > 1) mu0$hdeg[, k] <- colMeans(f_hat[high_deg, ])
    if (length(ld_lh) > 1) mu0$ld_lh[, k] <- colMeans(f_hat[ld_lh, ])
    if (length(ld_hh) > 1) mu0$ld_hh[, k] <- colMeans(f_hat[ld_hh, ])
    if (length(md_lh) > 1) mu0$md_lh[, k] <- colMeans(f_hat[md_lh, ])
    if (length(md_hh) > 1) mu0$md_hh[, k] <- colMeans(f_hat[md_hh, ])
    if (length(hd_lh) > 1) mu0$hd_lh[, k] <- colMeans(f_hat[hd_lh, ])
    if (length(hd_hh) > 1) mu0$hd_hh[, k] <- colMeans(f_hat[hd_hh, ])
    if (length(dist1) > 1) mu0$d1[, k] <- colMeans(f_hat[dist1, ])
    if (length(dist2) > 1) mu0$d2[, k] <- colMeans(f_hat[dist2, ])
    if (length(dist3) > 1) mu0$d3[, k] <- colMeans(f_hat[dist3, ])
    if (length(dist4) > 1) mu0$d4[, k] <- colMeans(f_hat[dist4, ])
    if (length(d1_lh) > 1) mu0$d1_lh[, k] <- colMeans(f_hat[d1_lh, ])
    if (length(d1_hh) > 1) mu0$d1_hh[, k] <- colMeans(f_hat[d1_hh, ])
    if (length(d2_lh) > 1) mu0$d2_lh[, k] <- colMeans(f_hat[d2_lh, ])
    if (length(d2_hh) > 1) mu0$d2_hh[, k] <- colMeans(f_hat[d2_hh, ])
    if (length(d3_lh) > 1) mu0$d3_lh[, k] <- colMeans(f_hat[d3_lh, ])
    if (length(d3_hh) > 1) mu0$d3_hh[, k] <- colMeans(f_hat[d3_hh, ])
    if (length(d4_lh) > 1) mu0$d4_lh[, k] <- colMeans(f_hat[d4_lh, ])
    if (length(d4_hh) > 1) mu0$d4_hh[, k] <- colMeans(f_hat[d4_hh, ])
  }

  # z = 1 
  z <- 1
  
  for (k in 1:n_g){
    
    g_val <- g_seq[k]
    X_m <- cbind(X, z, g_val)
    
    colnames(X_m)[(ncol(X_m) - 1):ncol(X_m)] <- c("Z", "G")
    out <- fit$tree_fit$tree_samples$coefs(as.matrix(t(X_m)))
    f_hat <- exp(out[1, , ] + log(r_mean))
    
    # mu(z,g)
    mu1$all[, k] <- colMeans(f_hat)
    if (length(low_deg) > 1) mu1$ldeg[, k] <- colMeans(f_hat[low_deg, ])
    if (length(med_deg) > 1) mu1$mdeg[, k] <- colMeans(f_hat[med_deg, ])
    if (length(high_deg) > 1) mu1$hdeg[, k] <- colMeans(f_hat[high_deg, ])
    if (length(ld_lh) > 1) mu1$ld_lh[, k] <- colMeans(f_hat[ld_lh, ])
    if (length(ld_hh) > 1) mu1$ld_hh[, k] <- colMeans(f_hat[ld_hh, ])
    if (length(md_lh) > 1) mu1$md_lh[, k] <- colMeans(f_hat[md_lh, ])
    if (length(md_hh) > 1) mu1$md_hh[, k] <- colMeans(f_hat[md_hh, ])
    if (length(hd_lh) > 1) mu1$hd_lh[, k] <- colMeans(f_hat[hd_lh, ])
    if (length(hd_hh) > 1) mu1$hd_hh[, k] <- colMeans(f_hat[hd_hh, ])
    if (length(dist1) > 1) mu1$d1[, k] <- colMeans(f_hat[dist1, ])
    if (length(dist2) > 1) mu1$d2[, k] <- colMeans(f_hat[dist2, ])
    if (length(dist3) > 1) mu1$d3[, k] <- colMeans(f_hat[dist3, ])
    if (length(dist4) > 1) mu1$d4[, k] <- colMeans(f_hat[dist4, ])
    if (length(d1_lh) > 1) mu1$d1_lh[, k] <- colMeans(f_hat[d1_lh, ])
    if (length(d1_hh) > 1) mu1$d1_hh[, k] <- colMeans(f_hat[d1_hh, ])
    if (length(d2_lh) > 1) mu1$d2_lh[, k] <- colMeans(f_hat[d2_lh, ])
    if (length(d2_hh) > 1) mu1$d2_hh[, k] <- colMeans(f_hat[d2_hh, ])
    if (length(d3_lh) > 1) mu1$d3_lh[, k] <- colMeans(f_hat[d3_lh, ])
    if (length(d3_hh) > 1) mu1$d3_hh[, k] <- colMeans(f_hat[d3_hh, ])
    if (length(d4_lh) > 1) mu1$d4_lh[, k] <- colMeans(f_hat[d4_lh, ])
    if (length(d4_hh) > 1) mu1$d4_hh[, k] <- colMeans(f_hat[d4_hh, ])
  }

  # estimate DE(g) = mu(1, g) - mu(0, g)
  DE <- list(
    all = mu1$all - mu0$all,
    ldeg = mu1$ldeg - mu0$ldeg,
    mdeg = mu1$mdeg - mu0$mdeg,
    hdeg = mu1$hdeg - mu0$hdeg,
    ld_lh = mu1$ld_lh - mu0$ld_lh,
    ld_hh = mu1$ld_hh - mu0$ld_hh,
    md_lh = mu1$md_lh - mu0$md_lh,
    md_hh = mu1$md_hh - mu0$md_hh,
    hd_lh = mu1$hd_lh - mu0$hd_lh,
    hd_hh = mu1$hd_hh - mu0$hd_hh,
    d1 = mu1$d1 - mu0$d1,
    d2 = mu1$d2 - mu0$d2,
    d3 = mu1$d3 - mu0$d3,
    d4 = mu1$d4 - mu0$d4,
    d1_lh = mu1$d1_lh - mu0$d1_lh,
    d1_hh = mu1$d1_hh - mu0$d1_hh,
    d2_lh = mu1$d2_lh - mu0$d2_lh,
    d2_hh = mu1$d2_hh - mu0$d2_hh,
    d3_lh = mu1$d3_lh - mu0$d3_lh,
    d3_hh = mu1$d3_hh - mu0$d3_hh,
    d4_lh = mu1$d4_lh - mu0$d4_lh,
    d4_hh = mu1$d4_hh - mu0$d4_hh
  )
  
  # IE(0, g) = mu(0, g) - mu(0, g_min)
  IE0 <- list(
    all = IndirectEffectCalc(mu0$all),
    ldeg = IndirectEffectCalc(mu0$ldeg), 
    mdeg = IndirectEffectCalc(mu0$mdeg),
    hdeg = IndirectEffectCalc(mu0$hdeg),
    ld_lh = IndirectEffectCalc(mu0$ld_lh),
    ld_hh = IndirectEffectCalc(mu0$ld_hh),
    md_lh = IndirectEffectCalc(mu0$md_lh),
    md_hh = IndirectEffectCalc(mu0$md_hh),
    hd_lh = IndirectEffectCalc(mu0$hd_lh),
    hd_hh = IndirectEffectCalc(mu0$hd_hh),
    d1 = IndirectEffectCalc(mu0$d1),
    d2 = IndirectEffectCalc(mu0$d2),
    d3 = IndirectEffectCalc(mu0$d3),
    d4 = IndirectEffectCalc(mu0$d4),
    d1_lh = IndirectEffectCalc(mu0$d1_lh),
    d1_hh = IndirectEffectCalc(mu0$d1_hh),
    d2_lh = IndirectEffectCalc(mu0$d2_lh),
    d2_hh = IndirectEffectCalc(mu0$d2_hh),
    d3_lh = IndirectEffectCalc(mu0$d3_lh),
    d3_hh = IndirectEffectCalc(mu0$d3_hh),
    d4_lh = IndirectEffectCalc(mu0$d4_lh),
    d4_hh = IndirectEffectCalc(mu0$d4_hh)
  )
  
  # IE(1, g) = mu(1, g) - mu(1, g_min)
  IE1 <- list(
    all = IndirectEffectCalc(mu1$all),
    ldeg = IndirectEffectCalc(mu1$ldeg), 
    mdeg = IndirectEffectCalc(mu1$mdeg),
    hdeg = IndirectEffectCalc(mu1$hdeg),
    ld_lh = IndirectEffectCalc(mu1$ld_lh),
    ld_hh = IndirectEffectCalc(mu1$ld_hh),
    md_lh = IndirectEffectCalc(mu1$md_lh),
    md_hh = IndirectEffectCalc(mu1$md_hh),
    hd_lh = IndirectEffectCalc(mu1$hd_lh),
    hd_hh = IndirectEffectCalc(mu1$hd_hh),
    d1 = IndirectEffectCalc(mu1$d1),
    d2 = IndirectEffectCalc(mu1$d2),
    d3 = IndirectEffectCalc(mu1$d3),
    d4 = IndirectEffectCalc(mu1$d4),
    d1_lh = IndirectEffectCalc(mu1$d1_lh),
    d1_hh = IndirectEffectCalc(mu1$d1_hh),
    d2_lh = IndirectEffectCalc(mu1$d2_lh),
    d2_hh = IndirectEffectCalc(mu1$d2_hh),
    d3_lh = IndirectEffectCalc(mu1$d3_lh),
    d3_hh = IndirectEffectCalc(mu1$d3_hh),
    d4_lh = IndirectEffectCalc(mu1$d4_lh),
    d4_hh = IndirectEffectCalc(mu1$d4_hh)
  )
  
  # save results
  results <- list(
    mu0 = mu0, 
    mu1 = mu1, 
    DE = DE, 
    IE0 = IE0, 
    IE1 = IE1
  )
  
  return(results)
}

guncBART <- function(sim.k, n_burn = 5000, n_post = 250, n_thin = 5){
  # Estimates direct and indirect effects using a log-linear BART regression 
  #   outcome model, with uncertainty propagation from the sulfate module to
  #    the effect estimates.
  # Input:
  #   sim.k: Multiple imputation sample index.
  #   n_burn: Number of burnin samples (default = 5000).
  #   n_post: Number of posterior samples to use in effect estimates (default = 250).
  #   n_thin: Thinning rate in MCMC sampler (default = 5).
  # Output:
  #   List containing posterior estimates of mu(z,g), DE(g), and IE(z,g) for
  #     the given theta sample.
  
  # set seed
  set.seed(95 + sim.k)
  
  # sample theta
  theta_k <- theta.samples[sim.k,]
  
  # initalize all other variables 
  outcome = out.df
  fac_df = fac.df
  proj_mat = proj.mat
  em_mat = em.mat
  adv_mats = adv.mats
  scrubber_vec = scrubber.vec
  key_assoc = key.assoc
  
  # upwind covariates:
  #   degree and log(upwind heat input)
  up_covars <- upwindCovar(
    theta = theta_k,
    mats = adv_mats,
    X_mat = em_mat,
    P_mat = proj_mat,
    s_vec = scrubber_vec,
    key = key_assoc,
    covars = fac_df[, "totHeatInput"],
    version = 2
  )
  
  outcome$degree <- up_covars[[1]]
  outcome$log_up_heat <- up_covars[[2]]
  
  # G_i estimate
  G_k <- gCalc(
    theta = theta_k, 
    mats = adv_mats, 
    X_mat = em_mat, 
    P_mat = proj_mat, 
    s_vec = scrubber_vec, 
    key = key_assoc, 
    version = 3
  )
  
  outcome$G <- G_k
  
  # remove (denom == 0)
  out_final <- outcome[outcome$Ped_16 > 0, ]

  # covariates
  X_mat <- as.matrix(out_final %>%
    st_drop_geometry() %>%
    mutate(Z = scrubbed) %>%
    dplyr::select(c(
      "Ped_16", "female", "ped_prop", "median_age", "white_prop",
      "black_prop", "hisp_prop", "hs_prop", "pov_prop", "log_income",
      "move_rate", "ins_prop", "renter_housing", "urban_prop",
      "log_density", "smokerate", "tmin", "tmax", "prcp", "vp", "rel_humid", "bc", 
      "log_heat", "log_optime", "pct_capacity", "dist_to_key", "log_up_heat", "degree",
      "Z", "G"
    )))
  
  # response
  y <- out_final %>% dplyr::pull(asthma)
  
  # offset
  offset_vec <-  log(out_final %>% dplyr::pull(Ped_16))
  
  # determine leaf prior hyperparameter
  rates <- y / exp(offset_vec)
  r_star <- quantile(rates, probs = 0.975)
  r_mean <- mean(rates)
  a_0 <- 0.5 * (log(r_star) - log(r_mean))
  
  ### Fit BART
  
  # Poisson regression with BART
  bart_res <- count_bart(
    y = y,
    x = X_mat,
    offset = exp(offset_vec + log(r_mean)),
    nburn = n_burn,
    nsim = n_post,
    nthin = n_thin,
    update_interval = 1000,
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
  
  # Estimate causal effects
  bartEstimates(
    fit = bart_res, n_post = n_post, r_mean = r_mean,
    X = X_mat[, -c(ncol(X_mat) - 1, ncol(X_mat))],
    degree_idx = which(colnames(X_mat) == "degree"),
    dist_idx = which(colnames(X_mat) == "dist_to_key"),
    up_heat_idx = which(colnames(X_mat) == "log_up_heat"),
    heat_idx = which(colnames(X_mat) == "log_heat")
  )
}
