# Functions for numeric calculations
# sezenismail@gmail.com 2015-10-28
#

# This function is similar to JAC2 function but
# up and down boundaries numerical calculation are different.
JAC3 <- function(b, c, dx, dy=dx) {
  
  if(!is.matrix(b) || !is.matrix(c)) stop("b and c both must be matrix")
  if(all(dim(b) != dim(c))) stop("Dimensions b and c must be equal.")
  
  l   <- nrow(b) # lines
  m   <- ncol(b) # columns
  JJ <- matrix(0,nrow=l, ncol=m)
  
  if(l<3) stop("Number of rows of b and c must be equal or greater than 3.")
  if(m<3) stop("Number of columns of b and c must be equal or greater than 3.")
  
  J_first_col <- function() {
    i <- 2:(l+1)
    rows <- c(l,1:l,1);cols <- 1:2
    b <- b[rows, cols]; c <- c[rows, cols]
    return((b[i,1]   + b[i+1,1]  - b[i,2]    - b[i+1,2]) * (c[i,1]   + c[i+1,1]) -
             (b[i-1,1] + b[i,1]    - b[i-1,2]  - b[i,2])   * (c[i-1,1] + c[i,1]) +
             (b[i+1,1] + b[i+1,2]  - b[i-1,1]  - b[i-1,2]) * (c[i,1]   + c[i,2]) +
             (b[i+1,1] - b[i,2])   * (c[i,1]   + c[i+1,2]) +
             (b[i,2]   - b[i-1,1]) * (c[i-1,2] + c[i,1]))
  }
  
  J_last_col <- function() {
    i <- 2:(l+1)
    rows <- c(l,1:l,1);cols <- (m-1):m
    b <- b[rows, cols]; c <- c[rows, cols]
    return((b[i,1]   + b[i+1,1] - b[i,2]   - b[i+1,2]) * (c[i,2]   + c[i+1,2]) -
             (b[i-1,1] + b[i,1]   - b[i-1,2] - b[i,2])   * (c[i-1,2] + c[i,2]) -
             (b[i+1,1] + b[i+1,2] - b[i-1,1] - b[i-1,2]) * (c[i,1]   + c[i,2]) -
             (b[i,1]   - b[i-1,2]) * (c[i-1,1] + c[i,2]) +
             (b[i,1]   - b[i+1,2]) * (c[i,2]   + c[i+1,1]))
  }
  
  J1 <- function(b,c,i,j) {
    return((b[i,j-1]   + b[i+1,j-1] - b[i,j+1]   - b[i+1,j+1]) * (c[i+1,j] - c[i,j]) +
             (b[i-1,j-1] + b[i,j-1]   - b[i-1,j+1] - b[i,j+1])   * (c[i,j]   - c[i-1,j]) +
             (b[i+1,j]   + b[i+1,j+1] - b[i-1,j]   - b[i-1,j+1]) * (c[i,j+1] - c[i,j]) +
             (b[i+1,j-1] + b[i+1,j]   - b[i-1,j-1] - b[i-1,j])   * (c[i,j]   - c[i,j-1]) +
             (b[i+1,j]   - b[i,j+1]) * (c[i+1,j+1] - c[i,j]) +
             (b[i,j-1]   - b[i-1,j]) * (c[i,j]     - c[i-1,j-1]) +
             (b[i,j+1]   - b[i-1,j]) * (c[i-1,j+1] - c[i,j]) +
             (b[i+1,j]   - b[i,j-1]) * (c[i,j]     - c[i+1,j-1]))
  }
  
  j <- 2:(m-1)
  JJ[1,j] <- J1(b[c(l,1,2),], c[c(l,1,2),], 2, j) # First Row
  JJ[l,j] <- J1(b[c((l-1),l,1),], c[c((l-1),l,1),], 2, j) # Last Row
  i <- 2:(l-1)
  JJ[i,j] <- J1(b,c,i,j) # inside grids
  JJ[1:l,1] <- J_first_col() # First Col
  JJ[1:l,m] <- J_last_col() # Last Col
  return(JJ/(12*dy*dx))
}

# This function is R port of FORTRAN algorithm
# implemented in JACOBIAN.FOR and NICKER.FOR.
JAC2 <- function(b, c, dx, dy=dx) {
  
  if(!is.matrix(b) || !is.matrix(c)) stop("b and c both must be matrix")
  if(all(dim(b) != dim(c))) stop("Dimensions b and c must be equal.")
  
  l   <- nrow(b) # lines
  m   <- ncol(b) # columns
  JJ <- matrix(0,nrow=l, ncol=m)
  
  if(l<3) stop("Number of rows of b and c must be equal or greater than 3.")
  if(m<3) stop("Number of columns of b and c must be equal or greater than 3.")
  
  J_first_col <- function() {
    i <- 2:(l+1)
    rows <- c((l-1),1:l,2); cols <- 1:2
    b <- b[rows, cols]; c <- c[rows, cols]
    return((b[i,1]   + b[i+1,1]  - b[i,2]    - b[i+1,2]) * (c[i,1]   + c[i+1,1]) -
           (b[i-1,1] + b[i,1]    - b[i-1,2]  - b[i,2])   * (c[i-1,1] + c[i,1]) +
           (b[i+1,1] + b[i+1,2]  - b[i-1,1]  - b[i-1,2]) * (c[i,1]   + c[i,2]) +
           (b[i+1,1] - b[i,2])   * (c[i,1]   + c[i+1,2]) +
           (b[i,2]   - b[i-1,1]) * (c[i-1,2] + c[i,1]))
  }
  
  J_last_col <- function() {
    i <- 2:(l+1)
    rows <- c((l-1),1:l,2); cols <- (m-1):m
    b <- b[rows, cols]; c <- c[rows, cols]
    return((b[i,1]   + b[i+1,1] - b[i,2]   - b[i+1,2]) * (c[i,2]   + c[i+1,2]) -
           (b[i-1,1] + b[i,1]   - b[i-1,2] - b[i,2])   * (c[i-1,2] + c[i,2]) -
           (b[i+1,1] + b[i+1,2] - b[i-1,1] - b[i-1,2]) * (c[i,1]   + c[i,2]) -
           (b[i,1]   - b[i-1,2]) * (c[i-1,1] + c[i,2]) +
           (b[i,1]   - b[i+1,2]) * (c[i,2]   + c[i+1,1]))
  }
  
  J1 <- function(b,c,i,j) {
    return((b[i,j-1]   + b[i+1,j-1] - b[i,j+1]   - b[i+1,j+1]) * (c[i+1,j] - c[i,j]) +
           (b[i-1,j-1] + b[i,j-1]   - b[i-1,j+1] - b[i,j+1])   * (c[i,j]   - c[i-1,j]) +
           (b[i+1,j]   + b[i+1,j+1] - b[i-1,j]   - b[i-1,j+1]) * (c[i,j+1] - c[i,j]) +
           (b[i+1,j-1] + b[i+1,j]   - b[i-1,j-1] - b[i-1,j])   * (c[i,j]   - c[i,j-1]) +
           (b[i+1,j]   - b[i,j+1]) * (c[i+1,j+1] - c[i,j]) +
           (b[i,j-1]   - b[i-1,j]) * (c[i,j]     - c[i-1,j-1]) +
           (b[i,j+1]   - b[i-1,j]) * (c[i-1,j+1] - c[i,j]) +
           (b[i+1,j]   - b[i,j-1]) * (c[i,j]     - c[i+1,j-1]))
  }

  j <- 2:(m-1)
  JJ[1,j] <- J1(b[c((l-1),1,2),], c[c((l-1),1,2),], 2, j) # First Row
  JJ[l,j] <- J1(b[c((l-1),l,2),], c[c((l-1),l,2),], 2, j) # Last Row
  
  
  i <- 2:(l-1)
  JJ[i,j] <- J1(b,c,i,j) # inside grids
  JJ[1:l,j] <- JJ[1:l,j]/(12*dy*dx)
  
  JJ[1:l,1] <- J_first_col()/(12*dy*dx) # First Col
  JJ[1:l,m] <- J_last_col()/(12*dy*dx) # Last Col
  return(JJ)
}

# This function has the same jacobian algorithm 
# with JACOBIAN.FOR and NICKER.FOR.
JAC1 <- function(a, b, c, dx, dy, l, m, l1, m1, l2, m2) {
  # k <- 1
  for(j in 2:m1) {
    dm <- 12*dy*dx[j]
    for(i in 1:l) {
      if((i-1)<=0) {
        im1 <- l1
        ip1 <- 2
      }else {
        if((i-l)<0) {
          im1 <- i-1
          ip1 <- i+1
        } else {
          im1 <- l1
          ip1 <- 2
        }
      }
#       if((im1)>i || (ip1)<i) {
#         rw <- if(im1>i) c(im1, im1+1, i, ip1) else c(im1, i, ip1-1 ,ip1)
#         cl <- c(j-1, j, j+1)
#         
#         mat <- b[rw, cl]
#         colnames(mat) <- cl
#         if(im1>i) rw[3] <- paste0(rw[3], "*") else rw[2] <- paste0(rw[2], "*")
#         rownames(mat) <- rw
#         print(mat)
#       }
      
      a[i,j] <- (b[i,j-1]   + b[ip1,j-1] - b[i,j+1]   - b[ip1,j+1]) * (c[ip1,j] - c[i,j]) +
                (b[im1,j-1] + b[i,j-1]   - b[im1,j+1] - b[i,j+1])   * (c[i,j]   - c[im1,j]) +
                (b[ip1,j]   + b[ip1,j+1] - b[im1,j]   - b[im1,j+1]) * (c[i,j+1] - c[i,j]) +
                (b[ip1,j-1] + b[ip1,j]   - b[im1,j-1] - b[im1,j])   * (c[i,j]   - c[i,j-1]) +
                (b[ip1,j]   - b[i,j+1]) * (c[ip1,j+1] - c[i,j]) +
                (b[i,j-1]   - b[im1,j]) * (c[i,j]     - c[im1,j-1]) +
                (b[i,j+1]   - b[im1,j]) * (c[im1,j+1] - c[i,j]) +
                (b[ip1,j]   - b[i,j-1]) * (c[i,j]     - c[ip1,j-1])
      #
      
      a[i,j] <- a[i,j]/dm
    }
  }
  # a <- matrix(0, nrow=l, ncol=m)
  #
  for(i in 1:l) {
    if((i-1)<1) {
      im1 <- l1
      ip1 <- 2
    } else {
      if((i-l)<0) {
        im1 <- i-1
        ip1 <- i+1
      } else {
        im1 <- l1
        ip1 <- 2
      }
    }
    # cat(im1, i, ip1, b[im1,1], b[i,1],b[ip1,1],"\n")
    dm <- 12*dy*dx[1]
    #
    a[i,1]  <- (b[i,1]   + b[ip1,1]  - b[i,2]    - b[ip1,2]) * (c[i,1]   + c[ip1,1]) -
               (b[im1,1] + b[i,1]    - b[im1,2]  - b[i,2])   * (c[im1,1] + c[i,1]) +
               (b[ip1,1] + b[ip1,2]  - b[im1,1]  - b[im1,2]) * (c[i,1]   + c[i,2]) +
               (b[ip1,1] - b[i,2])   * (c[i,1]   + c[ip1,2]) +
               (b[i,2]   - b[im1,1]) * (c[im1,2] + c[i,1])
    
    a[i,1] <- a[i,1]/dm
    dm     <- 12*dy*dx[m]
    
    a[i,m]  <- (b[i,m-1]   + b[ip1,m-1] - b[i,m]      - b[ip1,m]) * (c[i,m]   + c[ip1,m]) -
               (b[im1,m-1] + b[i,m-1]   - b[im1,m]    - b[i,m])   * (c[im1,m] + c[i,m]) -
               (b[ip1,m-1] + b[ip1,m]   - b[im1,m-1]  - b[im1,m]) * (c[i,m-1] + c[i,m]) -
               (b[i,m-1]   - b[im1,m])  * (c[im1,m-1] + c[i,m]) -
               (b[ip1,m]   - b[i,m-1])  * (c[i,m]     + c[ip1,m-1])
    
    a[i,m]  <- a[i,m]/dm 
  }
  # print(a)
  # stop()
  return(a)
}

LAPLAC <- function(t, dx) {
  da <- dx*dx
  n <- nrow(t); m <- ncol(t)
  difu <- matrix(0, nrow=n, ncol=m)
  n1 <- n-1; m1 <- m-1
  # Interior Points
  for(i in 2:n1) {
    for(j in 2:m1) {
      difu[i,j] <- (t[i,j+1] + t[i,j-1] + t[i+1,j] + t[i-1,j] - 4*t[i,j])/da
    }
  }
  # Side boundaries
  for(j in 2:m1) {
    difu[1,j] <- (2*t[2,j]  + t[1,j+1] + t[1,j-1] - 4*t[1,j])/da
    difu[n,j] <- (2*t[n1,j] + t[n,j+1] + t[n,j-1] - 4*t[n,j])/da
  }
  # Top and Bottom Boundaries
  for(i in 2:n1) {
    difu[i,1] <- (2*t[i,2]  + t[i-1,1] + t[i+1,1] - 4*t[i,1])/da                    
    difu[i,m] <- (2*t[i,m1] + t[i-1,m] + t[i+1,m] - 4*t[i,m])/da
  }
  db <- 2/da
  #Four Corners
  difu[1,1] <- db*(t[1,2]  + t[2,1]  - 2*t[1,1])                                  
  difu[1,m] <- db*(t[1,m1] + t[2,m]  - 2*t[1,m])                                 
  difu[n,m] <- db*(t[n1,m] + t[n,m1] - 2*t[n,m])                                
  difu[n,1] <- db*(t[n1,1] + t[n,2]  - 2*t[n,1])
  return(difu)
}

RELAX1 <- function(psi, eta, n, m, nscan, alpha){
  eps   <- 1e-03
  dx    <- 10
  nscan <- 0
  n1    <- n-1
  m1    <- m-1
  repeat {
    nscan <- nscan + 1
    rmax <- 0
    for(i in 2:n1) {
      for(j in 2:m1) {
        r1 <-  0.25*(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1])
        r2 <- -psi[i,j] - 0.25*(dx^2)*eta[i,j]
        r  <- r1+r2
        if(rmax < abs(r)){
          isave <- i
          jsave <- j
          rmax <- abs(r)
        }
        psi[i,j] <- psi[i,j] + alpha*r
      }
    }
    if(nscan == 300) {
      stop("No convergence achieved")
    }
    ps <- abs(psi[isave, jsave])
    if((rmax/ps-eps)<0) break;
  }
  return(psi)
}


