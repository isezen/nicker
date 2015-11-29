#
# This program simulates the evolution of a simple
# cloud model as described by nickerson (1965).
# The heat source is inserted at the bottom of the
# atmosphere and the fluid is initially stratified
# and at rest.
# The thermal begins to ascent,and a vortex circulation
# develops. The temperature within the thermal
# increases due to heating and  the thermal rises above the
# heat source. The maximum temperature reaches a peak
# then decreases when dissipation exceed buoyancy.
# The potential temperature field resembles a thin
# stemmed mushroom.
#

source("numeric.R")
source("plot.helper.R")
JAC <- JAC3

integral <- function(mat, dx, dy=dx) {
  n <- nrow(mat)-1; m <- ncol(mat)-1 # n, m must be even numbers.
  x0 <- y0 <- 1
  I <- vector(length=16)

                      I[1] <-   mat[x0,y0]
  i <- (1:(n/2-1))+1; I[2] <- 2*sum(mat[2*i,y0])
  i <- (1:(n/2))+1;   I[3] <- 4*sum(mat[2*i-1,y0])
                      I[4] <-   mat[n,y0]
                        I1 <-   sum(I[1:4])
  
  j <- (1:(m/2-1))+1;                     I[5] <-   sum(mat[   x0, 2*j])
  j <- (1:(m/2-1))+1; i <- (1:(n/2-1))+1; I[6] <- 2*sum(mat[  2*i, 2*j])
  j <- (1:(m/2-1))+1; i <- (1:(n/2))+1;   I[7] <- 4*sum(mat[2*i-1, 2*j])
  j <- (1:(m/2-1))+1;                     I[8] <-   sum(mat[    n, 2*j])
                                            I2 <- 2*sum(I[5:8])
  
  j <- (1:(m/2))+1;                      I[9] <-   sum(mat[   x0, 2*j-1])
  j <- (1:(m/2))+1; i <- (1:(n/2-1))+1; I[10] <- 2*sum(mat[  2*i, 2*j-1])
  j <- (1:(m/2))+1; i <- (1:(n/2))+1;   I[11] <- 4*sum(mat[2*i-1, 2*j-1])
  j <- (1:(m/2))+1;                     I[12] <-   sum(mat[    n, 2*j-1])
                                           I3 <- 4*sum(I[9:12])
  
                      I[13] <-   mat[x0,m]
  i <- (1:(n/2-1))+1; I[14] <- 2*sum(mat[2*i,m])
  i <- (1:(n/2))+1;   I[15] <- 4*sum(mat[2*i-1,m])
                      I[16] <-   mat[n,m]
                         I4 <-   sum(I[13:16])

  return((dx*dy/9)*(I1+I2+I3+I4))
}

# dx  : Grid resolution
# dt  : time-step
# change_dt_on : change dt value at this time stamp (seconds)
# Q0  : Time invariant non-adiabatic heating (deg /sec) coeffcient
# gnu : eddy kinematic coefficient of heat and momentum
# dx  : grid spacing. 10 meters in horizontal and vertical
# HEATING_IS_TIME_DEPENDENT : Make non-adiabatic heating is time dependent (Cool instead of heat in the middle of run time)
nicker <- function(dx=10, dt=1, change_dt_on=0,  Q0=4e-03, gnu=0.5, 
                   reduce.gradient=F, reduction.fact=1.5, HEATING_IS_TIME_DEPENDENT=F) {
  d_dx       <- function(m, dx) { i<-2:(n-1);(m[i+1,]-m[i-1,])/(2*dx) }
  indices_at <- function(x, min, max) {which((x>=min & x<=max), arr.ind=T)}
  # ----------------------- 
  dts      <- dt  # you can define different time intervals (seconds)
  END_TIME <- 600 # finish the run at this time
  
  x     <- seq(0, 380, dx) # x axis values
  z     <- seq(0, 600, dx) # z axis values
  n     <- length(x)       # vertical dimension
  m     <- length(z)       # horizontal dimensions
  tneut <- 300             # potential temperature of neutral environment
  g     <- 9.81            # accerelation of gravity
  alpha <- 1.89            # relaxation factor
  
  # INITIALIZATIONS
  # psi  : stream function 
  # q    : non-adiabatic heating (deg /sec)
  # u,w  : u (horizontal) and w (vertical) wind speeds
  # eta  : vorticity (per sec)
  # t    : temperature excess (deg celcius)
  #
  # #  initially u = w = psi = eta  = 0. 
  psi <- q <- u <- w <- array(0,c(n,m))
  eta <- t <- array(0,c(n, m, 2))
  
  # Holder for calculation results
  tp <- c(change_dt_on, END_TIME)
  LOOP_MAX <- sum((tp[2:length(tp)]-tp[1:(length(tp)-1)])/dts)
  l <- LOOP_MAX + 1
  res <- list(x=x, z=z, time=vector(length=l), 
              eta=array(0,c(n,m,l)), psi=array(0,c(n,m,l)), t=array(0,c(n,m,l)), 
              phi=array(0,c(n,m,l)), u=array(0,c(n,m,l)), w=array(0,c(n,m,l)))
  
  #  The boundary conditions used are those for freeslip, insulated
  #  surfaces.                                                 
  #  hort vel = 0. at x = 0 and x = 380                                            
  #  vert vel = 0. at z = 0 and z = 600                                            
  #                                            
  #  The temperature excess defined in the x-dir from 0 to 160m
  #  and in the z-direction from 100m to 300 m. Everywhere else
  #  is set to zero .The heating field is defined in x-dir from
  #  0 to 160 m and in z-dir from 80 to 120 m .It is set to 0
  #  elsewhere.
            i <- indices_at(x, 0, 160)   # The thermal having a maximum potential temperature excess of   
            j <- indices_at(z, 100, 300) # 0.5 deg celcius is inserted into a neutral environment initially at the rest.
   t[i, j, 1] <- outer(x[i], z[j], function(x,z) 0.5*cos(pi*x/320)*(cos(pi*(z-100)/400))^2)
   if(reduce.gradient) t[i, 10, 1] <- reduction.fact*t[i, 11, 1]  # decrease sharp gradiant at the bottom of buble.
       t[,,1] <- t[,,1]/tneut # Convert t to phi by dividing t by tneut (t is phi anymore.)
  
  # Save initial situation
  res$time[1]  <- 0
  res$phi[,,1] <- t[,,1]
  res$eta[,,1] <- eta[,,1]
  res$psi[,,1] <- psi
  res$t[,,1]   <- t[,,1]*tneut
  res$u[,,1]   <- u
  res$w[,,1]   <- w
  
  j <- indices_at(z, 80, 120) # Time invariant non-adiabatic heating (deg /sec) is defined as:
  q[i,j] <- Q0*outer(x[i], z[j], function(x,z) cos(pi*x/320)*(cos(pi*(z-100)/40))^2)
  filled.contour(x,z,q, color=colorRampPalette(c("white", "forestgreen","gold3","orange","firebrick4")), 
  xlab="x (meters)", ylab="z(meters)", key.title=title(line=1, main=expression(Q)), 
  main="Q Non-adiabatic Heating (degree/sec)")
  # ------------------------

  k <- c(2,1) # SWAP indices for eta[] and t[]
  time <- 0; dt <- 0; loop <- 0
  repeat {
    loop <- loop + 1
    k <- rev(k) # Swap indices
    if(HEATING_IS_TIME_DEPENDENT && time >= END_TIME/2) {
      q <- -q # make non-adiabatic heating time dependent
      HEATING_IS_TIME_DEPENDENT <- !HEATING_IS_TIME_DEPENDENT
    }
    
    if(time %in% change_dt_on) # You can have different dts for different time intervals.
      dt <- dts[which(time == change_dt_on, arr.ind = T)]
  
    # Routine JAC computes the horizantal advection of vorticity.
    # Routine LAPLAC computes the laplacian of the stream function.
  
    # Compute initial estimate of eta from vorticity equation. Equation (4)
    # FORWARD DIF. IN TIME (EXPLICIT)
    i<-2:(n-1); j<-2:(m-1)
    jac_psi_eta   <- JAC(psi, eta[,,k[1]], dx)[i,j]
    lap_eta       <- LAPLAC(eta[,,k[1]], dx)[i,j]
    eta[i,j,k[2]] <- eta[i,j,k[1]] + dt*(jac_psi_eta - g*d_dx(t[,j,k[1]],dx) + gnu*lap_eta)
    eta[,1,k[2]]  <- eta[,m,k[2]] <- 0 # Freeslip boundary conditions
    eta[1,,k[2]]  <- eta[n,,k[2]] <- 0
    
    # Calculate initial estimate of phi from the heat transfer equation. Equation (5)
    # FORWARD DIF. IN TIME (EXPLICIT)
    jac_psi_phi <- JAC(psi, t[,,k[1]], dx)
    lap_phi     <- LAPLAC(t[,,k[1]],dx)
    t[,,k[2]]   <- t[,,k[1]] + dt*(jac_psi_phi + q/tneut + gnu*lap_phi)
    
    # Relax eta to get psi.Routine RELAX1 solves the poisson equation. Equation (6)
    psi <- RELAX1(psi, eta[,,k[2]], n, m, alpha, dx)
    
    # Calculate final (corrector) estimates of eta and phi for this
    # time step from the predicted vorticity and temparature fields.
    # # Equation (4)
    i <- 2:(n-1); j <- 2:(m-1)
    jac_psi_eta   <- JAC(psi, eta[,,k[2]], dx)[i,j]
    lap_eta       <- LAPLAC(eta[,,k[2]],dx)[i,j]
    eta[i,j,k[2]] <- eta[i,j,k[1]] + dt*(jac_psi_eta - g*d_dx(t[,j,k[2]],dx) + gnu*lap_eta)
    
    # Relax final estimate of eta to get final psi
    # and consequently u and w fields. Equation (6)
    psi <- RELAX1(psi, eta[,,k[2]], n, m, alpha, dx)
    
    # Compute final estimate of phi for this time step. Equation (5)
    jac_psi_phi <- JAC(psi, t[,,k[2]], dx)
    lap_phi     <- LAPLAC(t[,,k[2]], dx)
    t[,,k[2]]   <- t[,,k[1]] + dt*(jac_psi_phi + q/tneut + gnu*lap_phi)
    
    # The horizontal velocity and the vertical velocity
    # are defined from the stream function.
    i <- 1:n;     j <- 1:(m-1); u[i,j] <-  (psi[i,j+1]-psi[i,j])/dx
    i <- 1:(n-1); j <- 1:m;     w[i,j] <- -(psi[i+1,j]-psi[i,j])/dx
    w[n,] <- w[(n-1),]; u[,m] <- u[,(m-1)] # set boundaries
    
    time <- time + dt
    res$time[loop+1]  <- time
    res$phi[,,loop+1] <- t[,,k[2]]
    res$eta[,,loop+1] <- eta[,,k[2]]
    res$psi[,,loop+1] <- psi
    res$t[,,loop+1]   <- t[,,k[2]]*tneut
    res$u[,,loop+1]   <- u
    res$w[,,loop+1]   <- w
  
    if(time >= END_TIME) break;
  }
  
  integ <- function(mat) apply(mat, 3, integral, dx)
  res$ws      <- sqrt((res$u)^2 + (res$w)^2)
  res$intg_ws <- integ(res$ws)
  res$sum_ws  <- colSums(res$ws, dims=2)
    
  res$k      <- 0.5*((res$u)^2 + (res$w)^2)
  res$intg_k <- integ(res$k)
  res$sum_k  <- colSums(res$k, dims=2)
  
  zmat <- array(matrix(res$z, nrow=n, ncol=m, byrow=T), c(n,m,l))
  
  res$k_gz_phi      <- (res$k - g*zmat*res$phi)
  res$intg_k_gz_phi <- integ(res$k_gz_phi)
  res$sum_k_gz_phi  <- colSums(res$k_gz_phi, dims=2)
  
  res$phi_square      <- (res$phi)^2
  res$intg_phi_square <- integ(res$phi_square)
  res$sum_phi_square  <- colSums(res$phi_square, dims=2)
  
  res$t_square      <- res$t^2
  res$intg_t_square <- integ(res$t_square)
  res$sum_t_square  <- colSums(res$t_square, dims=2)
  
  max_val <- function(mat) apply(mat, 3, max, na.rm=T)
  res$max_ws  <- max_val(res$ws)
  res$max_w   <- max_val(res$w)
  res$max_u   <- max_val(res$u)
  res$max_t   <- max_val(res$t)
  res$max_psi <- max_val(abs(res$psi))
  res$max_eta <- max_val(res$eta)
  
  min_val <- function(mat) apply(mat, 3, min, na.rm=T)
  res$min_ws  <- min_val(res$ws)
  res$min_w   <- min_val(res$w)
  res$min_u   <- min_val(res$u)
  res$min_t   <- min_val(res$t)
  res$min_psi <- min_val(abs(res$psi))
  res$min_eta <- min_val(res$eta)
  
  max_z <- function(mat) apply(mat, 3, function(t) {res$z[min(which(abs(t)==max(abs(t), na.rm=T), arr.ind=T)[,2])]})
  res$max_z_ws  <- max_z(res$ws)
  res$max_z_w   <- max_z(res$w)
  res$max_z_u   <- max_z(res$u)
  res$max_z_t   <- max_z(res$t)
  res$max_z_psi <- max_z(res$psi)
  res$max_z_eta <- max_z(res$eta)
  return(res)
}

run.nickerson <- function() {
  ptm  <- proc.time()
  runA <- nicker(reduce.gradient = T)
  runB <- nicker(Q0 = 0, gnu=0)
  runC <- nicker(Q0 = 0, gnu=0, reduce.gradient = T)
  # Krishnamurti run
  # runD <- nicker(dt=3, Q0=2*(4e-3)/pi, reduction.fact = 1, reduce.gradient=T, HEATING_IS_TIME_DEPENDENT=T)
  system("say Nickerson Run was completed")
  print(proc.time() - ptm)
  
  minutes <- runA$time/60
  nplot.energy.integrals(minutes, runA$intg_k, runA$intg_k_gz_phi, runC$intg_k_gz_phi, cols=c(2,2,4))
  nplot.phi_square(minutes, runA$intg_phi_square, runC$intg_phi_square)
  nplot.max.values(minutes, runA$max_t, runB$max_t, runC$max_t,ylab="theta - theta[0] (degree~C)", legend.pos = "topright")
  nplot.max.values(minutes, runA$max_eta, runB$max_eta, runC$max_eta, ylab="eta (s^-1)")
  nplot.max.values(minutes, runA$max_ws, runB$max_ws, runC$max_ws, ylab="\"Wind Speed\"~(m/sec)")
  nplot.max.values(minutes, runA$max_w, runA$max_u, ylab="\"Wind Speed Components\"~(m/sec)")
  nplot.max.values(minutes, runA$max_z_t, runA$max_z_w, runA$max_z_psi, ylab="\"z (meters)\"")
  
#   dx5dt0.5 <- list(runA=runA, runB=runB, runC=runC)
#   save(runA, runB, runC, file="dx5dt0.5.rdata")
  
}


