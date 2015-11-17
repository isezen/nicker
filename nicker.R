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
# Define variables.
#
# psi       : streamfunction
# difu      : eddy dissipation (m**2/sec)
#

source("numeric.R")
JAC <- JAC3

d_dx <- function(m, dx) {i<-2:(n-1);(m[i+1,]-m[i-1,])/(2*dx)}

# ----------------------- 
dts          <- c(3, 3, 3)                 # you can define different time intervals (seconds)
change_dt_on <- c(0, 360, 480)             # change dt value at this time stamp (seconds)
OUTPUT_TIMES <- c(0, 6, 12, 120, 360, 600) # Save outputs of the system at this time stamps
n <- 39; m <- 61 # vertical (n) and horizontal (m) dimensions
gnu   <- 0.5     # eddy kinematic coefficient of heat and momentum
dx    <- 10      # grid spacing. 10 meters in horizontal and vertical
tneut <- 300     # potential temperature of neutral environment
g     <- 9.81    # accerelation of gravity
alpha <- 1.89    # relaxation factor

result <- list(); index <- 1 # Holder for calculation results

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
              i <- 1:17  # The thermal having a maximum potential temperature excess of   
              j <- 11:31 # 0.5 deg celcius is inserted into a neutral environment initially rest.
     t[i, j, 1] <- outer(i, j, function(i,j) 0.5*cos(pi*(i-1)/32)*(cos(pi*((j-1)-10)/40))^2)
    t[i, 10, 1] <- t[i, 11, 1]  # decrease sharp gradiant at the bottom of buble.
         t[,,1] <- t[,,1]/tneut # Convert t to phi by dividing t by tneut (t is phi anymore.)

result[[index]] <- list(time=0, phi=t[,,1], t=t[,,1]*tneut, psi=psi, u=u, w=w) # Save initial situation

j <- 9:13 # Time invariant non-adiabatic heating (deg /sec) is defined as:
q[i,j] <- 2*outer(i, j, function(i,j) (4e-03)*cospi((i-1)/32)*(cospi(((j-1)-10)/4))^2)/pi
# ------------------------

time <- 0
k <- c(2,1) # SWAP indices for eta[] and t[]
dt <- 0
repeat {
  k <- rev(k) # Swap indices
  if(time == 306) q <- -q # make non-adiabatic heating time dependent
  
  if(time %in% change_dt_on) # You can have different dts for different time intervals.
    dt <- dts[which(time == change_dt_on, arr.ind = T)]

  # Compute initial estimate of eta from vorticity equation.
  # Routine JAC computes the horizantal advection of vorticity.
  # Routine LAPLAC computes the laplacian of the streamfunction.
  
  # Equation (4)
  i<-2:(n-1); j<-2:(m-1)
  jac_psi_eta <- JAC(psi, eta[,,k[1]], dx)[i,j]
  lap_eta     <- LAPLAC(eta[,,k[1]], dx)[i,j]
  eta[i,j,k[2]] <- eta[i,j,k[1]] + dt*(jac_psi_eta - g*d_dx(t[,j,k[1]],dx) + gnu*lap_eta)
  
  eta[,1,k[2]] <- eta[,m,k[2]] <- 0
  eta[1,,k[2]] <- eta[n,,k[2]] <- 0

  # Calculate initial estimate of phi from the heat transfer equation.
  # Equation (5)
  jac_psi_phi <- JAC(psi, t[,,k[1]], dx)
  lap_phi     <- LAPLAC(t[,,k[1]],dx)
  t[,,k[2]] <- t[,,k[1]] + dt * (jac_psi_phi + q/tneut + gnu*lap_phi)
  
  # Relax eta to get psi.Routine RELAX1 solves the poisson equation.
  # Equation (6)
  psi <- RELAX1(psi, eta[,,k[2]], n, m, alpha, dx)
  
  #
  # Calculate final (corrector) estimates of eta and phi for this
  # time step from the predicted vorticity and temparature fields.
  # # Equation (4)
  i <- 2:(n-1); j <- 2:(m-1)
  jac_psi_eta <- JAC(psi, eta[,,k[2]], dx)[i,j]
  lap_eta     <- LAPLAC(eta[,,k[2]],dx)[i,j]
  eta[i,j,k[2]] <- eta[i,j,k[1]] + dt*(jac_psi_eta - g*d_dx(t[,j,k[2]],dx) + gnu*lap_eta)
  
  #
  # Relax final estimate of eta to get final psi
  # and consequently u and w fields.
  # Equation (6)
  psi <- RELAX1(psi, eta[,,k[2]], n, m, alpha, dx)
  
  # Compute final estimate of phi for this time step
  # Equation (5)
  jac_psi_phi <- JAC(psi, t[,,k[2]], dx)
  lap_phi     <- LAPLAC(t[,,k[2]], dx)
  t[,,k[2]] <- t[,,k[1]] + dt*(jac_psi_phi + q/tneut + gnu*lap_phi)

  # The horizontal velocity and the vertical velocity
  # are defined from the stream function.
  i <- 1:n;     j <- 1:(m-1); u[i,j] <-  (psi[i,j+1]-psi[i,j])/dx
  i <- 1:(n-1); j <- 1:m;     w[i,j] <- -(psi[i+1,j]-psi[i,j])/dx
  w[n,] <- w[(n-1),]; u[,m] <- u[,(m-1)] # set boundaries
  
  time <- time + dt
  
  cat("time=", time, "\n")
  index <- index + 1
  result[[index]] <- list(time=time, phi=t[,,k[2]], t=t[,,k[2]]*tneut, psi=psi, u=u, w=w)

  if(time >= 600) break;
}

double_side <- F
x <- (0:(n-1))*dx
z <- (0:(m-1))*dx

len <- length(result)
tmax <- wmax <- psimax <- phi_square <- sum_k_gz_phi <- sum_k <- rep(0,len)
ztmax <- zwmax <- zpsimax <- rep(0,len)
time <- rep(0,len)
for(i in 1:len) {
  o <- result[[i]]
  time[i] <- o$time/60
  tmax[i]  <- max(o$t)
  ztmax[i] <- z[min(which(o$t==tmax[i], arr.ind = T)[,2])]

  wmax[i]  <- max(o$w)
  zwmax[i] <- z[min(which(o$w==wmax[i], arr.ind = T)[,2])]
  
  psimax[i]  <- min(o$psi)
  zpsimax[i] <- z[min(which(o$psi==psimax[i], arr.ind = T)[,2])]
  
  phi_square[i] <- sum((o$t/tneut)^2)*10000
  
  result[[i]]$k <- 0.5*(o$u^2 + o$w^2)
  sum_k[i] <- sum(result[[i]]$k)
  result[[i]]$k_gz_phi <- (result[[i]]$k - g*z*result[[i]]$phi)
  
  sum_k_gz_phi[i] <- sum(result[[i]]$k_gz_phi) 
}

plot(x=time, y=phi_square, type="l", xlab="t (minutes)", ylab="phi^2")
plot(x=time, y=sum_k, type="l", xlab="t (minutes)", ylab="sum_k")
plot(x=time, y=sum_k_gz_phi, type="l", xlab="t (minutes)", ylab="sum_k_gz_phi")

plot(x=time, y=ztmax, xlab="minutes", ylab="z (meters)", las=1, type="l", ylim=c(0,600))
lines(x=time, zwmax, col="red", lty=2)
lines(x=time, zpsimax, col="blue", lty=5)

for(i in 1:length(result)) {
  if(result[[i]]$time %in% OUTPUT_TIMES) {
    o <- result[[i]]
    time <- o$time
    psi <- o$psi; t <- o$t
    t[which(t<0)] <- 0 # remove negative values
    x2 <- x
    z2 <- z
    if(double_side) {
      psi <- rbind(-psi[nrow(psi):2,], psi)
      t   <- rbind(t[nrow(t):2,], t)
      x2  <- c(-x[length(x):2],x)
    }
    if(time==120){
      psi_levels <- c(40, 30, 20, 10, -10, -20, -30, -40)
    }else if(time==360){
      psi_levels <- c(120, 90, 60, 30, -30, -60, -90, -120)
    }else if(time==600){
      psi_levels <- c(180, 130, 90, 45, -45, -90, -130, -180)
    }else{
      zlim <- range(psi, finite=T)
      psi_levels <- pretty(zlim, 10)
      psi_levels <- psi_levels[-which(psi_levels==0, arr.ind = T)]
    }
    cat(time, ") ", max(t), "\n")
    # contour(x2, z2, t, lwd=2, lty=1)
    plot_time <- if(time<60) paste0("Time=",time, " sec.") else paste0("Time=",time/60, " min.")
    filled.contour(x2, z2, t, color = colorRampPalette(c("white", "blue", "green", "yellow","orange","red")), xlab="x", ylab="z",main=plot_time,
                   plot.axes={axis(1); axis(2); contour(x2, z2, psi, levels=psi_levels, lwd=2, lty=2, add = T)})
  }
}

stop()

require(animation)
saveGIF({
  for(i in 1:length(result)) {
    o    <- result[[i]]
    time <- o$time
    psi  <- o$psi; t <- o$t
    t[which(t<0)] <- 0 # remove negative values
    x2 <- x
    z2 <- z
    if(double_side) {
      psi <- rbind(-psi[nrow(psi):2,], psi)
      t   <- rbind(t[nrow(t):2,], t)
      x2  <- c(-x[length(x):2],x)
    }
    zlim <- range(psi, finite=T)
    psi_levels <- pretty(zlim, i)
    if(length(psi_levels)>1) psi_levels <- psi_levels[-which(psi_levels==0, arr.ind = T)]
    filled.contour(x2, z2, t, color =colorRampPalette(c("white", "blue", "green", "yellow","orange","red")), 
                   xlab="x", ylab="z",main=paste0("Time=",time, " sec."),
                   plot.axes={axis(1); axis(2); contour(x2, z2, psi, levels=psi_levels, lwd=1, lty=2, add = T)})
  }
},interval = 1)

