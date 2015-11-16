# c
# c  This program simulates the evolution of a simple
# c  cloud model as described by nickerson (1965).
# c  The heat source is inserted at the bottom of the
# c  atmosphere and the fluid is initially stratified
# c  and at rest.
# c  The thermal begins to ascent,and a vortex circulation
# c  develops . The temperature within the thermal
# c  increases due to heating and  the thermal rises above the
# c  heat source. The maximum temperature reaches a peak
# c  then decreases when dissipation exceed buoyancy.
# c  The potential temperature field resembles a thin
# c  stemmed mushroom.
# c
# c  Define variables.
# c
# c  n         : vertical dimension
# c  m         : horizontal dimension
# c  dt        : time step here dt=3 seconds
# c  loopt     : output interval equivalent to 2 minutes
# c  loop      : number of time steps. here max 200 = 10 minutes
# c  gnu       : eddy kinematic coefficient of heat and momentum
# c  dx        : grid spacing. 10 meters in horizontal and vertical
# c  tneut     : potential temperature of neutral environment
# c  alpha     : relaxation factor
# c  eta       : vorticity (per sec)
# c  t         : temperature excess (deg celcius)
# c  psi       : streamfunction
# c  difu      : eddy dissipation (m**2/sec)
# c

source("numeric.R")

JAC <- JAC2

# ----------------------- 
l <- n <- 39; m <- 61
n1    <- n-1
m1    <- m-1
dt    <- 1
loopt <- 2
LOOP_MAX <- 200
gnu   <- 0.5
dx    <- 10
tneut <- 300
g     <- 9.81
nold  <- 1
new   <- 2
alpha <- 1.89

result <- list()
index <- 0
x <- (0:(n-1))*dx
z <- (0:(m-1))*dx

d <- rep(dx, m)
# Initialize eta,phi(=t),psi,and q
psi <- q <- arak <- difu <- u <- w <- tmap <- wrk <- array(0,c(n,m))
eta <- t <- array(0,c(n, m, 2))

#  BASIC -----------------
#  This subroutine defines the initial state of the problem for                 
#  the stream function, vorticity, excess potential temperature,                
#  and the time invariant heating function.                                     
#  The boundary conditions used are those for freeslip, insulated
#  surfaces.                                                 
#  hort vel = 0. at x = 0 and x = 380                                            
#  vert vel = 0. at z = 0 and z = 600                                            
#  initially u = w = psi = eta  = 0.                                              
#  The temperature excess defined in the x-dir from 0 to 160m                        
#  and in the z-direction from 100m to 300 m. Everywhere else 
#  is set to zero .The heating field is defined in x-dir from 
#  0 to 160 m and in z-dir from 80 to 120 m .It is set to 0                                                
#  elsewhere.

# The thermal having a maximum potential temperature                               
# excess of 0.5 deg celcius is inserted into a                                  
# neutral environment initially rest.
i <- 1:17; j <- 11:31
     t[i, j, nold] <- outer(i, j, function(i,j) 0.5*cos(pi*(i-1)/32)*(cos(pi*((j-1)-10)/40))^2)
    t[i, 10, nold] <- t[i, 11, nold]
t[1:16,10:30,nold] <- t[1:16,10:30,nold]/tneut # Write output of initial state.
               wrk <- t[,,nold]
# The non-adiabatic heating (deg /sec) is defined as:
# Define the heating as a function of time
j=9:13
q[i,j] <- 2*outer(i, j, function(i,j) (4e-03)*cospi((i-1)/32)*(cospi(((j-1)-10)/4))^2)/pi
# ------------------------

time <- 0

index <- index + 1
result[[index]] <- list(time=time, wrk=wrk, psi=psi)

loop <- 0
repeat {
  if(loop == LOOP_MAX) break;
  loop <- loop + 1
  if(loop != 1) {
    nsave <- nold
    nold  <- new
    new   <- nsave
  }
  if(loop ==101) q <- -q
  
  #
  # Compute initial estimate of eta from vorticity equation.
  # Routine Jac computes the horizantal advection of vorticity.
  # Routine LAPLAC computes the laplacian of the streamfunction.
  #
  arak <- JAC(psi, eta[,,nold], dx)
  difu <- LAPLAC(eta[,,nold],dx)
  
  i=2:n1;j=2:m1
  eta[i,j,new] <- eta[i,j, nold] + dt* arak[i,j]+dt*gnu*difu[i,j]-dt*g*(t[i+1,j,nold]-t[i-1,j,nold])/(2*dx)
  
  i=1:n;    eta[i,1,new] <- eta[i,m,new] <- 0
  j <- 1:m; eta[1,j,new] <- eta[n,j,new] <- 0
  
  # Calculate initial estimate of phi from the heat transfer equation.
  arak <- JAC(psi, t[,,nold], dx)
  difu <- LAPLAC(t[,,nold],dx)
  
  i <- 1:n; j <- 1:m
  t[i,j,new] <- t[i,j,nold]+dt*(arak[i,j]+q[i,j]/tneut+gnu*difu[i,j])
  
  # Relax eta to get psi.Routine RELAX1 solves
  # the poisson equation.
  psi <- RELAX1(psi, eta[,,new], n, m, alpha, dx)
  
  #
  # The horizontal velocity and the vertical velocity
  # are defined from the streamfunction.
  #
  i <- 1:n; j <- 1:m1
  u[i,j] <- (psi[i,j+1]-psi[i,j])/dx
  
  i <- 1:n1; j <- 1:m
  w[i,j] <- -(psi[i+1,j]-psi[i,j])/dx
  
  #
  # Calculate final (corrector) estimates of eta and phi for this
  # time step from the predicted vorticity and temparature fields.
  #
  arak <- JAC(psi, eta[,,new], dx)
  difu <- LAPLAC(eta[,,new],dx)
  
  i <- 2:n1; j <- 2:m1
  eta[i,j,new] <- eta[i,j,nold]+dt*(arak[i,j]+gnu*difu[i,j]-g*(t[i+1,j,new]-t[i-1,j,new])/(2*dx))
  
  #
  # Relax final estimate of eta to get final psi
  # and consequently u and w fields.
  #
  psi <- RELAX1(psi, eta[,,new], n, m, alpha, dx)
  
  i <- 1:n; j <- 1:m1
  u[i,j] <- (psi[i,j+1]-psi[i,j])/dx
  
  i <- 1:n1; j <- 1:m
  w[i,j] <- -(psi[i+1,j]-psi[i,j])/dx
  
  #
  # Compute final estimate of phi for this time step
  #
  arak <- JAC(psi, t[,,new], dx)
  difu <- LAPLAC(t[,,new],dx)
  
  i <- 1:n; j <- 1:m
  t[i,j,new] <- t[i,j,nold]+dt*(arak[i,j]+q[i,j]/tneut+gnu*difu[i,j])
  
  #
  # Check if output is required.
  # loopt is the number time steps equivalent to 2 minutes
  # for output interval. the outputs can be displayed on the
  # screen by calling cloud with the fourth argument being 1.
  #
  time <- time + dt
  if((loop %% loopt) != 0) next
  #  cat("nscan=", nscan, "\n")
  time <- dt*loop
  cat("time=", time, "\n")
  #  call CLOUD (t(1,1,new),n,m,0)
  
  i <- 1:n; j <- 1:m
  tmap[i,j] <- t[i,j,new]*tneut
  
  shditv <- 0.1
  
  i <- 1:l; j <- 1:m
  wrk[i,j] <- tmap[i,j]
  
  #
  # Write output of temperature anomaly
  # for each required time.
  #
  #   print(time)
  #   print(wrk)
  index <- index + 1
  result[[index]] <- list(time=time, wrk=wrk, psi=psi)
  
  #  call CLOUD (psi,n,m,0)
  shditv <- time/60
  #
  # Write output of streamfunctions
  # for each required time.
  #
  
  #   write (20,119)time
  #   write (20,120)psi
  
  # call CLOUD (eta(1,1,new),n,m,0)
  shditv <- time/30000
  #  write (6,20)time
  w[n,1:m] <- w[n1,1:m]
  
  # call CLOUD (w,n,m,0)
  shditv <- time/1200
  # write (6,21)time
  u[1:n,m] <- u[1:n,m1]
  # call CLOUD (u,n,m,0)
}

stop()


require(animation)
saveGIF({
  for(i in 2:length(result)) {
    time <- result[[i]]$time
    psi <- result[[i]]$psi
    psi <- rbind(psi[nrow(psi):2,],psi)
    wrk <- abs(result[[i]]$wrk)
    wrk[which(wrk<1e-03)] <- 0
    wrk <- rbind(wrk[nrow(wrk):2,],wrk)
    x2 <- c(-x[length(x):2],x)
    z2 <- z
    filled.contour(x2, z2, wrk, color = heat.colors, xlab="x", ylab="z",main=paste0("Time=",time),
                   plot.axes={axis(1); axis(2); contour(x2, z2, psi, lwd=2, lty=1, add = T)})
  }
},interval = 0.1)


