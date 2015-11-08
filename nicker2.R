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


# This subroutine prints out the results.
# index = 1 for printing
# index = 0 no  printing
CLOUD <- function(rt, n, m, index) {
  rt    <- matrix(0, nrow=n, ncl=m)
  rag   <- matrix(0, nrow=39, ncol=61)
  const <- 1/10000000
  for(i in 1:n){
    for(j in 1:m){
      rag[i,j] <- rt[i,j]
    }
  }
  for(k in seq(1,n,10)){
    ik <- k
    k1 <-k+9
    if(k1 > n) k1 <- n
    if(index == 1) print()
  }
}


# ----------------------- 
l <- n <- 39; m <- 61
n1    <- n-1
m1    <- m-1
n2    <- n-2
m2    <- m-2
dt    <- 3
loopt <- 40
gnu   <- 0.5
dx    <- 10
e     <- dx
tneut <- 300
g     <- 9.81
nold  <- 1
new   <- 2
alpha <- 1.89

result <- list()
index <- 0

d <- rep(dx, m)
# Initialize eta,phi(=t),psi,and q
psi <- q <- arak <- difu <- u <- w <- tmap <- wrk <- array(0,c(n,m))
eta <- t <- array(0,c(n, m, 2))

#  BASIC
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
t[i, 10, nold]  <- t[i, 11, nold]
# The non-adiabatic heating (deg /sec) is defined as:
j <- 9:13
q[i, j] <- outer(i, j, function(i,j) 4e-03*cos(pi*(i-1)/32)*(cos(pi*((j-1)-10)/4))^2)

time <- 0

wrk <- t[,,nold]

index <- index + 1
result[[index]] <- list(time=time, wrk=wrk, psi=psi)

# Write output of initial state.
t[1:16,10:30,nold] <- t[1:16,10:30,nold]/tneut

loop <- 0

repeat {
  if(loop == 200) break;
  loop <- loop + 1
  if(loop != 1) {
    nsave <- nold
    nold  <- new
    new   <- nsave
  }
  # Define the heating as a function of time
  i=1:17;j=9:13
  q[i,j] <- 2*outer(i, j, function(i,j) (4e-03)*cos(pi*(i-1)/32)*(cos(pi*((j-1)-10)/4))^2)/pi
  if(loop >100 && loop <=200) q[i,j] <- -q[i,j]
  
  #
  # Compute initial estimate of eta from vorticity equation.
  # Routine Jac computes the horizantal advection of vorticity.
  # Routine LAPLAC computes the laplacian of the streamfunction.
  #
  arak <- JAC2(psi, eta[,,nold], dx)
  difu <- LAPLAC(eta[,,nold],dx)
  
  i=2:n1;j=2:m1
  eta[i,j,new] <- eta[i,j, nold] + dt* arak[i,j]+dt*gnu*difu[i,j]-dt*g*(t[i+1,j,nold]-t[i-1,j,nold])/(2*dx)
  
  i=1:n;    eta[i,1,new] <- eta[i,m,new] <- 0
  j <- 1:m; eta[1,j,new] <- eta[n,j,new] <- 0
  
  # Calculate initial estimate of phi from the
  # heat transfer equation.
  arak <- JAC2(psi, t[,,nold], dx)
  difu <- LAPLAC(t[,,nold],dx)
  
  i <- 1:n; j <- 1:m
  t[i,j,new] <- t[i,j,nold]+dt*arak[i,j]+dt*q[i,j]/tneut+dt*gnu*difu[i,j]
  
  # Relax eta to get psi.Routine RELAX1 solves
  # the poisson equation.
  psi <- RELAX1(psi, eta[,,new], n, m, nscan, alpha)
  
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
  arak <- JAC2(psi, eta[,,new], dx)
  difu <- LAPLAC(eta[,,new],dx)
  
  for(i in 2:n1)
    for(j in 2:m1)
      eta[i,j,new] <- eta[i,j,nold]+dt*arak[i,j]+dt*gnu*difu[i,j]-dt*g*(t[i+1,j,new]-t[i-1,j,new])/(2*dx)
  
  #
  # Relax final estimate of eta to get final psi
  # and consequently u and w fields.
  #
  psi <- RELAX1(psi, eta[,,new], n, m, nscan, alpha)
  
  i <- 1:n; j <- 1:m1
  u[i,j] <- (psi[i,j+1]-psi[i,j])/dx
  
  i <- 1:n1; j <- 1:m
  w[i,j] <- -(psi[i+1,j]-psi[i,j])/dx
  
  #
  # Compute final estimate of phi for this time step
  #
  arak <- JAC2(psi, t[,,new], dx)
  difu <- LAPLAC(t[,,new],dx)
  
  i <- 1:n; j <- 1:m
  t[i,j,new] <- t[i,j,nold]+dt*arak[i,j]+dt*q[i,j]/tneut+dt*gnu*difu[i,j]
  
  #
  # Check if output is required.
  # loopt is the number time steps equivalent to 2 minutes
  # for output interval. the outputs can be displayed on the
  # screen by calling cloud with the fourth argument being 1.
  #
  time <- time + dt
  if((loop %% loopt) != 0) next
  #  cat("nscan=", nscan, "\n")
  time <- 3*loop
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
  if(loop == 300) break
}

x <- (0:(n-1))*dx
z <- (0:(m-1))*dx
for(i in 1:6) {
  filled.contour(x, z, result[[i]]$psi, color = heat.colors, xlab="x", ylab="z",main=paste0("Time=",result[[i]]$time),
                 plot.axes={axis(1); axis(2); contour(x, z, result[[i]]$wrk, nlevels=5, lwd=2, lty=1, add = T)})
}

