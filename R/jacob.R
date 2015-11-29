# JACOBIAN TEST FUNTIONS
# sezenismail@gmail.com 2015-11-07
#
#

source("numeric.R")

create_test_data <- function(xrange=c(0,1800), yrange=c(0,3800), dx=200, dy=dx) {
  dx <- dy <- h
  x <- seq(min(xrange), max(xrange), by=h)
  y <- seq(min(yrange), max(yrange), by=h)
  l<-length(x); m<-length(y)
  
  yl   <-    pi/1000
  yk   <-  2*yl
  ykmx <-  matrix(yk*x, nrow=l, ncol=m)
  ylmy <-  matrix(yl*y, nrow=l, ncol=m, byrow = T)
  psi  <-  sin(ykmx) * sin(ylmy) + cos(ylmy)
  zta  <- -(yk^2+yl^2)*sin(ykmx)*sin(ylmy) - (yl^2)*cos(ylmy)
  #
  delpsi_delx <-  yk *  cos(ykmx)   * sin(ylmy)
  delpsi_dely <-  yl * (sin(ykmx)   * cos(ylmy) - sin(ylmy))
  delzta_delx <- -yk * (yk^2+yl^2)  * cos(ykmx) * sin(ylmy)
  delzta_dely <- -yl * ((yk^2+yl^2) * sin(ykmx) * cos(ylmy) + (yl^2) * sin(ylmy))
  jj_analytic <- delpsi_delx*delzta_dely - delpsi_dely*delzta_delx
  return(list(x=x, y=y, psi=psi, zta=zta, jj_analytic=jj_analytic))
}

# This function has the same algorithm with JACOBIAN.FOR.
JAC1_TEST <- function(h=200) {
  x_min <- 0; x_max <- 1800
  y_min <- 0; y_max <- 3800
  dy <- h
  
  l <- ((x_max-x_min)/h)+1; m <- ((y_max-y_min)/h)+1
  l1<-l-1; m1<-m-1
  l2<-l-2; m2<-m-2
  a <- zta <- psi <- matrix(0, nrow=l, ncol=m)
  x <- rep(0, l)
  dx <- y <- rep(0, m)
  
  pi <- 4*atan(1)
  for(j in 1:m) dx[j] <- dy
  yk <- 2*pi/1000
  yl <-   pi/1000
  
  x[1] <- y[1] <- 0
  for(i in 2:l) {
    im1 <- i-1
    for(j in 2:m) {
      jm1 <- j-1
      x[i] <- x[im1] + h
      y[j] <- y[jm1] + h
    }
  }
  
  sum <- 0
  for(i in 1:l) {
    for(j in 1:m) {
      psi[i,j] <- sin(yk*x[i]) * sin(yl*y[j]) + cos(yl*y[j])
      zta[i,j] <- -(yk^2+yl^2)*sin(yk*x[i])*sin(yl*y[j])-(yl^2) *cos(yl*y[j])
      a[1,j] <- zta[1,j]
      a[l,j] <- zta[l,j]
      a[i,1] <- zta[i,1]
      a[i,m] <- zta[i,m]
      sum <- sum + (zta[i,j]/(l*m))^2
    }
  }
  return(list(psi=psi, zta=zta, JJ=JAC1(a, psi, zta, dx, dy, l, m, l1, m1, l2, m2)))
}

JAC2_TEST <- function(h=200) {
  data <- create_test_data(dx=h)
  #
  # NUMERIC SOLUTION OF JACOBIAN
  J1 <- JAC2(data$psi, data$zta, h)
  J2 <- JAC3(data$psi, data$zta, h)
  # filled.contour(x, y, psi, xlab="x", ylab="z", nlevels=20)
  
  #   filled.contour(x, y, psi, color = terrain.colors, xlab="x", ylab="y",
  #                  plot.axes={axis(1); axis(2); contour(x, y, zta, nlevels=20, color=heat.colors, lwd=2, lty=2, add = T)})
  #   
  #   nlevels=20
  #   zlim <- range(a2, finite=T)
  #   cat("h=", h, "psi_range=", range(psi, finite=T), "\n")
  #   cat("h=", h, "zta_range=", range(zta, finite=T), "\n")
  #   cat("h=", h, "jj_analytic_range=", range(jj_analytic, finite=T), "\n")
  #   cat("h=", h, "zlim=", zlim, "\n")
  #   levels <- pretty(zlim, nlevels)
  #   print(levels)
  #   filled.contour(x, y, jj_analytic, color = heat.colors, xlab="x", ylab="y",
  #                  plot.axes={axis(1); axis(2); contour(x, y, a2, nlevels=20, color=heat.colors, lwd=2, lty=2, add = T)})
  
  indx <- min(dim(data$jj_analytic))
  s1 <- data$jj_analytic[indx,]*1e10
  s2 <- J1[indx,]*1e10
  s3 <- J2[indx,]*1e10
  ylim <- c(s1,s2,s3)
  ylim <- c(min(ylim), max(ylim))
  plot(s1, col="red", type="l",ylim=ylim);points(s1, pch=20, col="red")
  lines(s2);points(s2, pch=20)
  lines(s3,col="green");points(s3, pch=20, col="green")
  cbind(s1,s2)
  
  colMax <- function(data) apply(data, 1, max, na.rm = TRUE)
  
  colmax1 <- colMax(abs(data$jj_analytic-J1))
  colmax2 <- colMax(abs(data$jj_analytic-J2))
  lim <- range(c(colmax1, colmax2))
  plot(colmax1, type="l", ylim=lim); points(colmax1, pch=20)
  lines(colmax2, type="l",col="red"); points(colmax2, pch=20, col="red")
  
  return(list(psi=data$psi,zta=data$zta,jj_analytic=data$jj_analytic, jj_numeric=J1))
}

h <- 200
a1 <- JAC1_TEST(h)
a2 <- JAC2_TEST(h)
cat("a1$psi == a2$psi ? ", all(a1$psi == a2$psi), "\n")
cat("a1$zta == a2$zta ? ", all(a1$zta == a2$zta), "\n")
cat("a1$JJ == a2$JJ ? ", all(a1$JJ == a2$JJ), "\n") 