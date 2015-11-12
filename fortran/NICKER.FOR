      program NICKER
c
c  This program simulates the evolution of a simple
c  cloud model as described by nickerson (1965).
c  The heat source is inserted at the bottom of the
c  atmosphere and the fluid is initially stratified
c  and at rest.
c  The thermal begins to ascent,and a vortex circulation
c  develops . The temperature within the thermal
c  increases due to heating and  the thermal rises above the
c  heat source. The maximum temperature reaches a peak
c  then decreases when dissipation exceed buoyancy.
c  The potential temperature field resembles a thin
c  stemmed mushroom.
c
c  Define variables.
c
c  n         : vertical dimension
c  m         : horizontal dimension
c  dt        : time step here dt=3 seconds
c  loopt     : output interval equivalent to 2 minutes
c  loop      : number of time steps. here max 200 = 10 minutes
c  gnu       : eddy kinematic coefficient of heat and momentum
c  dx        : grid spacing. 10 meters in horizontal and vertical
c  tneut     : potential temperature of neutral environment
c  alpha     : relaxation factor
c  eta       : vorticity (per sec)
c  t         : temperature excess (deg celcius)
c  psi       : streamfunction
c  difu      : eddy dissipation (m**2/sec)
c
      parameter (l = 39, m = 61,n = 39)
      common d(61),tmap(39,61),wrk(39,61)
      common psi(39,61),eta(39,61,2),t(39,61,2),q(39,61),
     &       arak(39,61),difu(39,61),u(39,61),w(39,61)
c
      open(20,file='nicker.out',status='unknown')
c
   16 format(4x,'the number of scans required to obtained the',
     & 1x,'following stream function was ',i3,' in this time step')
   17 format(4x,'excess potential temperature at
     & time= ',f5.1,' seconds')
   18 format(4x,'stream function',10x,'elasped time= ',f5.1,
     & 'seconds')
   19 format(4x,'vorticity',15x,'time=',f5.1)
   20 format(4x,'vertical velocity',15x,'time',f5.1)
   21 format(4x,'horizontal velocity',15x,'time=',f5.1)
   25 format(4x,'the number of scans required for the first'
     & ,'relaxation in this time step was ',i3)
 1001 format(/)
      write(6,1001)
c
c  Define  constants
c
      n1        = n-1
      m1        = m-1
      n2        = n-2
      m2        = m-2
      dt        = 3.
      loopt     = 40
      gnu       = 0.5
      dx        = 10.
      e         = dx
      tneut     = 300.
      g         = 9.81
      pi        = 4.0*atan(1.0)
      nold      = 1
      new       = 2
      alpha     = 1.89

      do 6600 i = 1, m
 6600    d(i)   = dx
c
c Initialize eta,phi(=t),psi,and q
c
      call BASIC (eta(1,1,nold),t(1,1,nold),psi,
     &                     u,w,q,n,m,n1,m1,tneut)
      time      = 0.

      do 6602 i = 1, l
      do 6602 j = 1, m
 6602   wrk(i,j)= t(i,j,nold)
c
c  Write output of initial state .
c
      write (20,119) time
      write (20,120) wrk
  119 format(f12.2)
  120 format(6e13.6)

      do 6604 i = 1, 16
      do 6604 j = 10, 30
      t(i,j,nold) = t(i,j,nold)/tneut
 6604 continue

      loop      = 0
c
c  Begin new time loop
c
  999 continue
      if (loop.eq.200) go to 35
      loop      = loop+1
      if (loop.eq.1) go to 1
      nsave     = nold
      nold      = new
      new       = nsave
      do 6606 j = 1, m
      do 6606 i = 1, n
      t(i,j,nold) = t(i,j,nold)
 6606 continue
    1 continue
c
c  Define the heating as a function of time
c
      do 6608 j = 9,13
      do 6608 i = 1,17
         qq     = 4.e-03*cos(pi*(i-1)/32.)*
     &          (cos(pi*((j-1)-10.)/4.))**2
      q(i,j)    = 2.0*qq/pi
      if (loop.gt.100.and.loop.le.200) q(i,j) = -2.0*qq/pi
 6608 continue
c
c  Compute initial estimate of eta from vorticity equation .
c  Routine Jac computes the horizantal advection of vorticity.
c  Routine LAPLAC computes the laplacian of the streamfunction .
c
      call JAC    (arak,psi,eta(1,1,nold),d,e,n,m,n1,m1,n2,m2)
      call LAPLAC (eta(1,1,nold),difu,dx,n,m,n1,m1)

      do 6610 i = 2, n1
      do 6610 j = 2, m1
         eta(i,j,new) = eta(i,j,nold)+dt*arak(i,j)+dt*gnu*difu(i,j)
     &                  -dt*g*(t(i+1,j,nold)-t(i-1,j,nold))/(2.*dx)
 6610 continue

      do 6612 i = 1, n
         eta(i,1,new) = 0.
         eta(i,m,new) = 0.
 6612 continue
      do 6614 j = 1, m
         eta(1,j,new) = 0.
         eta(n,j,new) = 0.
 6614 continue
c
c  Calculate initial estimate of phi from the
c  heat transfer equation .
c
      call JAC    (arak,psi,t(1,1,nold),d,e,n,m,n1,m1,n2,m2)
      call LAPLAC (t(1,1,nold),difu,dx,n,m,n1,m1)
      do 6616 i = 1, n
      do 6616 j = 1, m
         t(i,j,new) = t(i,j,nold)+dt*arak(i,j)+dt*q(i,j)
     &                           /tneut+dt*gnu*difu(i,j)
 6616 continue
c
c  Relax eta to get psi.Routine RELAX1 solves
c  the poisson equation .
c
      call RELAX1 (psi,eta(1,1,new),n,m,nscan,alpha)
c
c  The horizontal velocity and the vertical velocity
c  are defined from the streamfunction.
c
      do 6618 i = 1, n
      do 6618 j = 1, m1
        u(i,j)  = (psi(i,j+1)-psi(i,j))/dx
 6618 continue
      do 6620 j = 1, m
      do 6620 i = 1, n1
         w(i,j) = -(psi(i+1,j)-psi(i,j))/dx
 6620 continue
c
c  Calculate final (corrector) estimates of eta and phi for this
c  time step from the predicted vorticity and temparature fields.
c
      call JAC    (arak,psi,eta(1,1,new),d,e,n,m,n1,m1,n2,m2)
      call LAPLAC (eta(1,1,new),difu,dx,n,m,n1,m1)
      do 6622 i = 2, n1
      do 6622 j = 2, m1
      eta(i,j,new) = eta(i,j,nold)+dt*arak(i,j)+dt*gnu*difu(i,j)
     &                 -dt*g*(t(i+1,j,new)-t(i-1,j,new))/(2.*dx)
 6622 continue
c
c  Relax final estimate of eta to get final psi
c  and consequently u and w fields.
c
      call RELAX1 (psi,eta(1,1,new),n,m,nscan,alpha)
      do 6624 i = 1, n
      do 6626 j = 1, m1
         u(i,j) = (psi(i,j+1)-psi(i,j))/dx
 6626 continue
 6624 continue
      do 6628 j = 1, m
      do 6628 i = 1, n1
         w(i,j) = -(psi(i+1,j)-psi(i,j))/dx
 6628 continue
c
c  Compute final estimate of phi for this time step
c
      call JAC    (arak,psi,t(1,1,new),d,e,n,m,n1,m1,n2,m2)
      call LAPLAC (t(1,1,new),difu,dx,n,m,n1,m1)
      do 6630 i = 1, n
      do 6630 j = 1, m
         t(i,j,new) = t(i,j,nold)+dt*arak(i,j)
     &       +dt*q(i,j)/tneut+dt*gnu*difu(i,j)
 6630 continue
c
c  Check if output is required.
c  loopt is the number time steps equivalent to 2 minutes
c  for output interval. the outputs can be displayed on the
c  screen by calling cloud with the fourth argument being 1.
c
      time      = time+dt
c      if (mod(loop,loopt).ne.0) go to 999
      write (6,16)nscan
      time      = 3.0*loop
      write (6,17)time
      call CLOUD (t(1,1,new),n,m,0)
      do 6632 i = 1, n
      do 6632 j = 1, m
      tmap(i,j) = t(i,j,new)*tneut
 6632 continue
      shditv    = 0.1

      do 6634 j = 1, m
      do 6634 i = 1, l
 6634  wrk(i,j) = tmap(i,j)
c
c  Write output of temperature anomaly
c  for each required time .
c
      write (20,119)time
      write (20,120)wrk
c
      write (6,18)time
      call CLOUD (psi,n,m,0)
      shditv    = time/60.
c
c  Write output of streamfunctions
c  for each required time .
c
      write (20,119)time
      write (20,120)psi

      write (6,19)time
      call CLOUD (eta(1,1,new),n,m,0)
      shditv    = time/30000.
      write (6,20)time
      do 6636 j = 1, m
         w(n,j) = w(n1,j)
 6636 continue
      call CLOUD (w,n,m,0)
      shditv    = time/1200.
      write (6,21)time
      do 6638 i = 1, n
         u(i,m) = u(i,m1)
 6638 continue
      call CLOUD (u,n,m,0)
      if (loop.eq.300) go to 35
      go to 999
   35 continue
      stop
      end

      subroutine BASIC (eta,t,psi,u,w,q,n,m,n1,m1,tneut)
c
c  This subroutine defines the initial state of the problem for
c  the stream function, vorticity, excess potential temperature,
c  and the time invariant heating function.
c  The boundary conditions used are those for freeslip, insulated
c  surfaces.
c  hort vel = 0. at x = 0 and x = 380
c  vert vel = 0. at z = 0 and z = 600
c  initially u = w = psi = eta  = 0.
c  The temperature excess defined in the x-dir from 0 to 160m
c  and in the z-direction from 100m to 300 m. Everywhere else
c  is set to zero .The heating field is defined in x-dir from
c  0 to 160 m and in z-dir from 80 to 120 m .It is set to 0
c  elsewhere.
c
      real eta(n,m),t(n,m),psi(n,m),q(n,m),u(n,m),w(n,m)
      do 6610 i = 1, n
      do 6610 j = 1, m
       psi(i,j) = 0.
       eta(i,j) = 0.
         t(i,j) = 0.
         q(i,j) = 0.
         u(i,j) = 0.
         w(i,j) = 0.
 6610 continue
      pi        = 4.0*atan(1.)
c
c  The thermal having a maximum potential temperature
c  excess of 0.5 deg celcius is inserted into a
c  neutral environment initially rest.
c
      do 6611 i = 1, 17
      do 6611 j = 11, 31
         t(i,j) = 0.5*cos(pi*(i-1)/32.)*(cos
     &              (pi*((j-1)-10.)/40.))**2
 6611 continue
      do 6612 i = 1, 17
         t(i,10)= t(i,11)
 6612 continue
c
c  The non-adiabatic heating (deg /sec) is defined as:
c
      do 6613 i = 1, 17
      do 6613 j = 9, 13
         q(i,j) = 4.e-03*cos(pi*(i-1)/32.)*
     &          (cos(pi*((j-1)-10.)/4.))**2
 6613 continue
      return
      end
      subroutine CLOUD (rt,n,m,index)
c
c  This subroutine prints out the results.
c  index = 1 for printing
c  index = 0 no  printing
c
      dimension rt (n,m),rag (39,61)
   20 format(5x,10i11)
   21 format(i10,10e11.3)
      const=1./10000000.
      do 6630 i = 1, n
      do 6630 j = 1, m
        rag(i,j)= rt(i,j)
 6630 continue
      do 6631 k = 1, n, 10
         ik     = k
         k1     = k+9
        if (k1.gt.n) k1 = n
        if (index.eq.1) print 20,(i,i=k,k1)
      do 6632 j = 1, m
         j1     = 62-j
      do 6633 i = ik,k1
         if (abs(rag(i,j1)).lt.const) rag(i,j1)=0.
 6633 continue
      if (index.eq.1) print 21,j1,(rag(i,j1),i=ik,k1)
 6632 continue
 6631 continue
      return
      end
      subroutine JAC (a,b,c,dx,dy,l,m,l1,m1,l2,m2)
c
c  This subroutine performs the so-called Arakawa
c  Jacobian .This scheme satisfies the integral
c  relations of invariance for the total kinetic
c  energy and the mean square vorticity.
c  The following code uses a nine-point
c  second order Arakawa Jacobian.
c
c variable definitions
c
c  input  : b and c
c  output : a
c  dx     : east - west grid distance.
c  dy     : north-south grid distance.
c  l,m    : first and second dimension of a,b,c
c

      real a(l,m,1), b(l,m,1), dx(m), c(l,m,1)
      k         = 1
      do 2310 j = 2, m1
      dm        = 12.*dy*dx(j)
      do 2310 i = 1, l
      if (i-1) 80,80,81
   80 im1       = l1
      ip1       = 2
      go to 83
   81 if (i-l) 82,80,80
   82 im1       = i-1
      ip1       = i+1
   83 continue
c
      a(i,j,k)  = (b(i,j-1,k)+b(ip1,j-1,k)-b(i,j+1,k)-b(ip1,j+1,k))
     &  *(c(ip1,j,k)-c(i,j,k))+(b(im1,j-1,k)+b(i,j-1,k)-b(im1,j+1,k
     &  )-b(i,j+1,k))*(c(i,j,k)-c(im1,j,k))+(b(ip1,j,k)+b(ip1,j+1,k
     &  )-b(im1,j,k)-b(im1,j+1,k))*(c(i,j+1,k)-c(i,j,k))+(b(ip1,j-1
     &  ,k)+b(ip1,j,k)-b(im1,j-1,k)-b(im1,j,k))*(c(i,j,k)-c(i,j-1,k
     &  ))+(b(ip1,j,k)-b(i,j+1,k))*(c(ip1,j+1,k)-c(i,j,k))+(b(i,j-1
     &  ,k)-b(im1,j,k))*(c(i,j,k)-c(im1,j-1,k))+(b(i,j+1,k)-b(im1,j
     &  ,k))*(c(im1,j+1,k)-c(i,j,k))+(b(ip1,j,k)-b(i,j-1,k))*(c(i,j
     &  ,k)-c(ip1,j-1,k))
c
      a(i,j,k)  = a(i,j,k)/dm
 2310 continue
      do 2311 i = 1, l
      if (i-1) 70,70,71
   70 im1       = l1
      ip1       = 2
      go to 73
   71 if (i-l) 72,70,70
   72 im1       = i-1
      ip1       = i+1
   73 continue
      dm        = 12.*dy*dx(1)
c
      a(i,1,k)  = (b(i,1,k)+b(ip1,1,k)-b(i,2,k)-b(ip1,2,k))*(c(i,1,k)+
     &    c(ip1,1,k))-(b(im1,1,k)+b(i,1,k)-b(im1,2,k)-b(i,2,k))*(c(im1
     &    ,1,k)+c(i,1,k))+(b(ip1,1,k)+b(ip1,2,k)-b(im1,1,k)-b(im1,2,k)
     &    )*(c(i,1,k)+c(i,2,k))+(b(ip1,1,k)-b(i,2,k))*(c(i,1,k)+c(ip1,
     &    2,k))+(b(i,2,k)-b(im1,1,k))*(c(im1,2,k)+c(i,1,k))
c
      a(i,1,k)  = a(i,1,k)/dm
      dm        = 12.*dy*dx(m)
c
      a(i,m,k)  = (b(i,m-1,k)+b(ip1,m-1,k)-b(i,m,k)-b(ip1,m,k))*(c(i,m,
     &    k)+c(ip1,m,k))-(b(im1,m-1,k)+b(i,m-1,k)-b(im1,m,k)-b(i,m,k))*
     &    (c(im1,m,k)+c(i,m,k))-(b(ip1,m-1,k)+b(ip1,m,k)-b(im1,m-1,k)-b
     &    (im1,m,k))*(c(i,m-1,k)+c(i,m,k))-(b(i,m-1,k)-b(im1,m,k))*(c(i
     &    m1,m-1,k)+c(i,m,k))-(b(ip1,m,k)-b(i,m-1,k))*(c(i,m,k)+c(ip1,m
     &    -1,k))
c
      a(i,m,k)  = a(i,m,k)/dm
 2311 continue
      return
      end
      subroutine LAPLAC (t,difu,dx,n,m,n1,m1)
c
c  This subroutine calculates the Laplacian operator
c  of array t.It uses the 5 point second order scheme
c  as defined in Chapt.2 but has special treatment of
c  boundaries.
c
      real t(n,m), difu(n,m)
      da        = dx*dx
c
c  interior points
c
      do 6640 i = 2, n1
      do 6640 j = 2, m1
      difu(i,j) = (t(i,j+1)+t(i,j-1)+t(i+1,j)+t(i-1,j)-4.*t(i,j))/da
 6640 continue
c
c  side boundaries
c
      do 6641 j = 2, m1
      difu(1,j) = (2.*t(2,j)+t(1,j+1)+t(1,j-1)-4.*t(1,j))/da
      difu(n,j) = (2.*t(n1,j)+t(n,j+1)+t(n,j-1)-4.*t(n,j))/da
 6641 continue
c
c  top and bottom boundaries
c
      do 6642 i = 2, n1
      difu(i,1) = (2.*t(i,2)+t(i-1,1)+t(i+1,1)-4.*t(i,1))/da
      difu(i,m) = (2.*t(i,m1)+t(i-1,m)+t(i+1,m)-4.*t(i,m))/da
 6642 continue
      db        = 2./da
c
c  four corners
c
      difu(1,1) = db*(t(1,2)+t(2,1)-2.*t(1,1))
      difu(1,m) = db*(t(1,m1)+t(2,m)-2.*t(1,m))
      difu(n,m) = db*(t(n1,m)+t(n,m1)-2.*t(n,m))
      difu(n,1) = db*(t(n1,1)+t(n,2)-2.*t(n,1))
      return
      end
      subroutine RELAX1(psi,eta,n,m,nscan,alpha)
c
c  This subroutine solves the poisson equation to
c  obtain the streamfunctions .
c  eps is the tolerance error in streamfunction(e-03)
c  nscan is the maximum iteration (300) and
c  r is the residual value .
c
      real psi(n,m),eta(n,m)
      eps       = 1.e-03
      dx        = 10.
      nscan     = 0
      n1        = n-1
      m1        = m-1
   89 continue
      nscan     = nscan+1
      rmax      = 0.
      do 6620 i = 2, n1
      do 6620 j = 2, m1
         r1     = .25*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))
         r2     = -psi(i,j)-.25*dx**2*eta(i,j)
         r      = r1+r2
         if (rmax.lt.abs(r)) go to 90
         go to 91
   90    isave  = i
         jsave  = j
         rmax   = abs(r)
   91 continue
      psi(i,j)  = psi(i,j)+alpha*r
 6620 continue
      if (nscan.eq.300) print 95
   95 format(1x,'no convergence achieved')
      if(nscan.eq.300)stop
      ps        = abs(psi(isave,jsave))
      if (rmax/ps-eps) 83,89,89
   83 continue
      return
      end
