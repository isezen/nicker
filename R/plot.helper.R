# This is a helper file to plot results.

nplot.contour <- function(x=seq(0, 1, length.out = nrow(z1)), 
                         y=seq(0, 1, length.out = ncol(z1)), t=0, z1, z2, 
                         z2.nlevels=20,
                         z1.lim=range(z1, finite=T), xlab="x (meters)", ylab="z (meters)", tlab="Time= %d sec.", 
                         main=expression(paste("Vorticity (", eta, ") and Streamlines (", psi, ")")),
                         z1.name="eta", z2.name="psi", 
                         color.palette=c("black", "purple", "blue", "white", "white", "forestgreen", "gold3", 
                                         "orange", "firebrick4")) {
  
  if(all(dim(z1) != dim(z2))) stop("mat1 and mat2 dims MUST be equal")
  if(length(dim(z1))!=2) stop("2D matrix is required")
  if(dim(z1)[1] != length(x)) stop("1th dim of matrices and length of x MUST be equal.")
  if(dim(z1)[2] != length(y)) stop("2th dim of matrices and length of x MUST be equal.")

  z2.lim      <- range(z2, finite=T)
  z2_levels <- pretty(z2.lim, z2.nlevels)
  if(length(z2_levels)>1) z2_levels <- z2_levels[-which(z2_levels==0, arr.ind = T)]
  eval(parse(text=paste0("legend.title <- expression(paste('(', ",z1.name,", ')'))")))
  
  filled.contour(x, y, z1, color=colorRampPalette(color.palette), zlim = z1.lim, 
                 xlab=xlab, ylab=ylab, key.title=title(line=1, main = legend.title),
                 plot.axes={axis(1); axis(2); contour(x, y, z2, levels=z2_levels, lwd=1, lty=2, add = T)})
  title(main=main,adj=0.5, line=3)
  title(sprintf(tlab, t), adj=0.5, line=1)
}

nplot.time.lapse <- function(x, y, t, mat1, mat2, interval=5, save=T, out_dir="animation", 
                    file="Rplot", size_factor=1.5, xlab= "x (meters)", ylab="z (meters)", 
                    tlab="Time= %d sec.", main=expression(paste("Vorticity (", eta, ") and Streamlines (", psi, ")")),
                    mat1.name="eta", mat2.name="psi",
                    color.palette=c("black", "purple", "blue", "white", "white", "forestgreen", "gold3", 
                                    "orange", "firebrick4")) {
  
  if(all(dim(mat1) != dim(mat2))) stop("mat1 and mat2 dims MUST be equal")
  if(length(dim(mat1))!=3) stop("3D matrix is required")
  if(dim(mat1)[1] != length(x)) stop("1th dim of matrices and length of x MUST be equal.")
  if(dim(mat1)[2] != length(y)) stop("2th dim of matrices and length of x MUST be equal.")
  if(dim(mat1)[3] != length(t)) stop("3th dim of matrices and length of t MUST be equal.")
  
  timestamps  <- t[t %% interval == 0]
  indices     <- which(t %in% timestamps, arr.ind = T)
  nlevels_max <- 40
  nlevels     <- vector(length = max(indices))
  nlevels[indices] <- ceiling((nlevels_max/length(indices))*(1:length(indices)))
  if(save) {
    dir.create(path = out_dir, showWarnings = F)
    ndigits <- nchar(as.character(max(seq_along(indices))))
    file <- sprintf("%s%%0%dd.png", file, ndigits)
    fnames <- vector(length = max(indices))
    fnames[indices] <- sprintf(file, seq_along(indices))
    pwidth  <- max(x)*size_factor
    pheigth <- max(y)*size_factor
  }
  
  for(i in indices) {
    time <- t[i];  m1 <- mat1[,,i]; m2  <- mat2[,,i]
    
    if(save) png(filename = file.path(out_dir, fnames[i]), width=pwidth, height=pheigth)

    plot.contour(x, y, time, m1, m2, nlevels[i], range(mat1),
                 xlab, ylab, tlab, main, mat1.name, mat2.name, color.palette);

    if(save) dev.off()
  }
}

nplot.psi.eta <- function(res, t=0, doubleside=F) {
  index.of.t <- which(res$time==t, arr.ind = T)
  z1 <- res$eta[,,index.of.t]
  z2 <- res$psi[,,index.of.t]
  x <- res$x
  z <- res$z
  if(doubleside) {
    z2 <- rbind(-z2[nrow(z2):2,], z2)
    z1 <- rbind(z1[nrow(z1):2,], z1)
    x  <- c(-x[length(x):2],x)
  }
  nplot.contour(x, z, t, z1, z2)
}

nplot.time.lapse.psi.eta <- function(res, interval=5, save=T, doubleside=F, out.dir="tl_psi_eta") {
  mat1 <- res$eta
  mat2 <- res$psi
  x <- res$x
  z <- res$z
  if(doubleside) {
    dim1 <- dim(mat2)[1]
    dim2 <- (dim(mat2)*c(2,1,1))-c(1,0,0)
    eta <- psi <- array(0,dim2)
    dim2 <- dim2[1]
    
    psi[1:(dim1-1),,] <- -mat2[dim1:2,,]
    psi[dim1:dim2,,] <- mat2[,,]
    
    eta[1:(dim1-1),,] <- mat1[dim1:2,,]
    eta[dim1:dim2,,] <- mat1[,,]
    
    mat1 <- t
    mat2 <- psi
    x   <- c(-x[length(x):2],x)
  }
  plot.time.lapse(x, z, res$time, mat1, mat2, interval, save, out.dir)
}

nplot.psi.t <- function(res, t=0, doubleside=F, t.lim=NULL, psi.nlevels=20) {
  index.of.t <- which(res$time==t, arr.ind = T)
  z1 <- res$t[,,index.of.t]
  z1[which(z1<0)] <- 0 # remove negative values for pot. temp. anomaly
  z2 <- res$psi[,,index.of.t]
  x <- res$x
  z <- res$z
  if(doubleside) {
    z2 <- rbind(-z2[nrow(z2):2,], z2)
    z1 <- rbind(z1[nrow(z1):2,], z1)
    x  <- c(-x[length(x):2],x)
  }
  if(is.null(t.lim)){
    t.lim <- range(z1, finite=T)
  }
  nplot.contour(x, z, t, z1, z2, z1.lim=t.lim, z2.nlevels=psi.nlevels,z1.name="theta - theta[0]",
               main=expression(paste("Pot.Temp.Anomaly (", theta,"-",theta[0] ,") and Streamlines (", psi, ")")),
               color.palette=c("white", "white","blue","purple", "forestgreen", "gold3", 
                                                "orange", "firebrick4"))
}

nplot.time.lapse.psi.t <- function(res, interval=5, save=T, doubleside=F, out.dir="tl_psi_t") {
  mat1 <- res$t
  mat1[which(mat1<0)] <- 0 # remove negative values
  mat2 <- res$psi
  x <- res$x
  z <- res$z
  if(doubleside) {
    dim1 <- dim(mat2)[1]
    dim2 <- (dim(mat2)*c(2,1,1))-c(1,0,0)
    t <- psi <- array(0,dim2)
    dim2 <- dim2[1]
    
    psi[1:(dim1-1),,] <- -mat2[dim1:2,,]
    psi[dim1:dim2,,] <- mat2[,,]
    
    t[1:(dim1-1),,] <- mat1[dim1:2,,]
    t[dim1:dim2,,] <- mat1[,,]
    
    mat1 <- t
    mat2 <- psi
    x   <- c(-x[length(x):2],x)
    
  }
  plot.time.lapse(x, z, res$time, mat1, mat2, interval, save, out.dir,
                  mat1.name="theta - theta[0]", mat2.name="psi",
                  main=expression(paste("Pot.Temp.Anomaly (", theta,"-",theta[0] ,") and Streamlines (", psi, ")")),
                  color.palette=c("white","blue","purple", "forestgreen", "gold3", 
                                  "orange", "firebrick4"))
}

ncreate.Movie <- function(run, run.name="run.A") {
  plot.time.lapse.psi.t(run, 1)
  plot.time.lapse.psi.eta(run, 1)
  system("./movie.sh")
}

nplot.energy.integrals <- function(x, ..., fact=10^5, labels=seq_along(list(...)),
                                  labels.x=c(7.6,8,8)[seq_along(list(...))], 
                                  labels.y=c(0.68, -0.35, -0.20)[seq_along(list(...))], 
                                  lty=seq_along(list(...)), lwd=2,
                                  cols=c(2,"forestgreen",seq_along(list(...))+3)) {
  args <- list(...)
  if(!all(sapply(args, length)==length(x))) stop("Length of arguments MUST be equal")
  len_args <- length(args)
  # http://stackoverflow.com/questions/16819275/obtain-names-of-variable-arguments-based-on-dot-dot-dot-in-function-r-deparse
  arg.names     <- setdiff(as.character(match.call()), as.character(match.call(expand.dots=F)))
  run.name      <- paste0("run", LETTERS)
  run.indices   <- apply(sapply(run.name, function(x) grepl(x, arg.names)),1, function(x) which(x==T, arr.ind = T))
  run.name      <- run.name[run.indices]
  cols2         <- cols[run.indices]
  
  arg.integrals <- arg.names
  for(r in paste0(run.name,"$")) arg.integrals <- gsub(r,"",arg.integrals, fixed=T)
  integral.names   <- c("intg_k", "intg_k_gz_phi")
  integral.expr    <- c("integral(integral(k,A))*da", "integral(integral((k-g*z*phi),A))*da")
  integral.indices <- sapply(arg.integrals, function(x) which(integral.names==x, arr.ind=T))
  run.integrals    <- paste0(run.name, " %->% ", integral.expr[integral.indices])
  mat <- matrix(c(...)/fact, ncol=len_args)
  matplot(x, mat, type="l", col=cols2, lty=lty, lwd=lwd, xlab = "t (minutes)", ylab=NA, las=1, xaxs="i")
  mtext(side = 2, bquote("Energy Integrals"~(10^.(floor(log10(fact)))~m^4~s^-2)), line = 2.5)
  # text(labels.x, labels.y, parse(text=run.integrals), col=cols2, cex=0.7)
  legend("topleft", parse(text=run.integrals), col=cols2, text.col = cols2, lty = 1:len_args,lwd=lwd, cex=0.7)
}

nplot.phi_square <- function(x, ..., fact=10^-2, labels=seq_along(list(...)),
                                  labels.x=c(7.6,8,8)[seq_along(list(...))], 
                                  labels.y=c(0.68, -0.35, -0.20)[seq_along(list(...))], 
                                  lty=seq_along(list(...)), lwd=2,
                                  cols=c(2,"forestgreen",seq_along(list(...))+3)) {
  args <- list(...)
  if(!all(sapply(args, length)==length(x))) stop("Length of arguments MUST be equal")
  len_args <- length(args)
  # http://stackoverflow.com/questions/16819275/obtain-names-of-variable-arguments-based-on-dot-dot-dot-in-function-r-deparse
  arg.names     <- setdiff(as.character(match.call()), as.character(match.call(expand.dots=F)))
  run.name      <- paste0("run", LETTERS)
  run.indices   <- apply(sapply(run.name, function(x) grepl(x, arg.names)),1, function(x) which(x==T, arr.ind = T))
  run.name      <- run.name[run.indices]
  cols2         <- cols[run.indices]
  
  arg.integrals <- arg.names
  for(r in paste0(run.name,"$")) arg.integrals <- gsub(r,"",arg.integrals, fixed=T)
  integral.names   <- c("intg_phi_square")
  integral.expr    <- c("integral(integral(phi^2,A))*da")
  integral.indices <- sapply(arg.integrals, function(x) which(integral.names==x, arr.ind=T))
  run.integrals    <- paste0(run.name, " %->% ", integral.expr[integral.indices])
  print(run.integrals)
  mat <- matrix(c(...)/fact, ncol=len_args)
  matplot(x, mat, type="l", col=cols2, lty=lty, lwd=lwd, xlab = "t (minutes)", ylab=NA, las=1, xaxs="i")
  mtext(side = 2, bquote("{"~phi^2~"}"~(10^.(floor(log10(fact)))~m^2)), line = 2.5)
  legend("topleft", parse(text=run.integrals), col=cols2, text.col = cols2, lty = 1:len_args,lwd=lwd, cex=0.7)
}

nplot.max.values <- function(x, ..., lty=seq_along(list(...)), lwd=2, legend.pos="topleft",
                            cols=c(2,"forestgreen",seq_along(list(...))+3), ylab="\"Max values\"") {
  args <- list(...)
  if(!all(sapply(args, length)==length(x))) stop("Length of arguments MUST be equal")
  len_args <- length(args)
  arg.names     <- setdiff(as.character(match.call()), as.character(match.call(expand.dots=F)))
  run.name      <- paste0("run", LETTERS)
  run.indices   <- unlist(apply(as.matrix(sapply(run.name, function(x) grepl(x, arg.names))),1, function(x) which(x==T, arr.ind = T)))
  run.name      <- run.name[run.indices]
  cols2         <- cols
  
  var.names   <- c("max_ws", "max_w", "max_u", "max_t", "max_psi", "max_eta",
                        "max_z_ws", "max_z_w", "max_z_u", "max_z_t", "max_z_psi", "max_z_eta")
  expr        <- paste0(c("V", "w", "u", "t", "psi","eta"),"[max]")
  expr        <- c(expr, paste0("Z[", expr ,"]" ))
  for(r in paste0(run.name,"$")) arg.names <- gsub(r,"",arg.names, fixed=T)
  var.indices <- sapply(arg.names, function(x) which(var.names==x, arr.ind=T))
  run.expr    <- paste0(run.name, " %->% ", expr[var.indices])
  
  mat <- matrix(c(...), ncol=len_args)
  matplot(x, mat, type="l", col=cols2, lty=lty, lwd=lwd, xlab = "t (minutes)", ylab=NA, las=1, xaxs="i", yaxs="i")
  mtext(side = 2, parse(text=ylab), line = 2.5)
  legend(legend.pos, parse(text=run.expr), col=cols2, text.col = cols2, lty = 1:len_args,lwd=lwd, cex=0.7)
}

