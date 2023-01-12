rm(list=ls())


# Github change


### DATA GENERATION ----

# Setup
n <- 500
Tmax <- 10

# Generating parameters
par <- data.frame(x0mn = -.8, y0mn = -.2, xAmn = .8, yAmn = 1,
         bx = -.35, by = -.25, gx = .1, gy = .2,
         x0v = .3, y0v = .6, xAv = .02, yAv = .08,
         x0xAcv = .046476, x0y0cv = .212132, x0yAcv = .030984,
         xAy0cv = .032863, y0yAcv = .131453, xAyAcv = .016,
         xmer = .1, ymer = .2, mercv = .02,
         xder = 0, yder = 0, dercv = 0)

# Run this for including dynamic error
par$xder <- .0436
par$yder <- .124
par$dercv <- sqrt(par$xder)*sqrt(par$yder)*0.2

# Initial conditions
# set.seed(10)
init <- with(par, 
             MASS::mvrnorm(n = n,
                           mu = c(x0mn, y0mn, xAmn, yAmn),
                           Sigma = matrix(c(x0v,    x0y0cv, x0xAcv, x0yAcv,
                                            x0y0cv, y0v,    xAy0cv, y0yAcv,
                                            x0xAcv, xAy0cv, xAv,    xAyAcv,
                                            x0yAcv, y0yAcv, xAyAcv, yAv), nrow=4)) )
colnames(init) <- c("x0", "y0", "xA", "yA")

# Measurement errors
I <- diag(1, Tmax)
m.err <- with(par, matrix(c(xmer, mercv, mercv, ymer), nrow=2))
psi <- kronecker(m.err, I)

m.errors <- MASS::mvrnorm(n = n,
                          mu = rep(0, Tmax*2),
                          Sigma = psi)
colnames(m.errors) <- c(paste0("X", 1:Tmax), paste0("Y", 1:Tmax))


# Dynamic errors
d.err <- with(par, matrix(c(xder, dercv, dercv, yder), nrow=2))
dyner <- kronecker(d.err, I)

d.errors <- MASS::mvrnorm(n = n,
                          mu = rep(0, Tmax*2),
                          Sigma = dyner)
colnames(d.errors) <- c(paste0("X", 1:Tmax), paste0("Y", 1:Tmax))

# Create matrices
xmat <- cbind(init[,1], matrix(NA, ncol=Tmax-1, nrow=n))
ymat <- cbind(init[,2], matrix(NA, ncol=Tmax-1, nrow=n))

# Generate latent scores
for(i in 1:(Tmax-1)){
  xmat[,(i+1)] <- xmat[,i] + init[,3] + par$bx*xmat[,i] + par$gx*ymat[,i] + d.errors[,i]
  ymat[,(i+1)] <- ymat[,i] + init[,4] + par$by*ymat[,i] + par$gy*xmat[,i] + d.errors[,(i+Tmax-1)]
}

latent <- cbind(xmat, ymat)
colnames(latent) <- c(paste0("X", 1:Tmax), paste0("Y", 1:Tmax))

# Generate observed scores
manifest <- latent + m.errors



### MODEL ESTIMATION ----
library(OpenMx)

# Variable names
x_manif <- paste0("X", 1:Tmax)
y_manif <- paste0("Y", 1:Tmax)
x_lat <- paste0("x", 1:Tmax)
y_lat <- paste0("y", 1:Tmax)
x_change <- paste0("x_ch_", 2:Tmax)
y_change <- paste0("y_ch_", 2:Tmax)
init <- c("x0", "xA", "y0", "yA")

BLCS <- mxModel("BLCS", 
                type="RAM",
                mxData(observed = manifest, type="raw"),
                manifestVars=c(x_manif, y_manif),
                latentVars=c(init, x_lat, x_change, y_lat, y_change),
                # From latent to manifest
                mxPath(from=c(x_lat, y_lat), to=c(x_manif, y_manif), arrows=1, 
                       free=FALSE, values=1),
                # From latent t-1 to latent t
                mxPath(from=x_lat[1:(Tmax-1)], to=x_lat[2:Tmax], arrows=1,
                       free=FALSE, values=1),
                mxPath(from=y_lat[1:(Tmax-1)], to=y_lat[2:Tmax], arrows=1,
                       free=FALSE, values=1),
                # From latent change scores to latent variables
                mxPath(from=c(x_change, y_change), to=c(x_lat[2:Tmax], y_lat[2:Tmax]), arrows=1, 
                       free=FALSE, values=1),
                # Self-feedback effects
                mxPath(from=x_lat[1:(Tmax-1)], to=x_change, arrows=1, 
                       free=TRUE, values=par$bx, labels="bx"),
                mxPath(from=y_lat[1:(Tmax-1)], to=y_change, arrows=1, 
                       free=TRUE, values=par$by, labels="by"),
                # Coupling effects
                mxPath(from=x_lat[1:(Tmax-1)], to=y_change, arrows=1, 
                       free=TRUE, values=par$gy, labels="gy"),
                mxPath(from=y_lat[1:(Tmax-1)], to=x_change, arrows=1, 
                       free=TRUE, values=par$gx, labels="gx"),
                # From initial conditions to first measurement occasion
                mxPath(from=c("x0", "y0"), to=c(x_lat[1], y_lat[1]), arrows=1, 
                       free=FALSE, values=1),
                # From additive component to latent change scores
                mxPath(from="xA", to=x_change[1:(Tmax-1)], arrows=1, 
                       free=FALSE, values=1),
                mxPath(from="yA", to=y_change[1:(Tmax-1)], arrows=1, 
                       free=FALSE, values=1),
                # Measurement error variances and covariance
                mxPath(from=x_manif, arrows=2, free=TRUE, values=par$xmer, labels="xmer"),
                mxPath(from=y_manif, arrows=2, free=TRUE, values=par$ymer, labels="ymer"),
                mxPath(from=x_manif, to=y_manif, arrows=2, free=TRUE, values=par$mercv, labels="mercv"),
                # Initial conditions and additive component variances and covariance
                mxPath(from=init, arrows=2, connect="unique.pairs", free=T, 
                       values=c(par$x0v, par$x0xAcv, par$x0y0cv, par$x0yAcv, par$xAv, 
                                par$xAy0cv, par$xAyAcv, par$y0v, par$y0yAcv, par$yAv),
                       labels = c("x0v", "x0xAcv", "x0y0cv", "x0yAcv", "xAv", 
                                  "xAy0cv", "xAyAcv", "y0v", "y0yAcv", "yAv")),
                # Mean structure
                mxPath(from="one", to=init, arrows=1, free=TRUE, 
                       values=c(par$x0mn, par$xAmn, par$y0mn, par$yAmn),
                       labels=c("x0mn", "xAmn", "y0mn", "yAmn")) )

# For estimating dynamic error
BLCS <- mxModel(BLCS,
                mxPath(from=x_change, arrows=2, free=T, values=0, labels="xder"),
                mxPath(from=y_change, arrows=2, free=T, values=0, labels="yder"),
                mxPath(from=x_change, to=y_change, arrows=2, free=T, values=0, labels="dercv"))


fit <- mxTryHard(BLCS)
summary(fit)
