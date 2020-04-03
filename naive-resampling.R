# Example implementation of naive resampling (and of its limitations) using simulated data (single-replicate and multi-replicate setting).
# This is code accompanying the paper "Semi-parametric resampling with extremes" by T. Opitz, D. Allard, G. Mariéthoz.
# Author: Thomas Opitz, BioSP-INRAE (thomas.opitz@inrae.fr)

# Note: we here use the known marginal distributions for resampling data and do not estimate them. 
# We here transform to exponential margins (as explained in the paper) to highlight problematic behavior. 


# clear R workspace 
# rm(list = ls())

library(fields)
library(RandomFields)

#set your working directory
setwd("~/research/output/lifting/")

# 1) Example with one single big spatial realization ####
# We simulate a centered Gaussian random field with exponential covariance.
# Support is a grid on the unit sphere. 
ngrid = 200
nsim = 1
xgrid = 1:ngrid/(ngrid + 1)
model = RMexp(var = 1, scale = .125)
set.seed(2)
x = RFsimulate(model, x=xgrid, y=xgrid, grid=TRUE, n=nsim) 
# plot the simulated image:
image.plot(z = matrix(x$variable1, nrow = ngrid), col = tim.colors(100))

# data transformed to standard exponential margins:
x.exp = qexp(pnorm(x$variable1))
# Histogram of data transformed to standard exponential margins:
hist(x.exp, breaks=25, prob=T, xlab="x", ylab="Density")
tmp=1:1000/50
lines(tmp, dexp(tmp), col="blue", lwd=3)
# zoom into the tail region:
hist(x.exp, breaks=30, xlim=c(3.5,13.5), ylim=c(0,0.045), prob=T, xlab="x", ylab="Density", main = "")
tmp=1:1000/50
lines(tmp, dexp(tmp), col="blue", lwd=3)
rug(qexp(pnorm(x$variable1)))

# put Gaussian data into a matrix:
x_mat=matrix(x$variable1,nrow=ngrid,ncol=ngrid)
# check image:
#filled.contour(xgrid, xgrid, x_mat, color.palette = tim.colors, nlevels=50)

# put data on exponential scale into matrix:
x.exp_mat = qexp(pnorm(x_mat))
#filled.contour(xgrid, xgrid, x.exp_mat,color.palette = tim.colors, nlevels=50,axes=F,asp=1)

# define training data vector:
x.train=as.numeric(x.exp_mat)
# extract order of training vectors:
ord.train=order(x.train)
# extract ranks of training vector:
rank.train=rank(x.train, ties.method="random")

# sample vector of same length with same marginal distributions:
x.resample=rexp(n=length(x.train))
hist(x.resample)
summary(x.train)
summary(x.resample)
hist(x.train,breaks=50)
hist(x.resample,breaks=50)

# use the ranks from the training data pattern to generate spatial dependence:
x.out=sort(x.resample)[rank.train]
x.out_mat=matrix(x.out,ncol=ngrid)

# show full images (original and resample): ####
zlim=range(x.exp_mat, x.out_mat)
par(mfrow=c(1,2),mar=rep(.5,4),oma=rep(.5,4))
image.plot(x.exp_mat, zlim=zlim, nlevel=20,asp=1,axes=F,cex.axis=2)
image.plot(x.out_mat,nlevel=20,asp=1,axes=F,cex.axis=2,cex.lab=2)
# dev.off()

# zoom in on what is going on around the global maximum:
tmp = which(x.out_mat == max(x.out_mat), arr.ind = TRUE)
tmp #index of global maximum
m = 8 # number of pixels around global maximum
i1max = (tmp[1,1]-m):(tmp[1,1] + m)
i2max = (tmp[1,2]-m):(tmp[1,2] + m)
zlim = range(x.exp_mat[i1max, i2max], x.out_mat[i1max, i2max])
par(mfrow = c(1, 2), mar = rep(.5, 4), oma = rep(.5, 4))
image.plot(x.exp_mat[i1max, i2max], zlim = zlim, nlevel = 20, asp=1, axes=F, cex.axis = 2)
image.plot(x.out_mat[i1max, i2max], zlim = zlim, nlevel = 20, asp = 1, axes = F, cex.axis = 2, cex.lab = 2)
# dev.off()

# Same example as before, but now with 100 many replicates of the Gaussian process ####
nsim = 100
ngrid = 100
xgrid = 1:ngrid/(ngrid+1)
model = RMexp(var = 1, scale = .125)
# use modus_operandi = "precise" to avoid approximation errors
x = RFsimulate(model, x = xgrid, y = xgrid, grid = TRUE, n = nsim, exactness = TRUE, modus_operandi = "precise") 
x = as.array(x)

x.exp = qexp(pnorm(as.numeric(x)))
par(mfrow = c(1,1))
hist(x.exp, breaks=25, prob=T, xlab="x", ylab="Density")
tmp = 1:1000/50
lines(tmp, dexp(tmp), col="blue", lwd=3)

# zoom into the tail:
hist(x.exp, breaks = 100, xlim = c(3.5,13.5), ylim = c(0,0.045), prob = T, xlab = "x", ylab = "Density", main = "")
tmp = 1:1000/50
lines(tmp, dexp(tmp), col = "blue", lwd = 3)
rug(x.exp)

x.exp_array = array(x.exp, dim = dim(x))

x.train = as.numeric(x.exp_array)
hist(x.train)
ord.train=order(x.train)
rank.train=rank(x.train, ties.method="random")

x.resample = rexp(n=length(x.train))
hist(x.resample)
summary(x.train)
summary(x.resample)
#hist(x.train,breaks=50)
#hist(x.resample,breaks=50)

x.out = sort(x.resample)[rank.train]
x.out_array = array(x.out, dim = dim(x))

# show behavior around the global maximum
tmp = which(max(x.out_array) == x.out_array, arr.ind=TRUE)
tmp
m = 8
i1max = (tmp[1, 1] - m):(tmp[1, 1] + m)
i2max = (tmp[1, 2] - m):(tmp[1, 2] + m)
zlim = range(x.exp_array[i1max,i2max, tmp[1, 3]], x.out_array[i1max, i2max,  tmp[1,3]])
range(x.exp_array[i1max, i2max, tmp[1, 3]])
range(x.out_array[i1max,i2max, tmp[1, 3]])
par(mfrow = c(1, 2), mar = rep(.5, 4), oma = rep(.5, 4))
image.plot(x.exp_array[i1max, i2max, tmp[1, 3]], zlim = zlim, nlevel = 20, asp = 1, axes = F, cex.axis = 2, cex.lab = 2)
image.plot(x.out_array[i1max, i2max, tmp[1, 3]], nlevel = 20, asp = 1, axes = F, cex.axis = 2,cex.lab = 2)
#dev.off()

# Simulate 1000 samples (each having 100 replicates of the spatial image) of maxima and spacings around the maximum ####
## First, with Gaussian random fields:
## careful, this takes quite long...

ngrid = 30
nsim = 100
xgrid = 1:ngrid/(ngrid+1)
model = RMexp(var = 1, scale = .2)

max.data = c()
max.sim = c()
space.data = c()
space.sim = c()
for(i in 1:1000){
  print(i)
  set.seed(i)
  # use RPdirect (Cholesky decomposition) to ensure good simulation quality
  x=RFsimulate(model = RPdirect(model), x=xgrid, y=xgrid, grid=TRUE, n=nsim, spConform=FALSE, max_variab=1000000, exactness = TRUE, modus_operandi = "neurotic")
  x.train=sort(qexp(pnorm(as.numeric(x))), decreasing = TRUE)
  x.resample=sort(rexp(length(x.train)), decreasing = TRUE)
  max.data[i] = x.train[1] # global maximum
  max.sim[i] = x.resample[1] # global maximum
  space.data[i] = diff(x.train[1:2]) # spacing after maximum
  space.sim[i] = diff(x.resample[1:2]) # spacing after maximum
}

# Maxima:
par(mfrow = c(1, 1), lwd = 2, mar = c(4.5,4.5,1,1), cex.lab = 1.25, cex.axis = 1.25)
xlim = range(max.data, max.sim) 
breaks = xlim[1] - 0.0001 + 0:ceiling(diff(xlim)/0.75)*0.75
xlim = xlim + c(-.75, .75)
tmp1 = hist(max.data, breaks = breaks, prob = TRUE)
tmp2 = hist(max.sim, breaks = breaks, prob = TRUE)
ylim = range(tmp1$density, tmp2$density)
hist(max.data, xlim = xlim, ylim = ylim, breaks = breaks, border = "black", prob = TRUE, xlab = "Sample maximum", main = "")
hist(max.sim, breaks = breaks, border = "red", add = TRUE, prob = TRUE)

# Spacings:
par(lwd = 2, mar = c(4.5,4.5,1,1), cex.lab = 1.25, cex.axis = 1.25)
xlim = range(space.data, space.sim)
breaks = xlim[2] + 0.0001 - 0:ceiling(diff(xlim)/0.5)*0.5
xlim = xlim + c(-.5, .5)
tmp1 = hist(space.data, breaks = breaks, prob = TRUE)
tmp2 = hist(space.sim, breaks = breaks, prob = TRUE)
ylim = range(tmp1$density, tmp2$density)
hist(space.data, xlim = xlim, ylim = ylim, breaks = breaks, border = "black", prob = TRUE, xlab = "Spacing", main = "")
hist(space.sim, breaks = breaks, border = "red", prob = TRUE, add = TRUE)

## Second, with t-random fields (df = 3): ####
ngrid = 30
nsim = 100
xgrid = 1:ngrid/(ngrid+1)
model = RMexp(var = 1, scale = .2)

nu = 3

max.data = c()
max.sim = c()
space.data = c()
space.sim = c()
for(i in 1:1000){
  print(i)
  set.seed(i)
  x=RFsimulate(model = RPdirect(model), x=xgrid, y=xgrid, grid=TRUE, n=nsim, 
               spConform=FALSE, max_variab=1000000, exactness = TRUE, modus_operandi = "neurotic")
  scales = sqrt(nu/rchisq(nsim, df = nu))
  for(j in 1:nsim){
    x[,,j] = scales[j] * x[,,j]
  }
  x_train=sort(qexp(pt(as.numeric(x), df = nu)), decreasing = TRUE)
  x_resample=sort(rexp(length(x_train)), decreasing = TRUE)
  max.data[i] = x_train[1]
  max.sim[i] = x_resample[1]
  space.data[i] = diff(x_train[1:2])
  space.sim[i] = diff(x_resample[1:2])
}

# plot a histogram of the last dataset (i=1000) transformed to standard exponential margins:
par(lwd = 2, mar = c(4.5,4.5,1,1), cex.lab = 1.25, cex.axis = 1.25)
hist(x.train, breaks=25, prob=T, xlab="x (Sample)", ylab="Density", main = "") 
tmp=1:1000/50
lines(tmp, dexp(tmp), col="blue", lwd=3)
#looks fine
# resample: 
par(lwd = 2, mar = c(4.5,4.5,1,1), cex.lab = 1.25, cex.axis = 1.25)
hist(x.resample, breaks=25, prob=T, xlab="x (Resample)", ylab="Density", main = "")
tmp=1:1000/50
lines(tmp, dexp(tmp), col="blue", lwd=3)
# looks fine, too

# Maxima:
par(lwd = 2, mar = c(4.5,4.5,1,1), cex.lab = 1.25, cex.axis = 1.25)
xlim = range(max.data, max.sim) 
breaks = xlim[1] - 0.0001 + 0:ceiling(diff(xlim)/0.75)*0.75
xlim = xlim + c(-.75, .75)
tmp1 = hist(max.data, breaks = breaks, prob = TRUE)
tmp2 = hist(max.sim, breaks = breaks, prob = TRUE)
ylim = range(tmp1$density, tmp2$density)
hist(max.data, xlim = xlim, ylim = ylim, breaks = breaks, border = "black", prob = TRUE, xlab = "Sample maximum", main = "")
hist(max.sim, breaks = breaks, border = "red", add = TRUE, prob = TRUE)

# Spacings:
par(lwd = 2, mar = c(4.5,4.5,1,1), cex.lab = 1.25, cex.axis = 1.25)
xlim = range(space.data, space.sim)
breaks = xlim[2] + 0.0001 - 0:ceiling(diff(xlim)/0.5)*0.5
xlim = xlim + c(-.5, .5)
tmp1 = hist(space.data, breaks = breaks, prob = TRUE)
tmp2 = hist(space.sim, breaks = breaks, prob = TRUE)
ylim = range(tmp1$density, tmp2$density)
hist(space.data, xlim = xlim, ylim = ylim, breaks = breaks, border = "black", prob = TRUE, xlab = "Spacing", main = "")
hist(space.sim, breaks = breaks, border = "red", prob = TRUE, add = TRUE)

