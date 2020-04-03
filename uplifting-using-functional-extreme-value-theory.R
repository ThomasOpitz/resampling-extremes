# Example implementation of lifting of extreme episodes using functional extreme-value theory with simulated data. 
# This is code accompanying the paper "Semi-parametric resampling with extremes" by T. Opitz, D. Allard, G. Mariéthoz.
# Author: Thomas Opitz, BioSP-INRAE (thomas.opitz@inrae.fr)

# Note: We here present an example usiing simulated t random fields on a grid since the input data used in the paper application (SAFRAN reanalyses) are not public.
# We use the median (taken over all pixels) for the risk functional r. 

library(extRemes)
library(fields)
library(RandomFields)

# 1) Generate synthetic data ####
ngrid = 100 # number of grid points in each direction in 2D
ns = ngrid * ngrid # number of pixels
m = 1000 # 1000 replicates
n.lift = 100 # number of lifted images to generate
xgrid = 1:ngrid/(ngrid + 1)
model = RMexp(var = 1, scale = .2)

nu = 3
x=RFsimulate(model = model, x = xgrid, y = xgrid, grid = TRUE, n = m, spConform = FALSE)
dim(x)
#generate scaling variables to obtain t-fields:
scales = sqrt(nu/rchisq(m, df = nu))
# rescale Gaussian fields:
for(j in 1:m){
  x[,,j] = scales[j] * x[,,j]
}
# show first field:
image.plot(xgrid, xgrid, x[,, 1], col = tim.colors(50), asp = 1)


# transform from array to matrix format:
loc = c() # grid locations
data_mat = matrix(NA, nrow = m, ncol = ns)
for(i in 1:ngrid){
  for(j in 1:ngrid){
    loc = rbind(loc, c(xgrid[i], xgrid[j]))
    data_mat[, (i-1)*ngrid+j] = x[i,j, ]  
  }
}
#plot(loc, pch = 19, cex = .05, xlab = "x", ylab = "y")

# 2) Set some global parameters ####

# p.marg: proportion of marginal data used to estimate the tail index
p.marg = 0.1
# u.marg.postprocess: threshold on Unif(-1,0)-scale below which to postprocess uniform data
# (can be equal to u.marg)
# (here we use the median)
u.marg.postprocess = -0.5
# p.aggr.extr: proportion of risk functional observations (here: median) to be considered as extreme 
p.aggr.extr=0.1
# We use this for estimating the r-extremal coefficient. 

# Lift spatial images to fixed return period of 1000 for the risk functional (median):
p.lift=1/1000 

# 3) Define some useful functions for the following analysis.

#density of generalized Pareto distribution:
dgp=function(x,sigma,xi){
  if(abs(xi)<10^{-4}){ # avoid numerical instabilities
    exp(-x/sigma)/sigma
  }else{
    result = ifelse(1+xi*x/sigma<=0,0,sigma^{-1}*(1+xi*x/sigma)^{-1/xi-1})
    ifelse(is.na(result) | (result == 0), .Machine$double.xmin, result)
  }
}

# cdf of generalized Pareto distribution ####
pgp = function(x,sigma,xi){
  if(abs(xi)<10^{-4}){ # avoid numerical instabilities
    1-exp(-x/sigma)
  }else{
    1- ifelse(1+xi*x/sigma<=0, 0, (1+xi*x/sigma)^{-1/xi})
  }
}

# quantile function of generalized Pareto distribution ####
qgp = function(x,sigma,xi){
  if(abs(xi)<10^{-4}){ # avoid numerical instabilities
    -sigma * log(1-x)
  }else{
    sigma/xi *((1-x)^{-xi}-1)
  }
}

# Calculate maximum likelihood estimate of tail index where the shape parameter sigma.fixed is given:
# We fix the maximum and minimum possible values of the tail index to sensible defaults (-1 and 5), and the initial value to 0.
get.tail.index.estimate = function(vec, p = 1-p.marg, sigma.fixed){
  thresh=quantile(vec,p)-10^{-10}
  is.exc=which(vec-thresh>0)
  fun2opt=function(par){
    -sum(log(dgp(x=vec[is.exc]-thresh,sigma=sigma.fixed,par)))
  }
  optim(par=0,fun2opt,method="Brent",lower=-1,upper=5)$par
}

# function to transform data to uniform distribution
# (kernel density estimate below the threshold, GP distribution above)
# vec is the data vector
# p is the exceedance probability used to fix the threshold
# n.dens number of points to discretize density below the threshold
get.uniform=function(vec,p,n.dens=2^12){
  thresh=quantile(vec,1-p)
  x.dens=density(vec,kernel="gaussian",n=n.dens) 
  idx.thresh=which.max(thresh-x.dens$x<0)-1
  dens.at.thresh=x.dens$y[idx.thresh] #density at the threshold
  sigma.gp=p/dens.at.thresh
  xi.gp=get.tail.index.estimate(vec,1-p,sigma.gp) 
  stepsize = diff(x.dens$x)[1] 
  Fx = cumsum(x.dens$y) * stepsize  # cdf values of kernel density estimate
  intervals = c(-Inf, x.dens$x + stepsize/2)
  intervals[length(intervals)] = Inf
  idx.vec = findInterval(vec[vec <= thresh], intervals)
  vec.unif = rep(NA, length(vec))
  vec.unif[vec <= thresh] = Fx[idx.vec]
  vec.unif[vec > thresh] = 1-p + p*pgp(vec[vec > thresh] - thresh, sigma = sigma.gp, xi = xi.gp)
  vec.unif
}


#function to backtransform lifted data from uniform scale to original marginal scale:
# (using kernel density estimate below the threshold, GP distribution above)
# (to avoid stocking intermediate results, we recalculate the original marginal distributions in this step)
# vec is the original data vector concatenated with the vector of lifted data on the uniform scale (to be backtransformed)
# idx.data gives the indices of the original data vector in vec
# idx.lift gives the indices of the lifted data on the uniform scale in vec
# p.marg is the exceedance probability used to fix the marginal threshold
# n.dens number of points to discretize density below the threshold
# data.orig is the original data vector
backtransform = function(vec, p.marg, idx.data, idx.lift, n.dens = 2^12){
  #first, establish the kernel density estimate
  thresh = quantile(vec[idx.data], 1 - p.marg) #threshold on original scal
  x.dens = density(vec[idx.data], kernel = "gaussian", n = n.dens)
  idx.thresh = which.max(thresh - x.dens$x < 0) - 1
  dens.at.thresh = x.dens$y[idx.thresh]
  sigma.gp = p.marg / dens.at.thresh
  xi.gp = get.tail.index.estimate(vec[idx.data], 1 - p.marg, sigma.gp)
  stepsize = diff(x.dens$x)[1] 
  Fx = cumsum(x.dens$y) * stepsize  # cdf values of kernel density estimate
  intervals = c(-.01, Fx)
  vec.lift = vec[idx.lift]
  idx.vec = findInterval(1 + vec.lift[vec.lift <= -p.marg], intervals)
  vec.orig = rep(NA, length(vec.lift))
  vec.orig[vec.lift <= -p.marg] = x.dens$x[idx.vec] 
  vec.orig[vec.lift > -p.marg] = thresh + qgp((vec.lift[vec.lift > -p.marg] + p.marg) / p.marg, sigma = sigma.gp, xi = xi.gp)
  vec.orig
}

# estimate the two GP parameters: 
# scale (defined through kernel density estimator)
# shape (maximum likelihood estimate)
# vec is the data vector
# p is the exceedance probability used to fix the marginal threshold
# n.dens number of points to discretize density below the threshold
get.gp.params = function(vec, p, n.dens = 2^12){
  thresh = quantile(vec, 1-p)
  x.dens = density(vec, kernel = "gaussian",n = n.dens)
  idx.thresh = which.max(thresh - x.dens$x < 0) - 1
  dens.at.thresh = x.dens$y[idx.thresh]
  sigma.gp = p / dens.at.thresh
  xi.gp = get.tail.index.estimate(vec, 1-p, sigma.gp)
  c(sigma.gp, xi.gp)
}

# 4) Transform data to uniform scale 
# (here we do it pixel by pixel, which is useful in case of spatial nonstationarities in practice)

# check if we get sensible values for tail indices:
# (this may take several minutes to calculate)
gp.params=apply(data_mat, 2, get.gp.params, p=p.marg)
# hist(gp.params[1,])  # scale parameter
# hist(gp.params[2,]) # tail index
# image.plot(matrix(gp.params[1,], nrow = ngrid), col = tim.colors(50))
# image.plot(matrix(gp.params[2,], nrow = ngrid), col = tim.colors(50))

# transform data to uniform scale at all pixels:
# (this may take several minutes)
unif_mat = apply(data_mat, 2, get.uniform, p = p.marg) - 1
# hist(unif_mat)

# get median value of standardized data for each replicate ####
medians = apply(unif_mat, 1, median)
# hist(medians)
# plot(medians,type="l")

# select extreme episodes:
u.r = quantile(medians, 1 - p.aggr.extr)
u.r
idx.epi.sel = which(medians > u.r)
# hist(medians[idx.epi.sel])

#generate n.lift lifted images: ####
#estimate the r-extremal coefficient: 
theta.r = p.aggr.extr / mean(medians > -p.aggr.extr)
theta.r
scale.new = -p.lift * theta.r # this is the new fixed median value
unif.lift = matrix(NA, nrow = n.lift, ncol = ns)
for(i in 1:length(idx.epi.sel)){
  idx.i = idx.epi.sel[i]
  eta.obs = unif_mat[idx.i,] / medians[idx.i]
  scale.obs = medians[idx.i]
  unif.lift[i,] = scale.new * eta.obs
  # find data to postprocess:
  data2rescale = unif_mat[idx.i, unif_mat[idx.i,] < u.marg.postprocess]
  scalingFactor = (1 + u.marg.postprocess * scale.new / scale.obs) / (1 + u.marg.postprocess)
  unif.lift[i, unif_mat[idx.i, ] < u.marg.postprocess] = -1 + scalingFactor * (1 + data2rescale)
}

# backtransform lifted images to original marginal scale:
data2backtransform=rbind(data_mat, unif.lift)
dim(data2backtransform)
idx.data = 1:m
idx.lift = (m+1):nrow(data2backtransform)
# (this may take several minutes)
data.lift = apply(data2backtransform, 2, backtransform, p.marg = p.marg, idx.data = idx.data, idx.lift = idx.lift)
# show first and last images, before and after lifting:
image.plot(matrix(data_mat[idx.epi.sel[1], ], nrow = ngrid), col = tim.colors(50))
image.plot(matrix(data.lift[1, ], nrow = ngrid), col = tim.colors(50))
image.plot(matrix(data_mat[idx.epi.sel[nrow(data.lift)], ], nrow = ngrid), col = tim.colors(50))
image.plot(matrix(data.lift[nrow(data.lift),], nrow = ngrid), col = tim.colors(50))

