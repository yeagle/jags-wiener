#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# analysis_easymodel.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-02
# last mod 2013-04-28 22:55 DW
#

# load all the required libraries, starting with rjags
library(rjags)

# load samples, which have been previously generated with the easymodel.R
# script
load("rdata/easymodel_50000.rdata")

# inspect raw chains
namelist <- list("a.mu","a.si","theta.mu","theta.s",
              "v.mu1","v.mu2","v.mu3","v.mu4","v.mu5",
              "v.si1","v.si2","v.si3","v.si4","v.si5","deviance")
# trace
for (i in 1:15) {
  x11()
  traceplot(samples[,i], main=namelist[[i]], xlim=c(1,2000))
}
# density
for (i in 1:15) {
  #x11()
  densplot(samples[,i], main=namelist[[i]])
}

# burn and thin the chains (preferably after inspection of the unchanged
# chains, to determine a reasonable burning and thinning interval)
burnin <- 1000
nsamp <- length(samples[[1]][,1])
thin <- 1000
iid_samples <- window(samples, burnin,nsamp,thin)

# density for iid_samples
for (i in 1:15) {
  #x11()
  densplot(iid_samples[,i], main=namelist[i])
}

# autocorrelation for chain 1
par(mfrow=c(3,5))
autocorr.plot(samples[[1]], auto.layout=F, ask=T, lag.max=50)
autocorr.plot(iid_samples[[1]], auto.layout=F, ask=T, lag.max=50)

# gelman statistic
gelman.diag(samples)
gelman.diag(iid_samples)
gelman.diag(samples[,9])

# model estimates for v.mu
round(summary(iid_samples)$statistic,4)

# deviance plot for good and bad chain
layout(matrix(c(1,1,2,3,3,4),2,3,T))
# bad chain from orig samples: deviance
traceplot(samples[,15], xlim=c(48000,48040),
          main="Original Deviance samples")
autocorr.plot(samples[,15][[1]], auto.layout=F, lag.max=20)
# good chain from iid samples: deviance
traceplot(iid_samples[,15],
          main="Deviance after thinning and burning")
autocorr.plot(iid_samples[,15][[1]], auto.layout=F, lag.max=20)

# v.si5
par(mfrow=c(2,1))
# bad chain from orig samples: deviance
traceplot(samples[,14], xlim=c(1000,1040),
          main="Original v.si(5) samples (Bad chain)")
# Good chain from iid samples: deviance
traceplot(iid_samples[,14],
          main="v.si(5) after thinning and burning (Good chain)")
