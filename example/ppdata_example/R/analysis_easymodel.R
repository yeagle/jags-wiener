#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# analysis_easymodel.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-02
# last mod 2013-03-20 10:57 DW
#

# NOTE: Everything looks fine after burning and thinning

library(rjags)
# load samples
load("rdata/easymodel_50000.rdata")

# note: there are a lot more elegant ways to reshape the data
a.mu <- samples[,1]
a.si <- samples[,2]
dev <- samples[,3]
theta.mu <- samples[,4]
theta.si <- samples[,5]
v.mu1 <- samples[,6]
v.mu2 <- samples[,7]
v.mu3 <- samples[,8]
v.mu4 <- samples[,9]
v.mu5 <- samples[,10]
v.si1 <- samples[,11]
v.si2 <- samples[,12]
v.si3 <- samples[,13]
v.si4 <- samples[,14]
v.si5 <- samples[,15]

# note: using the command window() for burning and thinnig is better
burnin <- 1000
nsamp <- length(samples[[1]][,1])
thin <- 1000
# --> 147 iid samples
iid_samples <- data.frame()
for (chain in 1:3) {
  iid_samples <- rbind(iid_samples,
    data.frame(a.mu=as.array.mcmc.list(a.mu)[seq(burnin,nsamp,thin),chain],
               a.si=as.array.mcmc.list(a.si)[seq(burnin,nsamp,thin),chain],
               theta.mu=as.array.mcmc.list(theta.mu)[seq(burnin,nsamp,thin),chain],
               theta.si=as.array.mcmc.list(theta.si)[seq(burnin,nsamp,thin),chain],
               v.mu1=as.array.mcmc.list(v.mu1)[seq(burnin,nsamp,thin),chain],
               v.mu2=as.array.mcmc.list(v.mu2)[seq(burnin,nsamp,thin),chain],
               v.mu3=as.array.mcmc.list(v.mu3)[seq(burnin,nsamp,thin),chain],
               v.mu4=as.array.mcmc.list(v.mu4)[seq(burnin,nsamp,thin),chain],
               v.mu5=as.array.mcmc.list(v.mu5)[seq(burnin,nsamp,thin),chain],
               v.si1=as.array.mcmc.list(v.si1)[seq(burnin,nsamp,thin),chain],
               v.si2=as.array.mcmc.list(v.si2)[seq(burnin,nsamp,thin),chain],
               v.si3=as.array.mcmc.list(v.si3)[seq(burnin,nsamp,thin),chain],
               v.si4=as.array.mcmc.list(v.si4)[seq(burnin,nsamp,thin),chain],
               v.si5=as.array.mcmc.list(v.si5)[seq(burnin,nsamp,thin),chain],
               dev=as.array.mcmc.list(dev)[seq(burnin,nsamp,thin),chain],
               chain=chain))
}
orig_samples <- rbind(cbind(as.data.frame(samples[[1]]),data.frame(chain=rep(1,nsamp))),
                      cbind(as.data.frame(samples[[2]]),data.frame(chain=rep(2,nsamp))),
                      cbind(as.data.frame(samples[[3]]),data.frame(chain=rep(3,nsamp))))

# Chains
samplelist <- list(a.mu,a.si,theta.mu,theta.si,
              v.mu1,v.mu2,v.mu3,v.mu4,v.mu5,
              v.si1,v.si2,v.si3,v.si4,v.si5,dev)
namelist <- list("a.mu","a.si","theta.mu","theta.s",
              "v.mu1","v.mu2","v.mu3","v.mu4","v.mu5",
              "v.si1","v.si2","v.si3","v.si4","v.si5","deviance")
for (i in 1:length(samplelist)) {
  x11()
  traceplot(samplelist[[i]], main=namelist[[i]], xlim=c(1,2000))
}

# Density
for (i in 1:length(samplelist)) {
  #x11()
  densplot(samplelist[[i]], main=namelist[[i]])
}
# Density for iid_samples
for ( i in 1:length(iid_samples[1,])) {
  #x11()
  plot(density(iid_samples[,i]), main=namelist[i])
}

# Autocorrelation: all in one
par(mfrow=c(3,5))
autocorr.plot(iid_samples[,1:15], auto.layout=F, ask=T, lag.max=50)
autocorr.plot(orig_samples[,1:15], auto.layout=F, ask=T, lag.max=50)

# Gelman Statistic
gelman.diag(samples)
gelman.diag(mcmc.list(as.mcmc(orig_samples$"v.mu[5]"[orig_samples$chain==1]),as.mcmc(orig_samples$"v.mu[5]"[orig_samples$chain==2]),as.mcmc(orig_samples$"v.mu[5]"[orig_samples$chain==3])))

# Deviance
traceplot(dev, main="Deviance for easymodel")
traceplot(dev, xlim=c(48000,50000), ylim=c(-10000,-4000), main="Deviance for easymodel")
plot(iid_samples$dev, type="l", lwd=2, col="red", main="Deviance")
iidN <- length(iid_samples$dev)
traceplot(mcmc.list(as.mcmc(iid_samples$dev[1:(iidN/3)]),
                    as.mcmc(iid_samples$dev[(iidN/3+1):(iidN/3*2)]),
                    as.mcmc(iid_samples$dev[(iidN/3*2+1):iidN])))

# Model estimates for v.mu
colMeans(orig_samples[,6:10])/10
sapply(orig_samples[,6:10],sd)/10
# For the burned and thinned dataset
colMeans(iid_samples[,5:9])/10
sapply(iid_samples[,5:9],sd)/10

length(iid_samples[,1])

# Plot for good and bad chain
layout(matrix(c(1,1,2,3,3,4),2,3,T))
# Bad chain from orig samples: deviance
#traceplot(dev, xlim=c(0,80), ylim=c(-5000,-1000), main="Original Deviance (Bad chain)")
traceplot(dev, xlim=c(40000,40000+iidN/3), ylim=c(-9600,-5300), 
          main="Original Deviance samples (Bad chain)")
autocorr.plot(orig_samples$dev, auto.layout=F, lag.max=20)
# Good chain from iid samples: deviance
traceplot(mcmc.list(as.mcmc(iid_samples$dev[1:(iidN/3)]),
                    as.mcmc(iid_samples$dev[(iidN/3+1):(iidN/3*2)]),
                    as.mcmc(iid_samples$dev[(iidN/3*2+1):iidN])),
          main="Deviance after thinning and burning (Good chain)",
          xlim=c(0,iidN/3), ylim=c(-9600,-5300))
autocorr.plot(iid_samples$dev, auto.layout=F, lag.max=20)

# v.si5
par(mfrow=c(2,1))
# Bad chain from orig samples: deviance
#traceplot(dev, xlim=c(0,80), ylim=c(-5000,-1000), main="Original Deviance (Bad chain)")
traceplot(v.si5, xlim=c(40000,40000+iidN/3),# ylim=c(-9000,-6000), 
          main="Original v.si(5) samples (Bad chain)")
# Good chain from iid samples: deviance
traceplot(mcmc.list(as.mcmc(iid_samples$v.si5[1:(iidN/3)]),
                    as.mcmc(iid_samples$v.si5[(iidN/3+1):(iidN/3*2)]),
                    as.mcmc(iid_samples$v.si5[(iidN/3*2+1):iidN])),
          main="v.si(5) after thinning and burning (Good chain)",
          xlim=c(0,iidN/3))
