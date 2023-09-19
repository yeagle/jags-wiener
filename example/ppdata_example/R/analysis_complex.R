#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# analysis_complex.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-02
# last mod 2013-04-28 15:19 DW
#

# NOTE: the chains for eta.mu and eta.si do not converge
# Except for that, everything looks fine after burning and thinning
# NOTE2: the presented way of handling the data is rather complex. See the
# methods used in analysis_easymodel.R for a neater way, as it basically
# does the same analysis in almost half the lines of code.

library(rjags)
# load samples
load("rdata/complexmodel_250000.rdata")

a.mu <- samples[,1]
a.si <- samples[,2]
chi.mu <-samples[,3]
chi.si <-samples[,4]
dev <- samples[,5]
eta.mu <- samples[,6]
eta.si <- samples[,7]
theta.mu <- samples[,8]
theta.si <- samples[,9]
v.mu1 <- samples[,10]
v.mu2 <- samples[,11]
v.mu3 <- samples[,12]
v.mu4 <- samples[,13]
v.mu5 <- samples[,14]
v.si1 <- samples[,15]
v.si2 <- samples[,16]
v.si3 <- samples[,17]
v.si4 <- samples[,18]
v.si5 <- samples[,19]

burnin <- 20000
nsamp <- length(samples[[1]][,1])
thin <- 2000
# --> 345 iid samples
iid_samples <- data.frame()
for (chain in 1:3) {
  iid_samples <- rbind(iid_samples,
    data.frame(a.mu=as.array.mcmc.list(a.mu)[seq(burnin,nsamp,thin),chain],
               a.si=as.array.mcmc.list(a.si)[seq(burnin,nsamp,thin),chain],
               chi.mu=as.array.mcmc.list(chi.mu)[seq(burnin,nsamp,thin),chain],
               chi.si=as.array.mcmc.list(chi.si)[seq(burnin,nsamp,thin),chain],
               eta.mu=as.array.mcmc.list(eta.mu)[seq(burnin,nsamp,thin),chain],
               eta.si=as.array.mcmc.list(eta.si)[seq(burnin,nsamp,thin),chain],
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
              chi.mu,chi.si,eta.mu,eta.si,
              v.mu1,v.mu2,v.mu3,v.mu4,v.mu5,
              v.si1,v.si2,v.si3,v.si4,v.si5,dev)
namelist <- list("a.mu","a.si","theta.mu","theta.s",
              "chi.mu","chi.si","eta.mu","eta.si",
              "v.mu1","v.mu2","v.mu3","v.mu4","v.mu5",
              "v.si1","v.si2","v.si3","v.si4","v.si5","deviance")
for (i in 1:length(samplelist)) {
  x11()
  traceplot(samplelist[[i]], main=namelist[[i]], xlim=c(1,10000))
}

# Density
for (i in 1:length(samplelist)) {
  x11()
  densplot(samplelist[[i]], main=namelist[[i]])
}
# Density for iid_samples
for ( i in 1:length(iid_samples[1,])) {
  x11()
  plot(density(iid_samples[,i]), main=namelist[i])
}

# Autocorrelation: all in one
par(mfrow=c(4,5))
autocorr.plot(iid_samples[,1:19], auto.layout=F, ask=T, lag.max=50)
autocorr.plot(orig_samples[,1:19], auto.layout=F, ask=T, lag.max=50)


# Gelman Statistic
gelman.diag(samples)

# Deviance
traceplot(dev, main="Deviance for complexmodel")
traceplot(dev, xlim=c(48000,50000), ylim=c(-10000,0), main="Deviance for complexmodel")
plot(iid_samples$dev, type="l", lwd=2, col="red", main="Deviance")
iidN <- length(iid_samples$dev)
traceplot(mcmc.list(as.mcmc(iid_samples$dev[1:(iidN/3)]),
                    as.mcmc(iid_samples$dev[(iidN/3+1):(iidN/3*2)]),
                    as.mcmc(iid_samples$dev[(iidN/3*2+1):iidN])))

# Model estimates for v.mu
colMeans(orig_samples[,10:14])/10
sapply(orig_samples[,10:14],sd)/10
# For the burned and thinned dataset
colMeans(iid_samples[,9:13])/10
sapply(iid_samples[,9:13],sd)/10

length(iid_samples[,1])

# Plot the joint posterior for eta.mu and eta.si with color for the
# different chains
plot(as.numeric(eta.mu[[1]])[1:10000],as.numeric(eta.si[[1]])[1:10000], type="b", pch="+",
     lwd=1, lty=1, col="black")
points(as.numeric(eta.mu[[2]])[1:10000],as.numeric(eta.si[[2]])[1:10000], type="b", pch="x",
       lwd=1, lty=2, col="red")
points(as.numeric(eta.mu[[3]])[1:10000],as.numeric(eta.si[[3]])[1:10000], type="b", pch="o",
       lwd=1, lty=3, col="green")


# Example traceplot
layout(matrix(c(1,1,2),1,3,T))
traceplot(mcmc.list(as.mcmc(iid_samples$v.mu1[1:(iidN/3)]),
                    as.mcmc(iid_samples$v.mu1[(iidN/3+1):(iidN/3*2)]),
                    as.mcmc(iid_samples$v.mu1[(iidN/3*2+1):iidN])),
          main="Chain for the first driftrate mean",
          xlim=c(0,iidN/3))
autocorr.plot(iid_samples$v.mu1, auto.layout=F, lag.max=20)
