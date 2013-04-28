#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# Rsample_var.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-01-07
# last mod 2013-04-28 15:50 DW
#

# Load rjags library and the wiener module
library(rjags)
load.module("wiener")

#1
## Sample and compare

sample_model1 <- textConnection("model {
  x ~ dwiener(1,.2,.3,2)
}")
sample_model2 <- textConnection("model {
  x ~ dwieners(1,.2,.3,2,1)
}")
sample_model3 <- textConnection("model {
  x ~ dwieners(.3,.2,.3,.6,0.3)
}")
samplemodel1 <- jags.model(sample_model1,n.chains=1,n.adapt=0)
samples1 <- coda.samples(samplemodel1,c("x"), n.iter=1000, thin=1)
samplemodel2 <- jags.model(sample_model2,n.chains=1,n.adapt=0)
samples2 <- coda.samples(samplemodel2,c("x"), n.iter=1000, thin=1)
samplemodel3 <- jags.model(sample_model3,n.chains=1,n.adapt=0)
samples3 <- coda.samples(samplemodel3,c("x"), n.iter=1000, thin=1)
x1 <- samples1[[1]]
x2 <- samples2[[1]]
x3 <- samples3[[1]]
t.test(x1,x2,alternative="two.sided",paired=F)
t.test(x2,x3,alternative="two.sided",paired=F)
plot(density(x1), xlim=c(-1.5,1.5), ylim=c(0,4))
points(density(x2), type="l", col="grey")
points(density(x3), type="l", col="red")

#2
## Sample and estimate

sample_model <- textConnection("model {
  x ~ dwieners(.2,.2,.3,0,0.1)
}")

wiener_model <- textConnection("model {
  for (i in 1:n) {
    x[i] ~ dwieners(alpha,tau,beta,delta,0.1)
  }
  alpha ~ dunif(-5,5)
  tau ~ dunif(0,1)
  beta ~ dunif(0,1)
  delta ~ dunif(-5,5)
}")

# Sample
samplemodel <- jags.model(sample_model,n.chains=1,n.adapt=0)
samples <- coda.samples(samplemodel,c("x"), n.iter=1000, thin=1)

# Now estimate parameters back
dat <- list(x=as.vector(samples[[1]][,1]), n=length(as.vector(samples[[1]][,1])))
inits <- list(alpha=1,tau=0.001,beta=0.5,delta=0)
estmodel <- jags.model(wiener_model,data=dat,inits=inits,n.chains=1,n.adapt=1)
estsamples <- coda.samples(estmodel,c("alpha","tau","beta","delta"),n.iter=1000,thin=1)
