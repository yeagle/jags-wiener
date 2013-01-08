#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# Rsample_var.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-01-07
# last mod 2013-01-07 17:00 DW
#
library(rjags)
load.module("wiener")

sample_model <- textConnection("model {
  x ~ dwieners(2,.2,.3,0,0.1)
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

# sample
samplemodel <- jags.model(sample_model,n.chains=1,n.adapt=0)
samples <- coda.samples(samplemodel,c("x"), n.iter=1000, thin=1)

# now estimate parameters back
dat <- list(x=as.vector(samples[[1]][,1]), n=length(as.vector(samples[[1]][,1])))
inits <- list(alpha=1,tau=0.001,beta=0.5,delta=0)
estmodel <- jags.model(wiener_model,data=dat,inits=inits,n.chains=1,n.adapt=1)
estsamples <- coda.samples(estmodel,c("alpha","tau","beta","delta"),n.iter=1000,thin=1)

