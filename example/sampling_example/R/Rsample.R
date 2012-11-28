#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# Rsample.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-12
# last mod 2012-11-27 18:40 DW
#

library(rjags)
load.module("wiener")

# sample
samplemodel <- jags.model("../sample_model.txt",n.chains=1,n.adapt=0)
samples <- coda.samples(samplemodel,c("x"), n.iter=1000, thin=1)

# now estimate parameters back
dat <- list(x=as.vector(samples[[1]][,1]), n=length(as.vector(samples[[1]][,1])))
inits <- list(alpha=1,tau=0.001,beta=0.5,delta=0)
estmodel <- jags.model("../wiener_model.txt",data=dat,inits=inits,n.chains=1,n.adapt=1)
estsamples <- coda.samples(estmodel,c("alpha","tau","beta","delta"),n.iter=1000,thin=1)

