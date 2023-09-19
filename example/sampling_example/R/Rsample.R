#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# Rsample.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-12
# last mod 2013-04-28 15:49 DW
#

# Load rjags library and the wiener module
library(rjags)
load.module("wiener")

# Draw random samples with JAGS
samplemodel <- jags.model("../sample_model.txt",n.chains=1,n.adapt=0)
samples <- coda.samples(samplemodel,c("x"), n.iter=1000, thin=1)

# Now estimate parameters back
dat <- list(x=as.vector(samples[[1]][,1]), n=length(as.vector(samples[[1]][,1])))
inits <- list(alpha=1,tau=0.001,beta=0.5,delta=0)
estmodel <- jags.model("../wiener_model.txt",data=dat,inits=inits,n.chains=1,n.adapt=1)
estsamples <- coda.samples(estmodel,c("alpha","tau","beta","delta"),n.iter=1000,thin=1)

