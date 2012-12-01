#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# Rsample.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-12
# last mod 2012-11-29 19:19 DW
#

library(rjags)
load.module("wiener")

# data
dat <- list(t=seq(0.4,4.3,0.1))

# sample
jmodel <- jags.model("../deterministic_node.txt",data=dat,n.chains=1,n.adapt=0)
samples <- coda.samples(jmodel,c("x","y"), n.iter=10, thin=1)

# now estimate parameters back
dat <- list(x=as.vector(samples[[1]][,1]), n=length(as.vector(samples[[1]][,1])))
inits <- list(alpha=1,tau=0.001,beta=0.5,delta=0)
estmodel <- jags.model("../wiener_model.txt",data=dat,inits=inits,n.chains=1,n.adapt=1)
estsamples <- coda.samples(estmodel,c("alpha","tau","beta","delta"),n.iter=1000,thin=1)

# plot
plot(samples[[1]][1,1:40], type="b", pch=4, col="red")
x11()
plot(samples[[1]][1,41:80], type="b", pch=4, col="red")
