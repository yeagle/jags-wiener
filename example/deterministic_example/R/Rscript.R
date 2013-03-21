#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# Rsample.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-12
# last mod 2013-03-20 11:04 DW
#

library(rjags)
load.module("wiener")

# data
dat <- list(t=seq(0.4,4.3,0.1))

# sample
jmodel <- jags.model("../deterministic_node.txt",data=dat,n.chains=1,n.adapt=0)
samples <- coda.samples(jmodel,c("x","y"), n.iter=10, thin=1)

# plot
plot(samples[[1]][1,1:40], type="b", pch=4, col="red")
x11()
plot(samples[[1]][1,41:80], type="b", pch=4, col="red")
