#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# complexmodel.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-01
# last mod 2012-12-06 14:20 DW
#

####### Number of Samples to draw #######
####### ######################### #######
NSAMP <- 250000
####### ######################### #######

library(rjags)
load.module("wiener")
load.module("dic")

# data and inits
dat <- read.jagsdata("../pp_data.txt")
cminits1 <- list(v.mu=c(0,0,0,0,0),v.si=c(1.5,1.5,1.5,1.5,1.5),
               eta.mu=0.1,eta.si=0.1,chi.mu=0.1,chi.si=0.1,
               a.mu=1.5,a.si=.5,theta.mu=.01,theta.si=.1)
cminits2 <- list(v.mu=c(1,2,-3,5,3),v.si=c(1.0,1.1,1.5,0.8,2.3),
               eta.mu=0.2,eta.si=0.15,chi.mu=0.1,chi.si=0.2,
               a.mu=1.0,a.si=.3,theta.mu=.02,theta.si=.08)
cminits3 <- list(v.mu=c(-1.-2.6,-1,3,-4,-5),v.si=c(2.1,2,0.5,1.5,2),
               eta.mu=0.2,eta.si=0.2,chi.mu=0.01,chi.si=0.01,
               a.mu=2.5,a.si=.9,theta.mu=.06,theta.si=.2)
cminits <- list(cminits1,cminits2,cminits3)

complexmodel <- jags.model("../complexmodel.txt", dat, cminits, 3, 0)

samples <- coda.samples(complexmodel, c("v.mu","v.si","a.mu","a.si",
                                     "theta.mu","theta.si", "eta.mu",
                                     "eta.si", "chi.mu", "chi.si",
                                     "deviance"), NSAMP, 1)

save(list=c("dat","cminits","complexmodel"),file="rdata/complexmodel.rdata")
save(samples,file=paste("rdata/complexmodel_", as.character(NSAMP), ".rdata", sep=""))
