#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# easymodel.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-01
# last mod 2012-12-06 14:20 DW
#

####### Number of Samples to draw #######
####### ######################### #######
NSAMP <- 50000 
####### ######################### #######

library(rjags)
load.module("wiener")
load.module("dic")

# data and inits
dat <- read.jagsdata("../pp_data.txt")
eminits1 <- list(v.mu=c(0,0,0,0,0),v.si=c(1.5,1.5,1.5,1.5,1.5),
               a.mu=1.5,a.si=.5,theta.mu=.01,theta.si=.1)
eminits2 <- list(v.mu=c(1,2,-3,5,-3),v.si=c(1.0,1.1,1.5,0.8,2.3),
               a.mu=1.9,a.si=.3,theta.mu=.02,theta.si=.01)
eminits3 <- list(v.mu=c(-4.5,-1,3,-4,-5),v.si=c(1.1,1,0.5,1.5,2),
               a.mu=2.5,a.si=.9,theta.mu=.06,theta.si=.04)
eminits <- list(eminits1,eminits2,eminits3)

easymodel <- jags.model("../easymodel.txt", dat, eminits, 3, 0)


samples <- coda.samples(easymodel, c("v.mu","v.si","a.mu","a.si",
                                     "theta.mu","theta.si", "deviance"), NSAMP, 1)

save(list=c("dat","eminits","easymodel"),file="rdata/easymodel.rdata")
save(samples,file=paste("rdata/easymodel_", as.character(NSAMP), ".rdata", sep=""))
