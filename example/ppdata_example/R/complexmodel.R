#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# complexmodel.R
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-11-01
# last mod 2013-05-09 14:22 DW
#

# This script demonstrates how to call JAGS from R, load the jags-wiener
# module, and fit a diffusion model with JAGS. The data set and model are
# described in the companion paper.

# First, choose the number of samples to draw
NSAMP <- 250000

# Then, load all the required libraries, starting with rjags
library(rjags)         # load the JAGS library
load.module("wiener")  # also load the jags-wiener library
load.module("dic")     # also load the DIC library

# Now, read in the data file (data from Vandekerckhove et al., 2007)
dat <- read.jagsdata("../pp_data.txt")

# We will run three MCMC chains, so make tree sets of initial values.
# The initial value sets could be randomly generated or deterministically
# set. Ideally, initial values are overdispersed compared to the posterior
# distribution (i.e., they should bracket most of its mass).
# It is not necessary to give initial values for all possible parameters,
# just those at the highest levels of the model. JAGS will automatically
# generate initial values for missing parameters. Note that JAGS does not
# generate initial values randomly, but deterministically based on the
# higher-level parameters, so in order to get different chains, you need
# to provide initial values for all nodes that have fixed priors.
# See model file ../complexmodel.txt for names of nodes.
cminits1 <- list(v.mu=c(0,0,0,0,0),v.si=c(1.5,1.5,1.5,1.5,1.5),
               eta.mu=0.1,eta.si=0.1,chi.mu=0.1,chi.si=0.1,
               a.mu=1.5,a.si=.5,theta.mu=.01,theta.si=.1)
cminits2 <- list(v.mu=c(1,2,-3,5,3),v.si=c(1.0,1.1,1.5,0.8,2.3),
               eta.mu=0.2,eta.si=0.15,chi.mu=0.1,chi.si=0.2,
               a.mu=1.0,a.si=.3,theta.mu=.02,theta.si=.08)
cminits3 <- list(v.mu=c(-1.-2.6,-1,3,-4,-5),v.si=c(2.1,2,0.5,1.5,2),
               eta.mu=0.2,eta.si=0.2,chi.mu=0.01,chi.si=0.01,
               a.mu=2.5,a.si=.9,theta.mu=.06,theta.si=.2)
cminits <- list(cminits1,cminits2,cminits3)  # concatenate inits lists

# Set up the analysis by passing the model file, data list, and inits
# to jags.model().
complexmodel <- jags.model("../complexmodel.txt", dat, cminits, 3, 0)

# Run the MCMC engine with coda.samples(), and request chains for all
# nodes of interest, including the special diagnostic node "deviance".
# Alternatively, one can use the command jags.samples(), which works the
# same but returns the samples in a different data type
#
# If you get a WARNING or NOTE that says 'Adaptation incomplete' or
# 'Stopping adaptation', ignore it. This message is a JAGS warning that the
# burn-in phase is not completed and the user will need to manually discard
# some of the initial samples (in practice, we rarely rely on JAGS'
# judgment regarding the burn-in phase; it is currently based solely on the
# number of samples drawn, with no regard for actual convergence quality)
samples <- coda.samples(complexmodel, c("v.mu","v.si","a.mu","a.si",
                                     "theta.mu","theta.si", "eta.mu",
                                     "eta.si", "chi.mu", "chi.si",
                                     "deviance"), NSAMP, 1)

# Write the chains to disk.
save(list=c("dat","cminits","complexmodel"),file="rdata/complexmodel.rdata")
save(samples,file=paste("rdata/complexmodel_", as.character(NSAMP), ".rdata", sep=""))

# Finally, inspect the chains and use posterior distributions for
# inference.
# Look at analysis_complex.R for an example analysis of the chains.
