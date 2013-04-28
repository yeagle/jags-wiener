# This script demonstrates how to call JAGS from R, load the jags-wiener
# module, and fit a diffusion model with JAGS.

# Load rjags library and the jags wiener module
library(rjags)
load.module("wiener")

# Read original data from the experiment
dat <- read.table("data/exp1_01_m23r_2012-10-31-110958.txt", skip=8, header=T)
dat <- dat[21:218,] # Throw out trainig trials
dat$RT[dat$Key == "a"] = -dat$RT[dat$Key == "a"]
dat <- data.frame(Id=dat$Id,Cond=as.numeric(dat$Stim), RT=dat$RT)

# Put data in a list for simple use in the jags.model() command
dat <- list(N=198, y=dat$RT, C=3, cond=dat$Cond)

# We will run three MCMC chains, so make tree sets of initial values.
# The initial value sets could be randomly generated or deterministically
# set. Ideally, initial values are overdispersed compared to the posterior
# distribution (i.e., they should bracket most of its mass).
inits1 <- list(alpha=1, tau=0.001, delta=c(0.3,.5,-2))
inits2 <- list(alpha=2, tau=0.003, delta=c(1,-1,0))
inits3 <- list(alpha=0.2, tau=0.004, delta=c(-2,0,1))
inits <- list(inits1,inits2,inits3)

# Set up the analysis by passing the model file, data list, and inits
# to jags.model().
easy_model <- jags.model("../easy_model.txt", dat, inits, n.chains=3)

# Run the MCMC engine with jags.samples(), and request chains for all
# nodes of interest, including the special diagnostic node "deviance".
# Alternatively, one can use the command coda.samples(), which works the
# same but returns the samples in a different data type
samples <- jags.samples(easy_model, n.iter=5000, 
                        variable.names=c("alpha", "tau", "delta"))

# Finally, inspect the chains and use posterior distributions for
# inference.

# Reshape chains as mcmc.lists
alpha <- samples$alpha[,,]
alpha <- mcmc.list(mcmc(alpha[,1]), mcmc(alpha[,2]), mcmc(alpha[,3]))
tau <- samples$tau[,,]
tau <- mcmc.list(mcmc(tau[,1]), mcmc(tau[,2]), mcmc(tau[,3]))
delta1 <- samples$delta[1,,]
delta1 <- mcmc.list(mcmc(delta1[,1]), mcmc(delta1[,2]), mcmc(delta1[,3]))
delta2 <- samples$delta[2,,]
delta2 <- mcmc.list(mcmc(delta2[,1]), mcmc(delta2[,2]), mcmc(delta2[,3]))
delta3 <- samples$delta[3,,]
delta3 <- mcmc.list(mcmc(delta3[,1]), mcmc(delta3[,2]), mcmc(delta3[,3]))

# Look at chain trace and posterior distribution
x11() ; par(mfrow=c(1,2))
traceplot(alpha, main="Chains", xlab="n", ylab="Value")
densplot(alpha, main="Posterior Distribution", xlab="Value", ylab="Density")
x11() ; par(mfrow=c(1,2))
traceplot(tau, main="Chains", xlab="n", ylab="Value")
densplot(tau, main="Posterior Distribution", xlab="Value", ylab="Density")
# Dark condition
x11() ; par(mfrow=c(1,2))
traceplot(delta1, main="Chains", xlab="n", ylab="Value")
densplot(delta1, main="Posterior Distribution", xlab="Value", ylab="Density")
# Bright condition
x11() ; par(mfrow=c(1,2))
traceplot(delta2, main="Chains", xlab="n", ylab="Value")
densplot(delta2, main="Posterior Distribution", xlab="Value", ylab="Density")
# Mixed condition
x11() ; par(mfrow=c(1,2))
traceplot(delta3, main="Chains", xlab="n", ylab="Value")
densplot(delta3, main="Posterior Distribution", xlab="Value", ylab="Density")
