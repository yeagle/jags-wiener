# load rjags lib and the jags wiener module
library(rjags)
load.module("wiener")

# read original data
dat <- read.table("data/exp1_01_m23r_2012-10-31-110958.txt", skip=8, header=T)
dat <- dat[21:218,] # Throw out trainig trials
dat$RT[dat$Key == "a"] = -dat$RT[dat$Key == "a"]
dat <- data.frame(Id=dat$Id,Cond=as.numeric(dat$Stim), RT=dat$RT)

# data and inits
dat <- list(N=198, y=dat$RT, C=3, cond=dat$Cond)
inits1 <- list(alpha=1, tau=0.001, delta.mu=c(1,-1,-2), delta.si=0.5)
inits2 <- list(alpha=2, tau=0.003, delta.mu=c(-1,1,0), delta.si=0.8)
inits3 <- list(alpha=0.2, tau=0.004, delta.mu=c(2,1.5,-2), delta.si=0.2)
inits <- list(inits1,inits2,inits3)

# create jags model and sample from it
easy_model <- jags.model("../easy_model.txt", dat, inits, n.chains=3)
samples <- jags.samples(easy_model, n.iter=5000, 
                        variable.names=c("alpha", "tau", 
                                         "delta.mu", "delta.si"))

# store chains as mcmc.lists
alpha <- samples$alpha[,,]
alpha <- mcmc.list(mcmc(alpha[,1]), mcmc(alpha[,2]), mcmc(alpha[,3]))
tau <- samples$tau[,,]
tau <- mcmc.list(mcmc(tau[,1]), mcmc(tau[,2]), mcmc(tau[,3]))
delta.mu1 <- samples$delta.mu[1,,]
delta.mu1 <- mcmc.list(mcmc(delta.mu1[,1]), mcmc(delta.mu1[,2]), mcmc(delta.mu1[,3]))
delta.mu2 <- samples$delta.mu[2,,]
delta.mu2 <- mcmc.list(mcmc(delta.mu2[,1]), mcmc(delta.mu2[,2]), mcmc(delta.mu2[,3]))
delta.mu3 <- samples$delta.mu[3,,]
delta.mu3 <- mcmc.list(mcmc(delta.mu3[,1]), mcmc(delta.mu3[,2]), mcmc(delta.mu3[,3]))
delta.si <- samples$delta.si[,,]
delta.si <- mcmc.list(mcmc(delta.si[,1]), mcmc(delta.si[,2]), mcmc(delta.si[,3]))

# Look at chains and posterior
x11() ; par(mfrow=c(1,2))
traceplot(alpha, main="Chains", xlab="n", ylab="Value")
densplot(alpha, main="Posterior Distribution", xlab="Value", ylab="Density")
x11() ; par(mfrow=c(1,2))
traceplot(tau, main="Chains", xlab="n", ylab="Value")
densplot(tau, main="Posterior Distribution", xlab="Value", ylab="Density")
# Dark condition
x11() ; par(mfrow=c(1,2))
traceplot(delta.mu1, main="Chains", xlab="n", ylab="Value")
densplot(delta.mu1, main="Posterior Distribution", xlab="Value", ylab="Density")
# Bright condition
x11() ; par(mfrow=c(1,2))
traceplot(delta.mu2, main="Chains", xlab="n", ylab="Value")
densplot(delta.mu2, main="Posterior Distribution", xlab="Value", ylab="Density")
# Mixed condition
x11() ; par(mfrow=c(1,2))
traceplot(delta.mu3, main="Chains", xlab="n", ylab="Value")
densplot(delta.mu3, main="Posterior Distribution", xlab="Value", ylab="Density")
x11() ; par(mfrow=c(1,2))
traceplot(delta.si, main="Chains", xlab="n", ylab="Value")
densplot(delta.si, main="Posterior Distribution", xlab="Value", ylab="Density")
