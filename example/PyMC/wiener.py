#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# wiener.py
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-01-28
# last mod 2013-01-28 17:14 DW
#

import pymc
import numpy as np
from math import pi, log, exp, ceil, floor, sqrt, pow, sin
from sys import float_info

def wiener_model(RT, Cond, N, C, inits=None, WIENER_ERR=0.0001):
  """ easy wiener process model
      inits takes initial values for alpha, tau and delta
  """
  if inits is None:
    inits = (1,0.001,[0 for i in range(0,C)])
  alpha = pymc.Uniform('alpha', lower=0.0001,upper=5, value=inits[0])
  tau = pymc.Uniform('tau', lower=0.0001,upper=1, value=inits[1])

  @pymc.deterministic
  def beta():
    return 0.5

  delta = np.empty(C, dtype=object)
  for i in range(0,C):
    delta[i] =  pymc.Normal('delta_%i' % i, mu=0, tau=1, value=inits[2][i])

  def Y_logp(value, alpha, tau, beta, delta):
    # make sure all variables are float:
    alpha = float(alpha)
    tau = float(tau)
    beta = float(beta)
    delta = float(delta)

    # extract RT
    if (value < 0):
      value = abs(value)
    else:
      beta = 1-beta;
      delta = -delta;
    
    value = value-tau # remove non-decision time from value
    value = value / (pow(alpha,2)) # convert t to normalized time tt

    if (value < 0):
      return -float_info.max

    # calculate number of terms needed for large t
    if (pi*value*WIENER_ERR<1): # if error threshold is set low enough
      kl = sqrt(-2*log(pi*value*WIENER_ERR)/(pow(pi,2)*value)) # bound
      if not(kl>1/(pi*sqrt(value))): # ensure boundary conditions met
        kl = 1/(pi*sqrt(value))
    else: # if error threshold set too high
      kl = 1/(pi*sqrt(value)) # set to boundary condition

    # calculate number of terms needed for small t
    if ((2*sqrt(2*pi*value)*WIENER_ERR)<1): # if error threshold is set low enough
      ks = 2+sqrt(-2*value*log(2*sqrt(2*pi*value)*WIENER_ERR)) # bound
      if not(ks>sqrt(value)+1): # ensure boundary conditions are met 
        ks = sqrt(value)+1
    else: # if error threshold was set too high
      ks = 2 # minimal kappa for that case

    # compute density: f(tt|0,1,beta)
    ans = 0 #initialize density
    if (ks<kl): # if small t is better (i.e., lambda<0)
      K = ceil(ks) # round to smallest integer meeting error
      for k in range(int(-floor((K-1)/2)), int(ceil((K-1)/2)+1)): # loop over k
        ans = ans+(beta+2*k)*exp(-(pow((beta+2*k),2))/2/value) # increment sum
      ans = log(pseudo_zero(ans))-0.5*log(2)-log(sqrt(pi))-1.5*log(value) # add constant term
    else: # if large t is better...
      K = ceil(kl) # round to smallest integer meeting error
      for k in range(0,int(K)):
        ans = ans+k*exp(-(pow(k,2))*(pow(pi,2))*value/2)*sin(k*pi*beta); # increment sum
      ans = log(pseudo_zero(ans))+2*log(sqrt(pi)) # add constant term

    # convert to f(t|delta,a,beta) and return result
    return ans+((-delta*alpha*beta -(pow(delta,2))*(value*pow(alpha,2))/2)-log(pow(alpha,2)))

  def Y_rand(alpha, tau, beta, delta):
    """ this function is not needed in this model
    """
    pass

  Y = np.empty(N, dtype=object)
  for i in range(0,N):
    Y[i] = pymc.Stochastic(logp=Y_logp,
                                doc = 'Wiener Model',
                                name = 'Y_%i' % i,
                                parents = {'alpha': alpha, 'tau': tau,
                                           'beta': beta, 'delta': delta[Cond[i]]},
                                random = Y_rand,
                                trace = True,
                                value = RT[i],
                                dtype= None,
                                rseed = 1.,
                                observed = True,
                                cache_depth = 2,
                                plot = False,
                                verbose = 0)
  
  return locals()

def pseudo_zero(x):
  if (x == 0):
    return float_info.min
  else:
    return x

def read_table(datafile, skip=0):
  """ read txt file containing data table
  """
  with open(datafile, "r") as f:
    rows = [raw.strip().split() for raw in f]
    cols = list(zip(*rows[skip:]))
    return cols

if __name__ == '__main__':
  # read data
  data_raw = read_table('../experiment_example/R/data/exp1_01_m23r_2012-10-31-110958.txt', 9)
  Cond = list()
  for cond in data_raw[1][20:]:
    if cond == 'brigh':
      Cond.append(0)
    elif cond == 'undef':
      Cond.append(1)
    elif cond == 'dark':
      Cond.append(2)
  C = 3
  RT = map(float, data_raw[4][20:])
  for i in range(0,len(RT)):
    if data_raw[2][20:][i] == 'a':
      RT[i] = - float(RT[i])
  N = len(RT)
  
  # Create model and sample from it
  M = pymc.MCMC(wiener_model(RT, Cond, N, C))
  M.sample(iter=5000, burn=0, thin=1)

  # plot results with R
  from rpy import r 
  r.x11()
  r('par(mfrow=c(3,1))')
  r.plot(M.trace('delta_0')[:], col="deepskyblue3", type="l", xlab="", ylab="")
  r.plot(M.trace('delta_1')[:], col="deepskyblue3", type="l", xlab="", ylab="")
  r.plot(M.trace('delta_2')[:], col="deepskyblue3", type="l", xlab="", ylab="")
  r.x11()
  r('par(mfrow=c(3,1))')
  r.plot(r.density(M.trace('delta_0')[:]), col="deepskyblue4", type="l", xlab="", ylab="", xlim=[-3,3])
  r.plot(r.density(M.trace('delta_1')[:]), col="deepskyblue4", type="l", xlab="", ylab="", xlim=[-3,3])
  r.plot(r.density(M.trace('delta_2')[:]), col="deepskyblue4", type="l", xlab="", ylab="", xlim=[-3,3])
