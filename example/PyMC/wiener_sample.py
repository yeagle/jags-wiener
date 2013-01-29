#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# wiener_sample.py
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-01-28
# last mod 2013-01-28 17:30 DW
#

import pymc
import numpy as np
from math import sqrt

def wiener_samplemodel(alpha, tau, beta, delta, N):
  """ easy wiener process model
      to create samples for given parameter values
  """

  @pymc.deterministic
  def alpha(alpha=alpha, N=N):
    return [alpha for i in range(0,N)]
  @pymc.deterministic
  def tau(tau=tau, N=N):
    return [tau for i in range(0,N)]
  @pymc.deterministic
  def beta(beta=beta, N=N):
    return [beta for i in range(0,N)]
  @pymc.deterministic
  def delta(delta=delta, N=N):
    return [delta for i in range(0,N)]
  @pymc.deterministic
  def N(N=N):
    return N

  def Y_logp(value, alpha, tau, beta, delta):
    """ this function is not needed in this model
    """
    return -1

  def Y_rand(alpha, tau, beta, delta, dt=0.0001, sigma=1):
    p = .5 * (1+((delta*sqrt(dt))/sigma))
    i = 0
    y = beta*alpha
    while (y < alpha and y > 0):
      if(pymc.random_number() <= p):
        y = y + sigma*sqrt(dt)
      else:
        y = y - sigma*sqrt(dt)
      i = i+1
    if(y >= alpha):
      t = (i*dt+tau)
    else:
      t = -(i*dt+tau)
    return t

  Y = np.empty(N, dtype=object)
  for i in range(0,N):
    Y[i] = pymc.Stochastic(logp=Y_logp,
                                doc = 'Wiener Model',
                                name = 'Y_%i' % i,
                                parents = {'alpha': alpha[i], 'tau': tau[i],
                                           'beta': beta[i], 'delta': delta[i]},
                                random = Y_rand,
                                trace = True,
                                value = 0,
                                dtype= None,
                                rseed = 1.,
                                observed = False,
                                cache_depth = 2,
                                plot = False,
                                verbose = 0)
  
  return locals()

if __name__ == '__main__':
  # Create model and sample from it
  N=100
  M = pymc.MCMC(wiener_samplemodel(1.0, 0.3, 0.5, 0, N))
  M.sample(iter=1, burn=0, thin=1)

  # the N samples:
  samples = [round(float(M.trace('Y_%i' %i)[:]),4) for i in range(0,N)]
  samples
