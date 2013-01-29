Alternative to the JAGS wiener module: PyMC
============================================

This is a short comparison, without being exhaustive, for people interested
in Python and/or a more flexible framework for MCMC sampling.


wiener.py
---------

The example shows, how to use the first hitting time distribution of a
diffusion process in Python, with the PyMC package.

The PyMC package offers MCMC functions, like JAGS, but with high
flexibility: As shown in the example, implementing the first hitting time
distribution of a wiener process is fairly simple (see the PyMC
documentation for details).

The generated samples from the PyMC package probably need more
burning and thinning, due to a less optimized sampler algorithm (at least,
at the time of writing this README). This could, however, be addressed by
writing a better sampler (as doing this in the PyMC framework should be
relatively simple).

Using the wiener module in JAGS is about 10 times faster for this easy
example. Furthermore, JAGS uses a slice sampler, which creates better
results.
Also note, that JAGS uses the well-known BUGS code for the model
definition, whereas the Python package uses general Python syntax.


wiener_sample.py
----------------

This simple example shows how to write a sampling function with PyMC.
