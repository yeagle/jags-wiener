JAGS wiener module
==================
The JAGS Wiener module is an extension for JAGS, which provides wiener
process distribution functions, mainly the Wiener first passage time
density. It allows to include stochastic nodes with the first hitting time
distribution of a diffusion process.

Example
-------
In a model file, you can use:

.. code:: R

  x ~ dwiener(alpha,tau,beta,delta)
 
With the parameters:
  - alpha being the boundary separation parameter
  - tau being the non-decision time
  - beta being the bias
  - delta being the driftrate

There is also a second distribution:

.. code:: R

  x ~ dwieners(alpha,tau,beta,delta,1)

with the last parameter as standard deviation of the drift process.

And there are two logical nodes:

  - dwiener(x, alpha,tau,beta,delta)
  - dlogwiener(x, alpha,tau,beta,delta)

to calculate the (log-) density at value x.

Please note
-----------
Copyright (C) 2012 Dominik Wabersich <dominik.wabersich@gmail.com>
and Joachim Vandekerckhove <joachim@uci.edu>

When using this module, please cite as: 
    Wabersich, D. and Vandekerckhove, J. (2014). Extending JAGS: 
    A tutorial on adding custom distributions to JAGS (with a diffusion
    model example)

Also note, that the tutorial was written for JAGS 3.3.0. The module code
has been updated to work with JAGS 4.2.0 and JAGS 4.3.1.

Known Issues
------------
* Loading of the JAGS Wiener Module does not work (so far reported only by mac users).
  4 reasons for this issue are possible:

  - In MAC OS X the Versions of JAGS and the JAGS Wiener Module have been
    compiled for different MAC Versions (e.g. Leopard and Mavericks)
  - File permissions of the library files are wrong.
  - The library files are installed in a different location than the JAGS
    installation.
  - The library files have been compiled with a different gcc version than
    the JAGS installation.
  
* The easiest way to solve this issue, is by removing both, the JAGS and
  JAGS Wiener Module installation. Afterwards, recompile and install first JAGS,
  then the JAGS Wiener Module. If the same gcc version and the same
  installation procedure have been used, it should work.

For more Information on compiling also check the JAGS installation manual.

License
-------
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundationa Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

Installation
------------

Linux and Mac
"""""""""""""
Use the tarball and install as usual: 

.. code:: sh

  ./configure && make && sudo make install

Note: It might be necessary to define a different prefix, depending on
where JAGS is copied (i.e. installed), for example:

.. code:: sh

  ./configure --prefix /usr && make && sudo make install

Mac
"""
See additional instructions here: 

https://github.com/kiante-fernandez/Rhddmjags/blob/main/jags_wiener_install.md

Windows
"""""""
For Windows we provide precompiled binaries, which come with an
installer.

Compiling from a cloned repository
----------------------------------
.. code:: sh

  # dependencies (on a clean ubuntu installation)
  sudo apt-get install autoconf automake libtool g++

  # creating all auxiliary files
  autoreconf -fvi

  # building
  ./configure
  make

  # or, if JAGS has been installed in a different location, e.g. /usr
  ./configure --prefix /usr
  make

  # install
  sudo make install

Windows Compiling 
"""""""""""""""""

**First, in Linux:**

- For building the module in Windows, it is easiest to use a tarball that
  was created in linux like this (starting from a source clone):

.. code:: sh

    autoreconf -fvi
    ./configure
    make dist-gzip

- As an alternative to building this tarball yourself from the github
  source, one can use the tarball that is available for the latest release.

- Copy the *.tar.gz file to your msys home directory and continue from
  there.

**Second, in Windows:**

- Use Rtools
  (https://cran.r-project.org/bin/windows/Rtools/).
  For Jags-4.2.0 use: Rtools33.exe
  For Jags-4.3.1 use: Rtools42.exe

- Start mingw included in Rtools, extract tarball in your home dir, cd into dir and do the following:

.. code:: sh

  export PATH=/x86_64-w64-mingw32.static.posix/bin:$PATH

  ./configure --host=x86_64-w64-mingw32.static.posix \
  LDFLAGS="-L/c/Progra~1/JAGS/JAGS-4.3.1/x64/bin" \
  CXXFLAGS="-I/c/Progra~1/JAGS/JAGS-4.3.1/include" 
  make win64-install

- Copy the win/win64/wiener.* 
  files to your JAGS modules directory to enable the module.
  For JAGS-4.3.1 this usually is: 
  x64: C:\\Program Files\\JAGS\\JAGS-4.3.1\\x64\\modules


- For the installer install NSIS 3.09 and do the following:

.. code:: sh

  export PATH=$PATH:/c/Program\ Files\ \(x86\)/NSIS
  make installer

- *All Windows commands in one* (for copy paste convenience)

.. code:: sh

  export PATH=/x86_64-w64-mingw32.static.posix/bin:$PATH
  export PATH=$PATH:/c/Program\ Files\ \(x86\)/NSIS
  ./configure --host=x86_64-w64-mingw32.static.posix \
  LDFLAGS="-L/c/Progra~1/JAGS/JAGS-4.3.1/x64/bin" \
  CXXFLAGS="-I/c/Progra~1/JAGS/JAGS-4.3.1/include" && \
  make win64-install && \
  make installer && \
  make clean
