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
    Wabersich, D. and Vandekerckhove, J. (in preparation). Extending JAGS: 
    A tutorial on adding custom distributions to JAGS (with a diffusion
    model example)

Also note, that the tutorial was written for JAGS 3.3.0. The module code
has been updated to work with JAGS 4.2.0.

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

- Install MinGW (MinGW installer including msys):
  http://www.mingw.org
  In MinGW get installer: msys-base is all you need.

- Install TDM-GCC Compiler Suite:
  http://tdm-gcc.tdragon.net 
  Use the exact same version of tdm-gcc as has been used to compile your
  JAGS Version on your system. In that case, also don't forget to untick
  the  "Check for updated files on the TDM-GCC server" 
  during the installation proccess of tdm-gcc. 
  For Jags-4.0.0 and Jags-4.1.0 use: 
  tdm64-gcc-4.6.1.exe
  For Jags-4.2.0 use: 
  tdm64-gcc-4.9.2-3.exe
  
- Delete all *.dll.a files int the TDM-GCC Installation, so the
  compiler will link to the static libraries.

- Change path in c:/mingw/msys/1.0/etc/fstab from
  this:     c:/mingw /mingw
  to:       c:/MinGW64 /mingw
  to use TDM-GCC compilers

- Alternatively to TDM-GCC, you could also use Rtools
  (https://cran.r-project.org/bin/windows/Rtools/).
  For Jags-4.2.0 use:
  Rtools33.exe

  Change path in c:/mingw/msys/1.0/etc/fstab
  to:       C:/Rtools/mingw_64 /mingw
  to use Rtools compilers

- --> Start msys, extract tarball in your home dir, 
      cd into dir and do the following:

- *Note*: the module needs the JAGS include files and
  the JAGS libraries.
  Therefore use the -L and -I options, as shown below.

- For building 32bit binaries:

.. code:: sh

  CXX="g++ -m32" \
  ./configure LDFLAGS="-L/c/Progra~1/JAGS/JAGS-4.2.0/i386/bin" CXXFLAGS="-I/c/Progra~1/JAGS/JAGS-4.2.0/include"

  make win32-install

- For building 64bit binaries:

.. code:: sh

  CXX="g++ -m64" \
  ./configure LDFLAGS="-L/c/Progra~1/JAGS/JAGS-4.2.0/x64/bin" CXXFLAGS="-I/c/Progra~1/JAGS/JAGS-4.2.0/include"

  make win64-install

- Copy the win/win32/wiener.* or win/win64/wiener.* 
  files to your JAGS modules directory to enable the module.
  For JAGS-4.2.0 this usually is: 
  x64: C:\\Program Files\\JAGS\\JAGS-4.2.0\\x64\\modules
  x32: C:\\Program Files\\JAGS\\JAGS-4.2.0\\x32\\modules

- *Note*: If you compile both 32bit and 64bit don't forget to make clean
  between the building.

- *All Windows commands in one* (for copy paste convenience)

.. code:: sh

  CXX="g++ -m32" \
  ./configure LDFLAGS="-L/c/Progra~1/JAGS/JAGS-4.2.0/i386/bin" \
  CXXFLAGS="-I/c/Progra~1/JAGS/JAGS-4.2.0/include" && \
  make win32-install && \
  make clean && \
  CXX="g++ -m64" \
  ./configure LDFLAGS="-L/c/Progra~1/JAGS/JAGS-4.2.0/x64/bin" \
  CXXFLAGS="-I/c/Progra~1/JAGS/JAGS-4.2.0/include" &&
  make win64-install && \
  make clean
