# Install MinGW (MinGW installer including msys) 
# and TDM-GCC Compiler Suite
# http://www.mingw.org
# http://tdm-gcc.tdragon.net 
# change path in c:/mingw/msys/1.0/etc/fstab from
# this:     c:/mingw /mingw
# to:       c:/MinGW64 /mingw
# to use TDM-GCC compilers

# Put the JAGS libs and include files in /local 
# on your Windows msys (or change paths accordingly)
# Also before building: Edit the Makefile in src, to link to the correct
# ljags-x and ljrmath-x
# Otherwise you'll probably get a linker error, 
# telling you that it cannot find ljags and ljrmath


# Start msys, extract tarball in your home dir, 
# cd into dir and do the following:

##############################################

# For building 32bit binaries
CXX="g++ -m32" \
./configure LDFLAGS="-L/local/lib -Wl,--enable-auto-import -I/local/include/JAGS"

# For building 64bit binaries
CXX="g++ -m64" \
./configure LDFLAGS="-L/local/lib -Wl,--enable-auto-import -I/local/include/JAGS"

# Copy the win/instxx/lib/modules-x/wiener.* files to your JAGS modules dir
# to enable the module