dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)

AC_INIT([JW],[1.0.0],[dominik.wabersich@gmail.com],[JAGS-WIENER-MODULE])
JW_MAJOR=4
JW_MINOR=0
AC_SUBST(JW_MAJOR)
AC_SUBST(JW_MINOR)

AC_CANONICAL_HOST
AC_CONFIG_SRCDIR([src])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

dnl libtool and ltdl configuration
LT_PREREQ(2.2.6)
LT_CONFIG_LTDL_DIR([libltdl])
LT_INIT([dlopen disable-static win32-dll])
LTDL_INIT([recursive])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AM_PROG_LEX
AC_PROG_YACC

dnl IEEE 754 arithmetic
AC_CHECK_HEADERS(ieeefp.h)
R_IEEE_754
if test "${ac_cv_func_finite}" = "yes"; then
R_FUNC_FINITE
fi
if test "${ac_cv_have_decl_isfinite}" = "yes"; then
R_FUNC_ISFINITE
fi

case "${host_os}" in
  mingw*)
    win=true ;;
  *)
    win=false ;;
esac
AM_CONDITIONAL(WINDOWS, test x$win = xtrue)

dnl fortran stuff
AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS
ACX_LAPACK([],[
	AC_MSG_ERROR("You need to install the LAPACK library")
])

dnl Rmath stuff

AC_DEFINE(MATHLIB_STANDALONE, 1, [Define if you have standalone R math library])

dnl Relevant parts lifted from R's configure.ac, version 2.3.0
AC_CHECK_LIB(m, sin)

## IEEE 754. We rely on this in e.g. the working log test.
##R_IEEE_754

R_CHECK_FUNCS([expm1 log1p], [#include <math.h>])
## Used to build src/include/Rmath.h.
## <NOTE>
## we don't use AC_CONFIG_HEADERS on Rmath.h.in because
## a) that would comment out #undef statements in Rmath.h.in and
## b) Rmath.h should be a self-contained file for standalone Rmath use.
## </NOTE>
if test "${ac_cv_have_decl_expm1}" = yes; then
  RMATH_HAVE_EXPM1="# define HAVE_EXPM1 1"
else
  RMATH_HAVE_EXPM1="# undef HAVE_EXPM1"
fi
AC_SUBST(RMATH_HAVE_EXPM1)
if test "${ac_cv_have_decl_log1p}" = yes; then
  RMATH_HAVE_LOG1P="# define HAVE_LOG1P 1"
else
  RMATH_HAVE_LOG1P="# undef HAVE_LOG1P"
fi
AC_SUBST(RMATH_HAVE_LOG1P)

R_FUNC_LOG
R_FUNC_LOG1P

dnl End Rmath stuff

jagsmoddir=${libdir}/JAGS/modules-${JAGS_MAJOR}
AC_SUBST(jagsmoddir)

AC_CONFIG_FILES([
  Makefile
  m4/Makefile
  libltdl/Makefile
  src/Makefile
  src/distributions/Makefile	
  src/jwmath/Makefile	
])
AC_OUTPUT
