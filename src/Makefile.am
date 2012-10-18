SUBDIRS = distributions jwmath 

jagsmod_LTLIBRARIES = wiener.la

wiener_la_SOURCES = wiener.cc

wiener_la_CPPFLAGS = -I$(top_srcdir)/src/include

wiener_la_LIBADD = distributions/wienerdist.la \
        jwmath/libjwmath.la	\
        $(top_builddir)/src/lib/libjags.la \
        @LAPACK_LIBS@ @BLAS_LIBS@ 

wiener_la_LDFLAGS = -module -avoid-version
if WINDOWS
wiener_la_LDFLAGS += -no-undefined
endif
