SUBDIRS = distributions

jagsmod_LTLIBRARIES = wiener.la

wiener_la_SOURCES = wiener.cc

wiener_la_CPPFLAGS = -I$(includedir)/JAGS

wiener_la_LIBADD = distributions/wienerdist.la \
									 -ljags \
									 @LAPACK_LIBS@ @BLAS_LIBS@ 


wiener_la_LDFLAGS = -module -avoid-version
if WINDOWS
wiener_la_LDFLAGS += -no-undefined
endif
