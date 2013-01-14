SUBDIRS = distributions functions

jagsmod_LTLIBRARIES = wiener.la

wiener_la_SOURCES = wiener.cc

wiener_la_CPPFLAGS = -I$(includedir)/JAGS

wiener_la_LIBADD = distributions/wienerdist.la \
									 functions/wienerfunc.la
if WINDOWS
wiener_la_LIBADD += -ljags-$(JAGS_MAJOR) -ljrmath-0
else
wiener_la_LIBADD += -ljags -ljrmath
endif

wiener_la_LDFLAGS = -module -avoid-version
if WINDOWS
wiener_la_LDFLAGS += -no-undefined
endif
