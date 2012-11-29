noinst_LTLIBRARIES = wienerfunc.la

wienerfunc_la_CPPFLAGS = -I$(top_srcdir)/src \
												 -I$(includedir)/JAGS

wienerfunc_la_LDFLAGS = -no-undefined -module -avoid-version

wienerfunc_la_SOURCES = WFunction.cc DFunction.cc DLogFunction.cc

noinst_HEADERS = WFunction.h DFunction.h DLogFunction.h
