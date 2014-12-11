noinst_PROGRAMS = heat

CLEANFILES=core.*

COMMON=$(srcdir)/../src
COMMON_FILES=$(COMMON)/main.c $(COMMON)/fail.c $(COMMON)/input.c $(COMMON)/output.c $(COMMON)/img.c

CPPFLAGS= -I $(COMMON)
heat_SOURCES = $(COMMON_FILES) compute.c
heat_CPPFLAGS = -openmp -O3 -vec-report3
heat_LDFLAGS = -openmp

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
#include $(srcdir)/../src/libcommon.mk


