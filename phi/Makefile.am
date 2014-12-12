noinst_PROGRAMS = heat heat_mic

CLEANFILES=core.*

COMMON=$(srcdir)/../src
COMMON_FILES=$(COMMON)/main.c $(COMMON)/fail.c $(COMMON)/input.c $(COMMON)/output.c $(COMMON)/img.c

CPPFLAGS= -I $(COMMON)
heat_SOURCES = $(COMMON_FILES) compute.c
heat_CPPFLAGS = -openmp -O3 -vec-report3
heat_LDFLAGS = -openmp

heat_mic_SOURCES = $(COMMON_FILES) compute.c
heat_mic_CPPFLAGS = -openmp -O3 -vec-report3 -mmic
heat_mic_LDFLAGS = -openmp -mmic


#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
#include $(srcdir)/../src/libcommon.mk


