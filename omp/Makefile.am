noinst_PROGRAMS = heat msort vecsort

heat_SOURCES = $(COMMON)/main.c compute.c

msort_SOURCES = merge.c
msort_LDADD = 

vecsort_SOURCES = vec.c
vecsort_LDADD = 

## Uncomment the following:
# noinst_PROGRAMS += matrix
# matrix_SOURCES = ...

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


