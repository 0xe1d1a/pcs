noinst_PROGRAMS = heat barrier pipesort locking 

heat_SOURCES = $(COMMON)/main.c compute.c
barrier_SOURCES = barrier.c
pipesort_SOURCES = pipesort.c
locking_SOURCES = locking.c
locking_CFLAGS = -fgnu-tm

## Uncomment the following:
# noinst_PROGRAMS += barrier_test
# barrier_test_SOURCES = ...

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


