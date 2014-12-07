noinst_PROGRAMS = heat

heat_SOURCES = $(COMMON)/main.c $(COMMON)/input.c $(COMMON)/fail.c compute.cu
heat_LINK = nvcc -g

.cu.o:
	nvcc -arch sm_21 -g -c -o $@ -I ../src $<

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force
include $(srcdir)/../src/libcommon.mk


# to allow build even when the reference is
# not available:
DISTCLEANFILES = ref1.c ref2.c
BUILT_SOURCES = ref1.c ref2.c
ref1.c ref2.c:
	touch $@
