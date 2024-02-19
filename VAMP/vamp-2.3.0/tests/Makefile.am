# Makefile.am --
# $Id: Makefile.am 319 2010-04-18 15:00:54Z ohl $
########################################################################

AM_FCFLAGS = -I$(top_builddir)/src
AM_LDFLAGS =

KINDS = $(top_builddir)/src/kinds.lo
LIBVAMP = $(top_builddir)/src/libvamp.la

TESTS = tao_test
EXTRA_PROGRAMS = tao_test

tao_test_SOURCES = tao_test.f90
tao_test_LDADD = $(KINDS) $(LIBVAMP)

TESTS += vamp_test
EXTRA_PROGRAMS += vamp_test

vamp_test_SOURCES = vamp_test.f90 constants.f90 kinematics.f90 coordinates.f90
vamp_test_LDADD = $(KINDS) $(LIBVAMP)

TESTS += vamp_test0
EXTRA_PROGRAMS += vamp_test0

vamp_test0_SOURCES = vamp_test0.f90
vamp_test0_LDADD = $(KINDS) $(LIBVAMP)

VPATH = $(srcdir):$(top_builddir)/src:$(top_srcdir)/src

SOURCE_FILES = \
    tao_test.f90 vamp_test0.f90 vamp_test.f90 \
    constants.f90 kinematics.f90 coordinates.f90

# Fortran90 module files are generated at the same time as object files
.$(OBJEXT).$(FC_MODULE_EXT):
	@:
#	touch $@

########################################################################
# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
########################################################################
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: $(SOURCE_FILES)
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.$(OBJEXT): /' \
	          -e 's/, *only:.*//' \
	          -e 's/, *&//' \
	          -e 's/, *.*=>.*//' \
	          -e 's/ *$$/.$(OBJEXT)/' ; \
	done > $@

DISTCLEANFILES = Makefile.depend

########################################################################
# noweb
########################################################################

NOTANGLE_IT = \
  $(NOTANGLE) -R'[[$(@F)]]' $(top_srcdir)/src/prelude.nw $^ \
    | $(CPIF) $@

SUFFIXES = .nw .f90

.nw.f90:
	$(NOTANGLE_IT)

tao_test.f90: tao_random_numbers.nw
	$(NOTANGLE_IT)

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	rm -f *.$(FC_MODULE_EXT) *.g90
	rm -f *.mod tmp.tao *.d *.grid *.grids
if FC_SUBMODULES
	-rm -f *.smod
endif

########################################################################
# The End.
########################################################################
