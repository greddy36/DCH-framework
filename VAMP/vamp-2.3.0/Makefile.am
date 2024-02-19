# Makefile.am --
# $Id: Makefile.am 319 2010-04-18 15:00:54Z ohl $
########################################################################
##
## Process this file with automake to produce Makefile.in
#
# Copyright (C) 1999-2016 by 
#     Wolfgang Kilian <kilian@physik.uni-siegen.de>
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
#     Juergen Reuter <juergen.reuter@desy.de>
#     Christian Speckner <cnspeckn@googlemail.com>
#
# WHIZARD is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by 
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# WHIZARD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
########################################################################

SUBDIRS = src share tests

## Use the same compiler and compiler flags for distcheck
DISTCHECK_CONFIGURE_FLAGS = FC=$(FC) FCFLAGS="$(FCFLAGS)" LIBS="$(LIBS)" \
	--enable-distribution

## Use the m4 directory for local Autoconf macros
ACLOCAL_AMFLAGS = -I m4

DISTCLEANFILES = \
  iso_c_binding_gfortran.$(FC_MODULE_EXT) \
  iso_fortran_env_2008.$(FC_MODULE_EXT) \
  query_kinds.$(FC_MODULE_EXT) \
  report_kinds.$(FC_MODULE_EXT)
