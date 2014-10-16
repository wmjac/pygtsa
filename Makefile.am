# Copyright (C) 2014 William M. Jacobs

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libgraph.la
pkgpyexec_LTLIBRARIES = cgraph.la
pkgpython_PYTHON = __init__.py calc_subgraphs.py histogram.py structure.py

noinst_HEADERS = graph.h random.h sample.h
dist_noinst_DATA = cgraph.pxd cgraph.pyx examples/

libgraph_la_SOURCES = graph.c bridges.c adjacents.c subgraphs.c random.c histogram_ev.c
libgraph_la_LIBADD = -lgsl -lgslcblas -lm

cgraph_la_SOURCES = cgraph.c
cgraph_la_CPPFLAGS = -I$(PYTHON_PREFIX)/include/python$(PYTHON_VERSION)
cgraph_la_LDFLAGS = -module
cgraph_la_LIBADD = libgraph.la
