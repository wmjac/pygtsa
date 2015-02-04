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

# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pygtsa], [0.1], [wmjacx@gmail.com])
AC_CONFIG_SRCDIR([bridges.c])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([gnu])
AC_PROG_INSTALL
LT_INIT

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([gsl], [gsl_vector_alloc])
AC_CHECK_LIB([gslcblas], [cblas_ddot])
AC_CHECK_LIB([m], [exp])

# Checks for header files.
AC_CHECK_HEADERS([float.h stddef.h stdlib.h string.h])
AC_CHECK_HEADERS([gsl/gsl_rng.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strchr])

# Search for Python and required packages
AM_PATH_PYTHON([3])
AC_SUBST(PYTHON_INCLUDES, `python3-config --includes`)
AC_PREFIX_DEFAULT(`python3-config --prefix`)
echo -ne "checking for numpy... "
AS_IF([$PYTHON -c 'import numpy'], [echo "yes"], [AC_MSG_FAILURE([cannot import numpy in Python$PYTHON_VERSION])])
echo -ne "checking for networkx... "
AS_IF([$PYTHON -c 'import networkx'], [echo "yes"], [AC_MSG_FAILURE([cannot import networkx in Python$PYTHON_VERSION])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
