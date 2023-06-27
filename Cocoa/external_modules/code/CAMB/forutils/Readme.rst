===================
ForUtils
===================
:ForUtils: Various fortran 2003/2008 utility classes.
:Version: 1.0
:Author: Antony Lewis
:Homepage: https://github.com/cmbant/forutils
:Code documentation: https://cosmologist.info/forutils/


Description
============

ForUtils is a package of Fortran classes and convenience functions for
Fortran 2003/2008 programs. The utilities comprise these areas:

* ArrayUtils - Find (minimal/maximal) index of an element in an array.
* FileUtils - Classes for handling file access, python-like loadtxt/savetxt functions
* IniObjects - Read/Write name=value configuration files with inheritance, array and default value support.
* MatrixUtils - Read/Write matrices and interface to some BLAS/LAPACK routines.
* MiscUtils - Utility functions for optional arguments.
* MpiUtils - Wrappers for MPI-routines to compile with(out) MPI library.
* ObjectLists - Lists of arbitrary objects including specializations for vectors.
* RandUtils - Some functions to generate random numbers.
* RangeUtils - Maintain sets of equally spaced intervals, e.g. for integration ranges.
* StringUtils - Utilities for strings, like concat of distinct types a.s.o.

For a class summary see the `class trees <https://cosmologist.info/forutils/classes/_index.html>`_.

Getting Started
================

Clone this git repository::

    $ git clone https://github.com/cmbant/forutils

Compile::

    $ make all

This generates the subdirectories Debug and Release and when an MPI library is
available also DebugMPI and ReleaseMPI. Each directory contains a libforutils.a
archive, which can be used directly for static linking on the compiler command
line by giving the absolute filename::

    forutils/<RTYPE>/libforutils.a

or by specifying it as a library::

    -Lforutils/<RTYPE> -lforutils

Specify one of Debug, DebugMPI, Release, or ReleaseMPI for <RTYPE>. The Debug
release types contain debug symbols and use no optimization, while the Release
types use a reasonable level of optimization.


Dependencies
=============
* Fortran 2008 compatible compiler - E.g., ifort 14+, gfortran 6 or higher.
* MPI library - Only when you want the MpiUtils fully functional. Without an MPI library MpiUtils compile, but the functions are merely no-ops and the makefile target DebugMPI and ReleaseMPI can not be built.
