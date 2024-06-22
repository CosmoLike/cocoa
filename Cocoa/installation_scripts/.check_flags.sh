#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${PYTHON3}" ]; then
  pfail "PYTHON3"; cdroot; return 1;
fi

if [ -z "${PYTHON_VERSION}" ]; then
  pfail "PYTHON_VERSION"; cdroot; return 1
fi

if [ -z "${CXX_COMPILER}" ]; then
  pfail 'CXX_COMPILER'; cdroot; return 1;
fi

if [ -z "${C_COMPILER}" ]; then
  pfail 'C_COMPILER'; cdroot; return 1;
fi

if [ -z "${PIP3}" ]; then
  pfail 'PIP3'; cdroot; return 1;
fi

if [ -z "${MPI_CXX_COMPILER}" ]; then
  pfail 'MPI_CXX_COMPILER'; cdroot; return 1;
fi

if [ -z "${MPI_CC_COMPILER}" ]; then
  pfail 'MPI_CC_COMPILER'; cdroot; return 1;
fi

if [ -z "${FORTRAN_COMPILER}" ]; then
  pfail 'FORTRAN_COMPILER'; cdroot; return 1;
fi 

if [ -z "${CMAKE}" ]; then
  pfail 'CMAKE'; cdroot; return 1;
fi

if [ -z "${GIT}" ]; then
  pfail 'GIT'; cdroot; return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------