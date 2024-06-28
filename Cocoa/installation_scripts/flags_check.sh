#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  cdroot
  unset -f error
  return 1;
}

if [ -z "${PYTHON3}" ]; then
  error "PYTHON3" || return 1
fi

if [ -z "${PYTHON_VERSION}" ]; then
  error "PYTHON_VERSION" || return 1
fi

if [ -z "${CXX_COMPILER}" ]; then
  error 'CXX_COMPILER' || return 1
fi

if [ -z "${C_COMPILER}" ]; then
  error 'C_COMPILER' || return 1
fi

if [ -z "${FORTRAN_COMPILER}" ]; then
  error 'FORTRAN_COMPILER' || return 1
fi 

if [ -z "${PIP3}" ]; then
  error 'PIP3' || return 1
fi

if [ -z "${MPI_CXX_COMPILER}" ]; then
  error 'MPI_CXX_COMPILER' || return 1
fi

if [ -z "${MPI_CC_COMPILER}" ]; then
  error 'MPI_CC_COMPILER' || return 1
fi

if [ -z "${CMAKE}" ]; then
  error 'CMAKE' || return 1
fi

if [ -z "${GIT}" ]; then
  error 'GIT' || return 1
fi

unset -f error

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------