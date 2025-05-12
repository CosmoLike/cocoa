#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell  
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v ECODEF FOLDER PACKDIR PRINTNAME PLIB
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error
    unset -f unset_env_funcs
    cdroot || return 1;
  }

  unset_all () {
    unset_env_vars
    unset_env_funcs
    unset -f unset_all
    cdroot || return 1;
  }
  
  error () {
    fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }
  
  unset_env_vars || return 1

  # ---------------------------------------------------------------------------
  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${CAMB_NAME:-"CAMB"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  ptop "COMPILING CAMB" || return 1

  cdfolder "${PACKDIR:?}" || return 1

  # ---------------------------------------------------------------------------  
  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/camb/__pycache__/"
  rm -f  "${PACKDIR:?}/camb/camblib.so"
  rm -rf "${PACKDIR:?}/forutils/Releaselib/"
  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  rm -rf  "${PLIB:?}/camb"
  rm -rf  "${PLIB:?}/camb"-*
  # ---------------------------------------------------------------------------
  
  if [ -z "${IGNORE_COSMOREC_CODE}" ] && [ -n "${IGNORE_HYREC_CODE}" ]; then

    (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
     export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
     RECOMBINATION_FILES="recfast cosmorec" COMPILER="${FORTRAN_COMPILER:?}" \
     F90C="${FORTRAN_COMPILER:?}" "${PYTHON3:?}" setup.py build \
     >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; })
  
  elif [ -n "${IGNORE_COSMOREC_CODE}" ] && [ -z "${IGNORE_HYREC_CODE}" ]; then

    (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
     export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
     RECOMBINATION_FILES="recfast hyrec" COMPILER="${FORTRAN_COMPILER:?}" \
     F90C="${FORTRAN_COMPILER:?}" "${PYTHON3:?}" setup.py build \
     >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; })
  
  elif [ -z "${IGNORE_COSMOREC_CODE}" ] && [ -z "${IGNORE_HYREC_CODE}" ]; then

    (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
     export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
     RECOMBINATION_FILES="recfast cosmorec hyrec" \
     COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
     "${PYTHON3:?}" setup.py build >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; })

  else

    (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
     export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
     COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
     "${PYTHON3:?}" setup.py build >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; })

  fi

  #prevent all compile_XXX.sh from using the internet (run @compute nodes)
  #FROM: https://github.com/pypa/pip/issues/12050
  #That is why we use --no-dependencies --no-index --no-build-isolation
  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
   env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
   ${PACKDIR:?} --no-dependencies --prefix="${ROOTDIR:?}/.local" --no-index \
   --no-build-isolation >${OUT1:?} 2>${OUT2:?} || { error "${EC3:?}"; return 1; })

  pbottom "COMPILING CAMB" || return 1

  cdfolder "${ROOTDIR}" || return 1;
  
  unset_all || return 1

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------