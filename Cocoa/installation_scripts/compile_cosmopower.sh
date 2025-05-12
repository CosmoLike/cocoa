#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_COSMOPOWER_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi
    
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

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------  
  # ---------------------------------------------------------------------------  
  
  unset_env_vars || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${COSMOPOWER_NAME:-"cosmopower"}"

  PACKDIR="${ECODEF:?}/emulators/${FOLDER:?}"

  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"

  ptop "COMPILING COSMOPOWER" || return 1

  # ---------------------------------------------------------------------------
  #rm -rf "${PACKDIR:?}/build/"
  #rm -rf "${PACKDIR:?}/euclidemu2.egg-info/"
  #rm -rf  "${PLIB:?}"/euclidemu2
  #rm -rf  "${PLIB:?}"/euclidemu2-*
  # ---------------------------------------------------------------------------  
 
  cdfolder "${PACKDIR}" || return 1

  #prevent all compile_XXX.sh from using the internet (run @compute nodes)
  #FROM: https://github.com/pypa/pip/issues/12050
  #That is why we use --no-dependencies --no-index --no-build-isolation
  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
   env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
   ${PACKDIR:?} --no-dependencies --prefix="${ROOTDIR:?}/.local" --no-index \
   --no-build-isolation >${OUT1:?} 2>${OUT2:?} || { error "${EC3:?}"; return 1; })

  pbottom "COMPILING COSMOPOWER" || return 1

  cdfolder "${ROOTDIR}" || return 1;

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------