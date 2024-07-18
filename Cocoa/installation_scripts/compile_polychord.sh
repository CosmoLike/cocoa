#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_SAMPLER_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi
  
  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
  
  unset_env_vars () {
    unset -v ECODEF FOLDER PACKDIR PLIB
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
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  # ---------------------------------------------------------------------------

  ptop 'COMPILING POLYCHORD' || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER=${POLY_NAME:-"PolyChordLite"}
  
  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  cdfolder "${PACKDIR}" || return 1

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation

  make clean >${OUT1:?} 2>${OUT2:?} || { error "${EC2:?}"; return 1; }
  
  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"

  rm -rf "${PLIB:?}"/pypolychord-*
  rm -rf "${PACKDIR:?}/lib/*.a"
  rm -rf "${PACKDIR:?}/lib/*.so"
  
  # ---------------------------------------------------------------------------

  make >${OUT1:?} 2>${OUT2:?} || { error "${EC7:?}"; return 1; }

  make -j $MNT all \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC7:?}"; return 1; }

  make -j $MNT pypolychord \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

  CC="${MPI_CC_COMPILER:?}" CXX="${MPI_CXX_COMPILER:?}" \
    "${PYTHON3:?}" setup.py install --prefix "${ROOTDIR:?}/.local" \
    >${OUT1:?} 2> ${OUT2:?} || { error "${EC9:?}"; return 1; }

 
  pbottom 'COMPILING POLYCHORD' || return 1

  # ----------------------------------------------------------------------------

  unset_all || return 1;
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------