#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${ROOTDIR} NOT DEFINED \033[0m"
    return 1
  fi
  
  # ----------------------------------------------------------------------------
  # Clean any previous compilation
  source "${ROOTDIR:?}/installation_scripts/clean_polychord.sh" || return 1
  # ----------------------------------------------------------------------------
  
  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
  
  unset_env_vars () {
    unset -v POLYF PACKDIR
    cdroot || return 1;
  }
  
  error () {
    fail_script_msg "compile_polychord.sh" "${1}"
    unset -f error
    unset_env_vars || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'COMPILING POLYCHORD' || return 1

  unset_env_vars || return 1;

  POLYF=${POLY_NAME:-"PolyChordLite"}
  
  PACKDIR="${ROOTDIR:?}/external_modules/code/${POLYF:?}"

  cdfolder "${PACKDIR}" || return 1

  make -j $MNT all >${OUT1:?} 2>${OUT2:?} || { error "${EC7:?}"; return 1; }

  make -j $MNT pypolychord >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

  CC="${MPI_CC_COMPILER:?}" CXX="${MPI_CXX_COMPILER:?}" \
    ${PYTHON3:?} setup.py install --prefix "${ROOTDIR:?}/.local" \
    >${OUT1:?} 2> ${OUT2:?} || { error "${EC9:?}"; return 1; }

  unset_env_vars || return 1;
  
  pbottom 'COMPILING POLYCHORD' || return 1
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------