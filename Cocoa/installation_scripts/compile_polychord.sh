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
  source "${ROOTDIR:?}/installation_scripts/clean_polychord.sh"
  # ----------------------------------------------------------------------------
  
  pfail() {
    echo -e \
    "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
    unset pfail
  }
  
  cdroot() {
    cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  
  if [ -z "${MPI_CXX_COMPILER}" ]; then
    pfail 'MPI_CXX_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${MPI_CC_COMPILER}" ]; then
    pfail 'MPI_CC_COMPILER'; cdroot; return 1;
  fi
    
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"; cdroot; return 1;
  fi
  
  unset_env_vars_comp_poly () {
    unset OUT1
    unset OUT2
    unset PMNT
    unset pfail
    unset POLYF
    unset unset_env_vars_comp_poly
    cdroot || return 1;
  }
  
  fail_comp_poly () {
    local MSG="\033[0;31m\t\t (compile_polychord.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_comp_poly
    unset_env_vars_comp_poly 
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export PMNT="${MAKE_NUM_THREADS:-1}"; 
    [[ ${PMNT} == +([0-9]) ]] || export PMNT=1
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export PMNT=1
  fi

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_comp_poly "CD FOLDER: ${1}"; return 1; }
  }
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'COMPILING POLYCHORD'

  export POLYF=${POLY_NAME:-"PolyChordLite"}
  export PACKDIR="${ROOTDIR:?}/external_modules/code/${POLYF:?}"

  cdfolder "${PACKDIR}" || return 1

  make all >${OUT1:?} 2>${OUT2:?} || { fail_comp_poly "MAKE ALL"; return 1; }

  make -j $PMNT pypolychord >${OUT1:?} 2>${OUT2:?} || 
    { fail_comp_poly "MAKE PYPOLYCHORD"; return 1; }

  CC="${MPI_CC_COMPILER:?}" CXX="${MPI_CXX_COMPILER:?}" "${PYTHON3:?}" \
    setup.py install --prefix "${ROOTDIR:?}/.local"  > ${OUT1} 2> ${OUT2} || 
    { fail_comp_poly "PYTHON3 SETUP INSTALL"; return 1; }

  unset_env_vars_comp_poly || return 1;
  
  pbottom 'COMPILING POLYCHORD'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------