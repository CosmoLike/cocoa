#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  if [ -z "${ROOTDIR}" ]; then
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${ROOTDIR} IS NOT DEFINED \033[0m"
    return 1
  fi
  # ---------------------------------------------------------------------------
  source "${ROOTDIR}/installation_scripts/clean_act.sh"
  # ---------------------------------------------------------------------------
  pfail() {
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'; cdroot; return 1;
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1;
  fi
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'; cdroot; return 1;
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"; cdroot; return 1;
  fi
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    pfail 'MAKE_NUM_THREADS'; cdroot; return 1;
  fi
  if [ -z "${ACT_NAME}" ]; then
    pfail 'ACT_NAME'; cdroot; return 1;
  fi
  unset_env_vars_comp_act () {
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_comp_act
    cdroot || return 1;
  }
  fail_comp_act () {
    local MSG="\033[0;31m\t\t (compile_act.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_comp_act
    unset_env_vars_comp_act
  }
  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_comp_act "CD FOLDER: ${1}"; return 1; }
  }
  if [ -z "${DEBUG_ACT_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------  
  ptop 'COMPILING ACT'

  cdfolder "${ROOTDIR}/external_modules/code/${ACT_NAME}"/ || return 1
 
  $PIP3 install . --prefix="${ROOTDIR}"/.local >${OUT1} 2>${OUT2} 
  if [ $? -ne 0 ]; then
    fail_comp_act "PIP3 INSTALL ."; return 1
  fi

  unset_env_vars_comp_act || return 1
  pbottom 'COMPILING ACT'
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------