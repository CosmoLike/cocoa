#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  
  pfail() {
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${1} NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cdroot; return 1;
  fi
  
  if [ -z "${ACT_NAME}" ]; then
    pfail 'ACT_NAME'
    cdroot; return 1;
  fi
  
  unset_env_vars_clean_act () {
    unset pfail
    unset OUT1
    unset OUT2
    unset cdfolder
    unset unset_env_vars_clean_act
    cdroot || return 1;
  }
  
  fail_clean_act () {
    local MSG="\033[0;31m\t\t (clean_act.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_clean_act
    unset_env_vars_clean_act
  }
  
  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_clean_act "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${DEBUG_ACT_OUTPUT}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------------------- 
  
  ptop 'CLEANING ACT'

  cdfolder "${ROOTDIR}/external_modules/code/${ACT_NAME}"/ || return 1

  rm -rf ./build/
  rm -rf ./pyactlike.egg-info/

  $PYTHON3 setup.py clean >${OUT1} 2>${OUT2} ||
    { fail_clean_act "PYTHON SETUP CLEAN"; return 1; }

  unset_env_vars_clean_act || return 1
  
  pbottom 'CLEANING ACT'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------