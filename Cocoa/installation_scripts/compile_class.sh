#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${ROOTDIR} NOT DEFINED \033[0m"
    return 1
  fi
  
  # ---------------------------------------------------------------------------
  # Clean any previous compilation
  source "${ROOTDIR}/installation_scripts/clean_class.sh"
  # ---------------------------------------------------------------------------
  
  pfail() {
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${1} NOT DEFINED \033[0m"
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
  
  if [ -z "${CLASS_NAME}" ]; then
    pfail 'CLASS_NAME'; cdroot; return 1
  fi
  
  unset_env_vars_comp_class () {
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_comp_class
    cdroot || return 1;
  }
  
  fail_comp_class () {
    local MSG="\033[0;31m\t\t (compile_class.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_comp_class
    unset_env_vars_comp_class
  }
  
  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_comp_class "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'COMPILING CLASS'

  cdfolder "${ROOTDIR}/external_modules/code/${CLASS_NAME}/" || return 1

  # ---------------------------------------------------------------------------
  # historical: workaround Cocoa .gitignore entry on /include 
  # also good because class compilation changes files on /include
  # in case this script is called twice
  # ---------------------------------------------------------------------------  
  rm -rf ./include

  # ---------------------------------------------------------------------------
  # historical: workaround Cocoa .gitignore entry on /include
  # also good because class compilation changes files on /include
  # ---------------------------------------------------------------------------
  cp -r  ./include2 ./include >${OUT1} \ 
    2>${OUT2} || { fail_comp_class "CP INCLUDE2 FOLDER"; return 1; }


  CC=$C_COMPILER PYTHON=$PYTHON3 make all >${OUT1} \
    2>${OUT2} || { fail_comp_class "MAKE ALL"; return 1; }
   
  cdfolder ./python || return 1

  CC=$C_COMPILER $PYTHON3 setup.py build >${OUT1} \
    2>${OUT2} || { fail_comp_class "PYTHON3 SETUP.PY BUILD"; return 1; }

  unset_env_vars_comp_class || return 1

  pbottom 'COMPILING CLASS'

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------