#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_EUCLID_EMULATOR_V2_CODE}" ]; then
  
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

  # ---------------------------------------------------------------------------

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${EE2_NAME:-"euclidemu2"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="EUCLID EMULATOR V2"

  ptop "COMPILING ${PRINTNAME:?}" || return 1

  cdfolder "${PACKDIR}" || return 1

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation

  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/euclidemu2.egg-info/"
  
  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"

  rm -rf  "${PLIB:?}"/euclidemu2
  rm -rf  "${PLIB:?}"/euclidemu2-*
  
  # ---------------------------------------------------------------------------  
 
  #WE NEED TO PREVENT ALL COMPILE COMMANDS FROM USING THE INTERNET
  #FROM: https://github.com/pypa/pip/issues/12050
  #So, you're installing setuptools in your current environment. 
  #Then you're installing the current project with --no-index. 
  #But the current project needs setuptools (and wheel) installed. 
  #Since pip 23.1 we don't use setup.py install for legacy projects, 
  #we do a build using an isolated environment. It's that build 
  #(which is a subprocess) that needs setuptools and can't get it 
  #because of the --no-index.
  #You can avoid this by using --no-build-isolation, or you can make 
  #setuptools and wheel available (via --find-links or similar) 
  #for the isolated enviornment creation.

  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
    ${PACKDIR:?} --no-dependencies --prefix="${ROOTDIR:?}/.local" \
    --no-index --no-deps --no-build-isolation \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC3:?}"; return 1; }

  pbottom "COMPILING ${PRINTNAME:?}" || return 1

  # ---------------------------------------------------------------------------

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------