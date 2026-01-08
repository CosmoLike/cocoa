#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FASTPT_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER TFILE TFILEP AL
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile
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

  cpfolder() {
    cp -r "${1:?}" "${2:?}"  \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ---------------------------------------------------------------------------
  
  URL="${FASTPT_URL:?}"

  URL2="${FASTPT_WRAPPER_URL:?}"

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${FASTPT_NAME:-"FAST-PT"}"

  FOLDER2="${FASTPT_WRAPPER_NAME:-"PyFAST-PT"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  PACKDIR2="${ECODEF:?}/${FOLDER2:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="PyFAST-PT"

  ptop "INSTALLING ${PRINTNAME:?}" || return 1;

  # ----------------------------------------------------------------------------
  # In case this script is called twice ----------------------------------------
  # ----------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_FASTPT_CODE}" ]; then
    rm -rf "${PACKDIR:?}";
    rm -rf "${PACKDIR2:?}";
  fi
  
  # FAST-PT
  if [[ ! -d "${PACKDIR:?}" ]]; then
    
    # --------------------------------------------------------------------------
    # Clone from original FAST-PT repo -----------------------------------------
    # --------------------------------------------------------------------------
    cdfolder "${ECODEF:?}" || { cdroot; return 1; }

    "${GIT:?}" clone "${URL:?}" "${FOLDER:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    cdfolder "${PACKDIR}" || { cdroot; return 1; }

    if [ -n "${FASTPT_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${FASTPT_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
      
  fi

  # FAST-PT cobaya theory wrapper
  if [[ ! -d "${PACKDIR2:?}" ]]; then
    
    # --------------------------------------------------------------------------
    # Clone from original FAST-PT theory wrapper repo --------------------------
    # --------------------------------------------------------------------------
    cdfolder "${ECODEF:?}" || { cdroot; return 1; }

    "${GIT:?}" clone "${URL2:?}" "${FOLDER2:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    cdfolder "${PACKDIR2}" || { cdroot; return 1; }

    if [ -n "${FASTPT_WRAPPER_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${FASTPT_WRAPPER_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

  fi

  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "INSTALLING ${PRINTNAME:?}" || return 1
  
  # ----------------------------------------------------------------------------

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
