#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v COB CCCOB COBLIKE URL FOLDER PACKDIR PRINTNAME ECODEF URL 
    unset -v COBTH
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile cppatch cppatchfolder
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

  cpfile() {
    cp "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

  COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

  CCCOB="${CCIL:?}/cobaya_changes" # CC = CoCoA, COB = Cobaya (Cocoa Cobaya)
  
  COBLIKE="cobaya/likelihoods"     # COB = Cobaya, LIKE = likelihoods
  
  COBTH="cobaya/theories"          # COB = Cobaya, TH = theories
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop "SETUP SIMONS OBSERVATORY SYSLIBRARY" || return 1;

  URL="${SO_SYSLIB_URL:-"https://github.com/simonsobs/syslibrary.git"}"
    
  FOLDER="${SO_SYSLIB_NAME:-"syslibrary"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  if [[ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE}" ]]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [[ ! -d "${PACKDIR:?}" ]]; then
    cdfolder "${ECODEF}" || return 1;

    "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:?} \
      --recursive "${FOLDER:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

    cdfolder "${PACKDIR}" || return 1;

    if [ -n "${SO_SYSLIB_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${SO_SYSLIB_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  cdfolder "${ROOTDIR}" || return 1

  pbottom "SETUP SIMONS OBSERVATORY SYSLIBRARY" || return 1;
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop "SETUP SIMONS OBSERVATORY MFLIKE" || return 1;

  URL="${SO_MFLIKE_URL:-"https://github.com/simonsobs/LAT_MFLike.git"}"
  
  FOLDER="mflike"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  if [[ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE}" ]]; then    
    rm -rf "${PACKDIR:?}"
  fi

  if [[ ! -d "${PACKDIR:?}" ]]; then
    cdfolder "${ECODEF}" || return 1;

    "${GIT:?}" clone "${URL:?}" --recursive "${FOLDER:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

    cdfolder "${PACKDIR}" || return 1;

    if [ -n "${SO_MFLIKE_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${SO_MFLIKE_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

    # PATCH MKFLIKE
    cdfolder "${PACKDIR:?}/mflike"|| return 1;

    cp "${CCCOB:?}/${COBLIKE:?}/mflike/mflike.patch" "${PACKDIR:?}/mflike" \
      2>"/dev/null" || 
      { error "CP FILE ${CCCOB:?}/${COBLIKE:?}/mflike/mflike.patch on ${PACKDIR:?}/mflike"; return 1; }

    patch -u "mflike.py" -i "mflike.patch" >${OUT1:?} \
        2>${OUT2:?} || { error "${EC17:?} (mflike.patch)"; return 1; }

    # --------------------------------------------------------------------------
    # need to create symlinks for this likelihood to work w/ COBAYA.
    # we copied the code below to start_cocoa.sh shell script as well
    # we added corresponding code to stop_cocoa.sh that delete these symlinks
    # --------------------------------------------------------------------------
    if [[ -d "${PACKDIR:?}/mflike" ]]; then
      if [[ -L "${COB:?}/${COBLIKE:?}/mflike" ]]; then
        rm -f "${COB:?}/${COBLIKE:?}/mflike"
      fi
      ln -s "${PACKDIR:?}/mflike" "${COB:?}/${COBLIKE:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }

      if [[ -L "${COB:?}/${COBTH:?}/mflike" ]]; then
        rm -f "${COB:?}/${COBTH:?}/mflike"
      fi
      ln -s "${PACKDIR:?}/mflike" "${COB:?}/${COBTH:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
    fi
  fi
  
  cdfolder "${ROOTDIR}" || return 1

  pbottom "SETUP SIMONS OBSERVATORY MFLIKE" || return 1;

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  cdfolder "${ROOTDIR}" || return 1

  unset_all || return 1;
  
fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------