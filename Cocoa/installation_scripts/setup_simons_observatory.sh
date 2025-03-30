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
  unset -v IGNORE_FGSPECTRA_CODE
  ( source "${ROOTDIR:?}/installation_scripts/setup_fgspectra.sh" ) || return 1;

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

  COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

  CCCOB="${CCIL:?}/cobaya_changes"  # CC = CoCoA, COB = Cobaya (Cocoa Cobaya)
  
  COBLIKE="cobaya/likelihoods"      # COB = Cobaya, LIKE = likelihoods
  
  cppatch() {
    cp "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" \
      2>"/dev/null" || 
      { error "CP FILE ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  URL="${SO_SYSLIB_URL:-"https://github.com/simonsobs/syslibrary.git"}"
    
  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${SO_SYSLIB_NAME:-"SOSYSLIB"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="SIMONS OBSERVATORY SYSLIBRARY"

  ptop "SETUP ${PRINTNAME:?}" || return 1;

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE}" ]; then

    rm -rf "${PACKDIR:?}"

  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    
    # ---------------------------------------------------------------------------
    # clone from original repo
    # ---------------------------------------------------------------------------
    cdfolder "${ECODEF}" || return 1;

    "${CURL:?}" -fsS "${URL:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

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

  pbottom "SETUP ${PRINTNAME:?}" || return 1;
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  URL="${SO_MFLIKE_URL:-"https://github.com/simonsobs/LAT_MFLike.git"}"
  
  FOLDER="SOMKFLIKE"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="SIMONS OBSERVATORY MFLIKE"

  ptop "SETUP ${PRINTNAME:?}" || return 1;

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"
  rm -rf "${COB:?}/${COBLIKE}/mflike"

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF}" || return 1;

  "${CURL:?}" -fsS "${URL:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

  "${GIT:?}" clone "${URL:?}" --recursive "${FOLDER:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR}" || return 1;

  if [ -n "${SO_MFLIKE_GIT_COMMIT}" ]; then
    "${GIT:?}" checkout "${SO_MFLIKE_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

  # ---------------------------------------------------------------------------
  # move likelihood to cobaya
  
  mv "${PACKDIR}/mflike" "${COB:?}/${COBLIKE}"

  # ---------------------------------------------------------------------------
  # ERASE TMP FOLDER
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # PATCH MKFLIKE
  
  cppatch "${COBLIKE:?}/mflike" "mflike.patch" || return 1

  cdfolder "${COB:?}/${COBLIKE:?}/mflike"|| return 1;

  patch -u "mflike.py" -i "mflike.patch" >${OUT1:?} \
      2>${OUT2:?} || { error "${EC17:?} (mflike.patch)"; return 1; }

  cdfolder "${ROOTDIR}" || return 1

  pbottom "SETUP ${PRINTNAME:?}" || return 1;

  #-----------------------------------------------------------------------------
  
  unset_all || return 1;
  
fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------