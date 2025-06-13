#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

unset_env_vars () {
  unset -v PRINTNAME FOLDER URL PACKAGE_VERSION PACKAGE_BRANCH
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error gitact gitact1 gitact2
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

gitact() {  
  local PROJECT="${ROOTDIR}/projects"

  local PACKDIR="${PROJECT:?}/${1:?}"

  cdfolder "${PROJECT:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_PRIVATE_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  # ---------------------------------------------------------------------------
  # clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  if [ ! -d "${PACKDIR:?}" ]; then
    "${GIT:?}" clone "${2:?}" "${1:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

    cdfolder "${1}" || return 1;

    "${GIT:?}" checkout ${3} \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR}" || return 1;
}

gitact2() {  
  local PROJECT="${ROOTDIR}/projects"

  local PACKDIR="${PROJECT:?}/${1:?}"

  cdfolder "${PROJECT:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_PRIVATE_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  # ---------------------------------------------------------------------------
  # clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  if [ ! -d "${PACKDIR:?}" ]; then
    "${GIT:?}" clone "${2:?}" "${1:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR}" || return 1;
}


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------
# -------------------------------- AXIONS -----------------------------------
# ----------------------------------------------------------------------------

if [ -n "${INSTALL_PRIVATE_AXIONS_PROJECT}" ]; then 
  
  echo "what?"

  # Name to be printed on this shell script messages
  PRINTNAME="AXIONS PROJECT"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${AXIONS_PROJECT_NAME:-"axions"}"

  URL="${AXIONS_PROJECT_URL:?}"

  if [ -n "${AXIONS_PROJECT_COMMIT}" ]; then
    gitact "${FOLDER:?}" "${URL:?}" "${AXIONS_PROJECT_COMMIT:?}"  || return 1
  else
    gitact2 "${FOLDER:?}" "${URL:?}"  || return 1
  fi

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------