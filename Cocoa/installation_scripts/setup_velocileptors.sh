#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_VELOCILEPTORS_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR
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

ptop 'INSTALLING VELOCILEPTORS' || { unset_all; return 1; }

URL="${VELOCILEPTORS_URL:-"https://github.com/sfschen/velocileptors.git"}"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

FOLDER=${VELOCILEPTORS_NAME:-"velocileptors"}

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# ---------------------------------------------------------------------------
# In case this script is called twice ---------------------------------------
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_VELOCILEPTORS_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}"

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  # ---------------------------------------------------------------------------
  # Clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { unset_all; return 1; }

  if [[ -n "${VELOCILEPTORS_GIT_COMMIT:-}" ||
        -n "${VELOCILEPTORS_GIT_BRANCH:-}" ||
        -n "${VELOCILEPTORS_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${VELOCILEPTORS_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${VELOCILEPTORS_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${VELOCILEPTORS_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout "${VELOCILEPTORS_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${VELOCILEPTORS_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${VELOCILEPTORS_GIT_TAG:?}" -b "${VELOCILEPTORS_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

fi
  
cdfolder "${ROOTDIR}" || { unset_all; return 1; }

pbottom 'INSTALLING VELOCILEPTORS' || { unset_all; return 1; }

#-------------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------