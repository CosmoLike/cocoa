#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_DERIVKIT_CODE:-}" ]; then
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

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

URL="${DERIVKIT_URL:-"https://github.com/nikosarcevic/derivkit.git"}"

FOLDER="${DERIVKIT_NAME:-"derivkit"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

ptop "INSTALLING DERIVKIT DERIVATIVE CALCULATOR" || { unset_all; return 1; }

if [ -n "${OVERWRITE_EXISTING_DERIVKIT_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}"

fi

if [ ! -d "${PACKDIR:?}" ]; then
  env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
    'numdifftools==0.9.41' \
    --no-cache-dir --prefer-binary --use-pep517 --prefix="${ROOTDIR:?}/.local" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

  if [[ -n "${DERIVKIT_GIT_COMMIT:-}" ||
        -n "${DERIVKIT_GIT_BRANCH:-}" ||
        -n "${DERIVKIT_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${DERIVKIT_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${DERIVKIT_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${DERIVKIT_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${DERIVKIT_GIT_BRANCH:?}" "origin/${DERIVKIT_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${DERIVKIT_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${DERIVKIT_GIT_TAG:?}" -b "${DERIVKIT_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
 
fi

cdfolder "${ROOTDIR}" || { unset_all; return 1; }

pbottom "INSTALLING DERIVKIT DERIVATIVE CALCULATOR" || { unset_all; return 1; }

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