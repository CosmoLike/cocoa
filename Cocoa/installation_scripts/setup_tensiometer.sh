#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_TENSIOMETER_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR PIPCP PIPCP_HASH SENTINEL_PIPCP 
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

URL="${TENSIOMETER_URL:-"https://github.com/mraveri/tensiometer.git"}"

FOLDER="${TENSIOMETER_NAME:-"tensiometer"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

ptop "INSTALLING TENSIOMETER" || { unset_all; return 1; }

if [ -n "${OVERWRITE_EXISTING_TENSIOMETER_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}"
  rm -f  "${ROOTDIR:?}/.local/pip_tensiometer_"*

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  PIPCP=(
    'autograd==1.8.0'
    'pymanopt==0.2.5'
  )

  PIPCP_HASH=$(
    {
      printf '%s ' "${PIPCP[@]}";
      printf 'PYTHON=%s\n' "${PYTHON3:-}";
      printf 'PIP=%s\n'    "${PIP3:-}";
      printf 'MPICC=%s\n'  "${MPI_CC_COMPILER:-}";
    } | md5sum | cut -d' ' -f1
  )

  SENTINEL_PIPCP="${ROOTDIR:?}/.local/pip_tensiometer_${PIPCP_HASH:?}"

  if [ ! -f "${SENTINEL_PIPCP:?}" ]; then
    
    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install "${PIPCP[@]}" \
      --no-cache-dir --prefer-binary --use-pep517 \
      --prefix="${ROOTDIR:?}/.local" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC13:?}"; return 1; }

    touch "${SENTINEL_PIPCP:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "SENTINEL PIP TEN"; return 1; }

  fi
  
  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

  if [[ -n "${TENSIOMETER_GIT_COMMIT:-}" ||
        -n "${TENSIOMETER_GIT_BRANCH:-}" ||
        -n "${TENSIOMETER_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${TENSIOMETER_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${TENSIOMETER_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${TENSIOMETER_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${TENSIOMETER_GIT_BRANCH:?}" "origin/${TENSIOMETER_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${TENSIOMETER_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${TENSIOMETER_GIT_TAG:?}" -b "${TENSIOMETER_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
  
fi

cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }

pbottom "INSTALLING TENSIOMETER" || { unset_all; return 1; }

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