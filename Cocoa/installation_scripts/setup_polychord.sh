#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${IGNORE_POLYCHORD_SAMPLER_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER TFILE TFILEP AL PRINTNAME
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

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# ----------------------------------------------------------------------------

URL="${POLY_URL:-"https://github.com/PolyChord/PolyChordLite.git"}"
  
CHANGES="${CCIL:?}/polychord_changes"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

FOLDER=${POLY_NAME:-"PolyChordLite"}

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# Name to be printed on this shell script messages
PRINTNAME="POLYCHORD CODE"

ptop "INSTALLING ${PRINTNAME:?}" || { unset_all; return 1; }

# ---------------------------------------------------------------------------
# in case this script is called twice
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_POLYCHORD_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}"

fi

# ---------------------------------------------------------------------------
# clone from original repo
# ---------------------------------------------------------------------------
if [ ! -d "${PACKDIR:?}" ]; then

  cdfolder "${ECODEF}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { unset_all; return 1; }

  if [[ -n "${POLYCHORD_GIT_COMMIT:-}" ||
        -n "${POLYCHORD_GIT_BRANCH:-}" ||
        -n "${POLYCHORD_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${POLYCHORD_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${POLYCHORD_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${POLYCHORD_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${POLYCHORD_GIT_BRANCH:?}" "origin/${POLYCHORD_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${POLYCHORD_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${POLYCHORD_GIT_TAG:?}" -b "${POLYCHORD_GIT_TAG:?}TMP" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
  
  # ---------------------------------------------------------------------------
  # We patch the files below so they use the right compilers ------------------
  # ---------------------------------------------------------------------------
  declare -a TFOLDER=("" 
                      ""
                      "" ) # If nonblank, path must include /
  
  # T = TMP
  declare -a TFILE=("Makefile" 
                    "setup.py"
                    "Makefile_gnu")

  #T = TMP, P = PATCH
  declare -a TFILEP=("Makefile.patch" 
                     "setup.patch"
                     "Makefile_gnu.patch")

  # AL = Array Length
  AL=${#TFOLDER[@]}

  for (( i=0; i<${AL}; i++ ));
  do
    cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || { unset_all; return 1; }

    cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
      2>>${OUT2:?} || { unset_all; return 1; }

    patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done

fi

cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }

pbottom "INSTALLING ${PRINTNAME:?}" || { unset_all; return 1; }

#-------------------------------------------------------------------------------

unset_all || return 1;

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------