#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${INSTALL_AXIE_CAMB_V2:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER 
  unset -v TFILE TFILEP AL PRINTNAME
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

# ---------------------------------------------------------------------------

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "SETUP AXIECAMB" || { unset_all; return 1; }

URL="${AXIE_CAMB_URL:-"https://github.com/SBU-COSMOLIKE/NewTestingAxieCAMB.git"}"

CHANGES="${CCIL:?}/axiecamb_changes"

FOLDER="${AXIE_CAMB_NAME:-"AXIECAMB"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

if [ -n "${OVERWRITE_EXISTING_AXIECAMB_CODE:-}" ]; then
  rm -rf "${PACKDIR:?}"
fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { unset_all; return 1; }

  if [[ -n "${AXIE_CAMB_GIT_COMMIT:-}" ||
        -n "${AXIE_CAMB_GIT_BRANCH:-}" ||
        -n "${AXIE_CAMB_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${AXIE_CAMB_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${AXIE_CAMB_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${AXIE_CAMB_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${AXIE_CAMB_GIT_BRANCH:?}" "origin/${AXIE_CAMB_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${AXIE_CAMB_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${AXIE_CAMB_GIT_TAG:?}" -b "${AXIE_CAMB_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
  
  # --------------------------------------------------------------------------
  # We patch the files below so they use the right compilers -----------------
  # --------------------------------------------------------------------------
  # PREFIX: T = TMP, P = PATCH, AL = Array Length
  declare -a TFOLDER=("camb/" 
                      "fortran/" 
                      "forutils/"
                      "fortran/") # If nonblank, path must include /
  declare -a TFILE=("_compilers.py" 
                    "Makefile" 
                    "Makefile_compiler"
                    "Makefile_main")
  declare -a TFILEP=("_compilers.patch" 
                     "Makefile.patch" 
                     "Makefile_compiler.patch"
                     "Makefile_main.patch")
  case "$(uname -s)" in
    Linux)
      declare -a TFILEP=("_compilers.patch" 
                         "Makefile.patch" 
                         "Makefile_compiler.patch"
                         "Makefile_main.patch")
      ;;
    Darwin)
      declare -a TFILEP=("_compilers.patch" 
                         "Makefile.patch" 
                         "Makefile_compiler.patch"
                         "Makefile_main_osx.patch")
      ;;
  esac
  
  AL=${#TFOLDER[@]}

  for (( i=0; i<${AL}; i++ ));
  do
    cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || { unset_all; return 1; }

    cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
      2>>${OUT2:?} || { unset_all; return 1; }

    patch -l -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done

fi

pbottom "SETUP AXIECAMB" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "SETUP AXION HMCODE" || { unset_all; return 1; }

URL="${AXION_HMCODE_URL:-"https://github.com/SophieMLV/axionHMcode.git"}"

FOLDER="${AXION_HMCODE_NAME:-"axionHMcode"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

if [ -n "${OVERWRITE_EXISTING_AXION_HMCODE_CODE:-}" ]; then
  rm -rf "${PACKDIR:?}"
fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { unset_all; return 1; }

  if [[ -n "${AXION_HMCODE_GIT_COMMIT:-}" ||
        -n "${AXION_HMCODE_GIT_BRANCH:-}" ||
        -n "${AXION_HMCODE_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${AXION_HMCODE_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${AXION_HMCODE_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${AXION_HMCODE_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${AXION_HMCODE_GIT_BRANCH:?}" "origin/${AXION_HMCODE_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${AXION_HMCODE_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${AXION_HMCODE_GIT_TAG:?}" -b "${AXION_HMCODE_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

fi

pbottom "SETUP AXION HMCODE" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

cdfolder "${ROOTDIR}" || { unset_all; return 1; }

unset_all || return 1

#-----------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------