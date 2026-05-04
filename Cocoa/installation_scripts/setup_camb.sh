#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_CAMB_CODE:-}" ]; then
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

ptop "SETUP CAMB" || { unset_all; return 1; }

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

URL="${CAMB_URL:-"https://github.com/cmbant/CAMB"}"

CHANGES="${CCIL:?}/camb_changes"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

FOLDER="${CAMB_NAME:-"CAMB"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

if [ -n "${OVERWRITE_EXISTING_CAMB_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}"
<<<<<<< HEAD
  rm -rf "${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/camb*"
=======
>>>>>>> main

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { unset_all; return 1; }

  if [[ -n "${CAMB_GIT_COMMIT:-}" ||
        -n "${CAMB_GIT_BRANCH:-}" ||
        -n "${CAMB_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${CAMB_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${CAMB_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${CAMB_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${CAMB_GIT_BRANCH:?}" "origin/${CAMB_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${CAMB_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${CAMB_GIT_TAG:?}" -b "${CAMB_GIT_TAG:?}" \
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

<<<<<<< HEAD
    patch -l -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
=======
    patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
>>>>>>> main
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done

fi

pbottom "SETUP CAMB" || { unset_all; return 1; }

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