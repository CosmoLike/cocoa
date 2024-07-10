#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_HYREC_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
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

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ---------------------------------------------------------------------------

  URL="${HYREC_URL:-"https://github.com/nanoomlee/HYREC-2.git"}"

  CHANGES="${CCIL:?}/hyrec_changes"

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${HYREC_NAME:-"hyrec2"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="HYREC2 RECOMBINATION CODE"

  ptop "INSTALLING ${PRINTNAME:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script is called twice ---------------------------------------
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # Clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || { cdroot; return 1; }

  "${CURL:?}" -fsS "${URL:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

  "${GIT:?}" clone "${URL:?}" --recursive "${FOLDER:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { cdroot; return 1; }

  if [ -n "${HYREC_GIT_COMMIT}" ]; then
    "${GIT:?}" checkout "${HYREC_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

  # ---------------------------------------------------------------------------
  # We patch the files below so they use the right compilers ------------------
  # ---------------------------------------------------------------------------
  # T = TMP
  declare -a TFOLDER=("" 
                     ) # If nonblank, path must include /
  
  # T = TMP
  declare -a TFILE=("Makefile")

  #T = TMP, P = PATCH
  declare -a TFILEP=("Makefile.patch" 
                    )

  # AL = Array Length
  AL=${#TFOLDER[@]}

  for (( i=0; i<${AL}; i++ ));
  do
    cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || return 1

    cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
      2>${OUT2:?} || return 1;

    patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" >${OUT1:?} \
      2>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done
  
  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "INSTALLING ${PRINTNAME:?}" || return 1
  
  # ---------------------------------------------------------------------------

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------