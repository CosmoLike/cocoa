#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_MGCAMB_CODE}" ]; then

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

  URL="${MGCAMB_URL:-"https://github.com/sfu-cosmo/MGCobaya.git"}"

  CHANGES="${CCIL:?}/camb_changes"

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${MGCAMB_NAME:-"MGCAMB"}"

  # T = TMP
  TFOLDER="TMPMGCAMB"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="MGCAMB"

  ptop "INSTALLING ${PRINTNAME:?}" || return 1;

  # ----------------------------------------------------------------------------
  # In case this script is called twice ----------------------------------------
  # ----------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_MGCAMB_CODE}" ]; then
    rm -rf "${ECODEF:?}/${TFOLDER:?}"
    rm -rf "${PACKDIR:?}"
  fi

  # ----------------------------------------------------------------------------
  # Clone from original repo ---------------------------------------------------
  # ----------------------------------------------------------------------------
  if [ ! -d "${PACKDIR:?}" ]; then

    cdfolder "${ECODEF:?}" || { cdroot; return 1; }

    "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:?} \
      --recursive "${TFOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    cdfolder "${ECODEF:?}/${TFOLDER:?}" || { cdroot; return 1; }

    if [ -n "${MGCAMB_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${MGCAMB_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
    
    # ---------------------------------------------------------------------------
    # move MGCAMB folders
    # ---------------------------------------------------------------------------
    mv "${ECODEF:?}/${TFOLDER:?}/MGCAMB" "${ECODEF:?}"
    
    if [ "MGCAMB/" != "${FOLDER:?}" ] && [ "MGCAMB" != "${FOLDER:?}" ] ; then
      mv "${ECODEF:?}/MGCAMB" ${FOLDER:?}
    fi

    rm -rf "${ECODEF:?}/${TFOLDER:?}"

    cdfolder "${PACKDIR}" || { cdroot; return 1; }

    # ---------------------------------------------------------------------------
    # We patch the files below so they use the right compilers ------------------
    # ---------------------------------------------------------------------------
    # T = TMP
    declare -a TFOLDER=("camb/" 
                        "fortran/" 
                        "forutils/") # If nonblank, path must include /
    
    # T = TMP
    declare -a TFILE=("_compilers.py" 
                      "Makefile" 
                      "Makefile_compiler")

    #T = TMP, P = PATCH
    declare -a TFILEP=("_compilers.patch" 
                       "Makefile.patch" 
                       "Makefile_compiler.patch")

    # AL = Array Length
    AL=${#TFOLDER[@]}

    for (( i=0; i<${AL}; i++ ));
    do
      cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || return 1

      cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . 2>>${OUT2:?} || return 1;

      patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" >${OUT1:?} \
      2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
    done
  
  fi
  
  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "INSTALLING ${PRINTNAME:?}" || return 1
  
  unset_all || return 1
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------