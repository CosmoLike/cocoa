#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_COSMOREC_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER 
    unset -v TFILE TFILEP AL PRINTNAME URL_BASE
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

  FILENAME="${COSMOREC_CODE_FILE:-"CosmoRec.v2.0.3b"}"

  EXT="${COSMOREC_CODE_FILE_EXT:-"tar.gz"}"

  FILE="${FILENAME:?}.${EXT:?}"

  URL_BASE="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"

  URL="${COSMOREC_URL:-"${URL_BASE:?}"}/${FILE:?}"

  CHANGES="${CCIL:?}/cosmorec_changes"

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${COSMOREC_NAME:-"cosmorec"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="COSMOREC RECOMBINATION CODE"

  ptop "INSTALLING ${PRINTNAME:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script is called twice ---------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOREC_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
    rm -f  "${ECODEF:?}/${FILE}"
    rm -rf "${ECODEF:?}/${FILENAME}"
  fi
  # ----------------------------------------------------------------------------
  # Clone from original repo ---------------------------------------------------
  # ----------------------------------------------------------------------------
  if [ ! -d "${PACKDIR:?}" ]; then
  
    cdfolder "${ECODEF:?}" || { cdroot; return 1; }

    "${WGET:?}" "${URL:?}" -q --show-progress --progress=bar:force || { error "${EC24:?}"; return 1; }
    
    tar -zxvf "${FILE:?}" >>${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

    mv "${FILENAME:?}" "${FOLDER:?}" >>${OUT1:?} 2>>${OUT2:?} || { error "${EC30:?}"; return 1; }

    rm -f  "${ECODEF:?}/${FILE}"
    rm -rf "${ECODEF:?}/${FILENAME}"

    # --------------------------------------------------------------------------
    # We patch the files below so they use the right compilers -----------------
    # --------------------------------------------------------------------------
    # PREFIX: T = TMP, P = PATCH, AL = Array Length
    declare -a TFOLDER=("") # If nonblank, path must include /
    
    declare -a TFILE=("Makefile.in")

    case "$(uname -s)" in
      Linux)
        declare -a TFILEP=("Makefile.in.patch")
        ;;
      Darwin)
        declare -a TFILEP=("MakefileOSX.in.patch")
        ;;
    esac

    AL=${#TFOLDER[@]}

    for (( i=0; i<${AL}; i++ ));
    do
      cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || return 1

      cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
        2>${OUT2:?} || return 1;

      patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" >${OUT1:?} \
        2>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
    done

  fi
  
  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "INSTALLING ${PRINTNAME:?}" || return 1
  
  unset_all || return 1
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------