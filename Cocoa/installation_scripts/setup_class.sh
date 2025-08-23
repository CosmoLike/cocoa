#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi
  
  # Parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v URL CHANGES ECODEF FOLDER PACKDIR TFOLDER TFILE TFILEP AL PRINTNAME
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
    cp -r "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
    
  unset_env_vars || return 1

  # ---------------------------------------------------------------------------

  URL="${CLASS_URL:-"https://github.com/lesgourg/class_public.git"}"
  
  CHANGES="${ROOTDIR:?}/../cocoa_installation_libraries/class_changes"

  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER=${CLASS_NAME:-"class_public"}
  
  PACKDIR="${ECODEF:?}/${FOLDER:?}/"

  # Name to be printed on this shell script messages
  PRINTNAME="CLASS"

  ptop "SETUP ${PRINTNAME:?}" || return 1;

  # ---------------------------------------------------------------------------
  # in case this script is called twice ---------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_CLASS_CODE}" ]; then
    
    rm -rf "${PACKDIR:?}"
  
  fi

  # ---------------------------------------------------------------------------
  # clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  if [ ! -d "${PACKDIR:?}" ]; then

    cdfolder "${ECODEF}" || return 1;
    
    "${GIT:?}" clone --depth ${GIT_CLONE_MAXIMUM_DEPTH:?} "${URL:?}" \
      --recursive "${FOLDER:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    cdfolder "${PACKDIR}" || return 1;

    if [ -n "${CLASS_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${CLASS_GIT_COMMIT:?}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16\?}"; return 1; }
    fi
    
    # --------------------------------------------------------------------------
    # historical reasons (we used to save class_python on Cocoa Branch) --------
    # historical: Workaround Cocoa .gitignore entry on /include ----------------
    # --------------------------------------------------------------------------
    mv ./include ./include2/ \
      >>${OUT1:?} 2>>${OUT2:?} || 
      { error "MKDIR FROM INCLUDE TO INCLUDE2"; return 1; }
    
    # --------------------------------------------------------------------------
    # We patch the files below so they use the right C compiler ----------------
    # --------------------------------------------------------------------------
    # T = TMP
    declare -a TFOLDER=("" 
                        "python/") # Must include
    
    # T = TMP
    declare -a TFILE=("Makefile" 
                      "setup.py")

    #T = TMP, P = PATCH
    declare -a TFILEP=("Makefile.patch" 
                       "setup.patch")

    # AL = Array Length
    AL=${#TFOLDER[@]}

    for (( i=0; i<${AL}; i++ ));
    do
      cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || return 1;

      cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
        2>>${OUT2:?} || return 1;

      patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" >>${OUT1:?} \
        2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
    done

  fi
  
  cdfolder "${ROOTDIR}" || return 1;

  pbottom "SETUP ${PRINTNAME:?}" || return 1

  unset_all || return 1
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------