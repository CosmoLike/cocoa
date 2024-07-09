#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_XXX_DATA}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () { 
    unset -v EDATAF URL PACKDIR FOLDER FILE EXT PRINTNAME
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

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data"

  URL="${XXX_DATA_URL:-"https://website/XXX"}"

  # FOLDER = the directory name of the dataset
  FOLDER="XXX"

  declare -a FILE=( "filename1" 
                    "filename2" 
                    "filename3"
                  )

  declare -a EXT=( "tar.gz" 
                   "tar.gz" 
                   "tar.gz"
                  )

  # Name to be printed on this shell script messages
  PRINTNAME="XXX"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  # ---------------------------------------------------------------------------

  ptop "SETUP/UNXV ${PRINTNAME:?} DATA" || return 1
  
  # ----------------------------------------------------------------------------
  # note: in case script run >1x

  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  
  mkdir ${PACKDIR:?}

  cdfolder "${PACKDIR:?}"

  for (( i=0; i<${#FILE[@]}; i++ ));
  do

    "${WGET:?}" "${URL}/${FILE[$i]:?}.${EXT[$i]:?}" -q --show-progress \
      --progress=bar:force:noscroll || { error "${EC24:?}"; return 1; }

    if [ "${EXT:?}" == "tar.gz" ]; then
      
      tar -xvzf "${FILE[$i]:?}.${EXT[$i]:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
    
    elif [ "${EXT:?}" == "tar.xz" ]; then
    
      tar -xf "${FILE[$i]:?}.${EXT[$i]:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
    
    else
      
      error "${EC29:?}"; return 1;
  
    fi

  done

  # ----------------------------------------------------------------------------
  # erase compressed files 
  
  for (( i=0; i<${#FILE[@]}; i++ ));
  do
    rm -f "${PACKDIR:?}/${FILE[$i]:?}.${EXT[$i]:?}"
  done 
  
  # ---------------------------------------------------------------------------
    
  pbottom "SETUP/UNXV ${PRINTNAME:?} DATA" || return 1

  unset_all || return 1;
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------