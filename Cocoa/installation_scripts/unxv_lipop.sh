#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_LIPOP_CMB_DATA}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () { 
    unset -v EDATAF URL LPDVS FOLDER PRINTNAME TMP
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

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data"

  URL="${LIPOP_DATA_URL:-"https://portal.nersc.gov/cfs/cmb/planck2020/likelihoods"}"

  LPDVS=${LIPOP_DATA_VERSION:-"4.2"}

  # T = TMP
  declare -a FILE=( "planck_2020_hillipop_EE_v${LPDVS:?}" 
                    "planck_2020_hillipop_TE_v${LPDVS:?}" 
                    "planck_2020_hillipop_TT_v${LPDVS:?}"
                    "planck_2020_hillipop_TTTEEE_v${LPDVS:?}"
                    "planck_2020_lollipop")

  ptop "SETUP/UNXV LIPOP DATA" || return 1
    
  if [ -n "${OVERWRITE_EXISTING_LIPOP_CMB_DATA}" ]; then
    rm -rf "${EDATAF:?}/planck/planck_2020"
    rm -rf "${EDATAF:?}/planck/hillipop"  
    rm -rf "${EDATAF:?}/planck/lollipop"
    if [ -n "${REDOWNLOAD_EXISTING_LIPOP_CMB_DATA}" ]; then
      for (( i=0; i<${#FILE[@]}; i++ ));
      do
        rm -f  "${EDATAF:?}/planck/${FILE[$i]:?}.tar.gz"
      done 
    fi
  fi 
  
  if [[ ! -d "${EDATAF:?}/planck/hillipop" ||  
        ! -d "${EDATAF:?}/planck/lollipop" ]]; then
    cdfolder "${EDATAF:?}/planck"

    for (( i=0; i<${#FILE[@]}; i++ ));
    do
      if [ ! -e "${EDATAF:?}/planck/${FILE[$i]:?}.tar.gz" ]; then
        "${WGET:?}" "${URL}/${FILE[$i]:?}.tar.gz" \
          -q \
          --show-progress \
          --progress=bar:force:noscroll || { error "${EC24:?}"; return 1; }
      fi

      tar -xvzf "${FILE[$i]:?}.tar.gz" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?} (${FILE[$i]:?})"; return 1; }
    done

    # note: by default, LIPOP saves the likelihoods on planck_2020
    mv "${EDATAF:?}/planck/planck_2020/hillipop" "${EDATAF:?}/planck"
    mv "${EDATAF:?}/planck/planck_2020/lollipop" "${EDATAF:?}/planck"
    rm -rf "${EDATAF:?}/planck/planck_2020"
  fi
    
  pbottom "SETUP/UNXV LIPOP DATA" || return 1

  unset_all || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------