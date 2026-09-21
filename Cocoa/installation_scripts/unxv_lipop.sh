#!/bin/bash
# ------------------------------------------------------------------------------
# unxv_lipop.sh: unpack the Planck LoLLiPoP/HiLLiPoP likelihoods into
# external_modules/data.
#
# Sourced by setup_cocoa.sh, never executed directly. Datasets are downloaded
# or unpacked at setup time only, never at compile or run time.
#
# Skipped (return 99) when IGNORE_LIPOP_CMB_DATA is set. Reruns are safe:
# existing results are kept unless OVERWRITE_EXISTING_LIPOP_CMB_DATA is set.
# On failure error() names the failing step and unsets everything defined
# here; on success the script returns 55, the value the runner caches so a
# finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_LIPOP_CMB_DATA:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell 
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
  
# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () { 
  unset -v EDATAF URL LPDVS FOLDER PRINTNAME TMP FILE
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder error
  unset -f unset_env_funcs
  cdroot || return 1;
}

# unset_all: the single cleanup entry point: variables, then functions, then
# itself.
unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

# error: print which step failed (fail_script_msg names this script and the
# message) and run the FULL cleanup. A call site whose helper already routes
# through error() therefore uses plain `|| return 1`; adding another unset_all
# there would run the cleanup twice.
error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

# cdfolder: cd that reports the target and cleans up on failure.
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

ptop "SETUP/UNXV LIPOP DATA" || { unset_all; return 1; }
  
if [ -n "${OVERWRITE_EXISTING_LIPOP_CMB_DATA:-}" ]; then
  rm -rf "${EDATAF:?}/planck/planck_2020"
  rm -rf "${EDATAF:?}/planck/hillipop"  
  rm -rf "${EDATAF:?}/planck/lollipop"
  
  if [ -n "${REDOWNLOAD_EXISTING_LIPOP_CMB_DATA:-}" ]; then
    
    for (( i=0; i<${#FILE[@]}; i++ ));
    do
      rm -f  "${EDATAF:?}/planck/${FILE[$i]:?}.tar.gz"
    done 
  
  fi

fi 

if [[ ! -d "${EDATAF:?}/planck/hillipop" ||  
      ! -d "${EDATAF:?}/planck/lollipop" ]]; then
  
  cdfolder "${EDATAF:?}/planck" || return 1;

  for (( i=0; i<${#FILE[@]}; i++ ));
  do
    if [ ! -e "${EDATAF:?}/planck/${FILE[$i]:?}.tar.gz" ]; then
    
      "${WGET:?}" "${URL}/${FILE[$i]:?}.tar.gz" -q --show-progress \
        --no-check-certificate --progress=bar:force:noscroll \
        --timeout=30 --tries=2 --waitretry=0 --retry-connrefused \
        --read-timeout=30 || { error "${EC24:?}"; return 1; }
    
    fi

    tar -xvzf "${FILE[$i]:?}.tar.gz" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?} (${FILE[$i]:?})"; return 1; }
  done

  # note: by default, LIPOP saves the likelihoods on planck_2020
  mv "${EDATAF:?}/planck/planck_2020/hillipop" "${EDATAF:?}/planck" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "MV HILLIPOP DATA"; return 1; }
  
  mv "${EDATAF:?}/planck/planck_2020/lollipop" "${EDATAF:?}/planck" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "MV LOLLIPOP DATA"; return 1; }
  
  rm -rf "${EDATAF:?}/planck/planck_2020"
fi
  
pbottom "SETUP/UNXV LIPOP DATA" || { unset_all; return 1; }

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
