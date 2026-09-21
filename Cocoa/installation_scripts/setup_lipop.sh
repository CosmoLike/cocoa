#!/bin/bash
# ------------------------------------------------------------------------------
# setup_lipop.sh: download the Planck LoLLiPoP/HiLLiPoP likelihoods.
#
# Sourced by setup_cocoa.sh, never executed directly. Setup scripts MAY use
# the internet (git clone, wget, pip index); the matching compile step never
# does, so anything that must be downloaded happens here.
#
# Skipped (return 99) when IGNORE_LIPOP_LIKELIHOOD_CODE is set. Reruns are
# safe: existing results are kept unless OVERWRITE_EXISTING_LIPOP_CMB_CODE is
# set. On failure error() names the failing step and unsets everything defined
# here; on success the script returns 55, the value the runner caches so a
# finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_LIPOP_LIKELIHOOD_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v COB CCCOB COBLIKE URL TFILE TFOLDER ECODEF CCIL
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder error fhilipop flolipop
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

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods" # COB = Cobaya, LIKE = likelihoods

fhilipop() {
  local TF="${ECODEF:?}/${1:?}"
  
  if [ -n "${OVERWRITE_EXISTING_LIPOP_CMB_CODE:-}" ]; then
  
    rm -rf "${TF:?}"
  
  fi

  if [ ! -d "${TF:?}" ]; then
    cdfolder "${ECODEF:?}" || return 1;

    "${GIT:?}" clone "${3:?}" --recursive "${1:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
    cdfolder "${TF:?}" || return 1;

    "${GIT:?}" reset --hard "${4:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC23:?}"; return 1; }

    cp "${CCIL:?}/${2:?}_changes/init.patch" "${TF:?}/${2:?}" \
      2>>"/dev/null" || { error "CP FILE init.patch"; return 1; }

    cdfolder "${TF:?}/${2:?}" || return 1;

    patch -u '__init__.py' -i 'init.patch' \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?}"; return 1; }
  fi

  cdfolder "${ROOTDIR:?}" || return 1;
}

flolipop() {
  local TF="${ECODEF:?}/${1:?}"
  
  if [ -n "${OVERWRITE_EXISTING_LIPOP_CMB_CODE}" ]; then
  
    rm -rf "${TF:?}"
  
  fi

  if [ ! -d "${TF:?}" ]; then
    
    cdfolder "${ECODEF:?}" || return 1;

    "${GIT:?}" clone "${3:?}" --recursive  "${1:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
    cdfolder "${TF:?}" || return 1;

    "${GIT:?}" reset --hard "${4:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC23:?}"; return 1; }

    cp "${CCIL:?}/${2:?}_changes/init.patch" "${TF:?}/${2:?}" \
      2>>"/dev/null" || { error "CP FILE init.patch"; return 1; }

    cp "${CCIL:?}/${2:?}_changes/lollipop.patch" "${TF:?}/${2:?}" \
      2>>"/dev/null" || { error "CP FILE lollipop.patch"; return 1; }

    cdfolder "${TF:?}/${2:?}" || return 1;

    patch -u '__init__.py' -i 'init.patch' \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?}"; return 1; }

    patch -u 'lollipop.py' -i 'lollipop.patch' \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?}"; return 1; }
  fi

  cdfolder "${ROOTDIR:?}" || return 1;
}

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

ptop "SETUP HILLIPOP LIKELIHOOD" || { unset_all; return 1; }

TFOLDER="${PL2020_HILLIPOP_NAME:-"planck_2020_hillipop"}"

URL="${HILLIPOP_URL:-"https://github.com/planck-npipe/hillipop.git"}"

fhilipop "${TFOLDER:?}" "planck_2020_hillipop" "${URL:?}" \
  "${HILLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;

pbottom "SETUP HILLIPOP LIKELIHOOD" || { unset_all; return 1; }

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

ptop "SETUP LOLLIPOP LIKELIHOOD" || { unset_all; return 1; }

TFOLDER="${PL2020_LOLLIPOP_NAME:-"planck_2020_lollipop"}"

URL="${LOLLIPOP_URL:-"https://github.com/planck-npipe/lollipop.git"}"

flolipop "${TFOLDER:?}" "planck_2020_lollipop" "${URL:?}" \
  "${LOLLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;

pbottom "SETUP LOLLIPOP LIKELIHOOD" || { unset_all; return 1; }

#-----------------------------------------------------------------------------

cdfolder "${ROOTDIR:?}" || return 1;

#-------------------------------------------------------------------------------

unset_all || return 1;

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
