#!/bin/bash
# ------------------------------------------------------------------------------
# setup_cobaya.sh: download Cobaya (the Bayesian analysis framework Cocoa
# drives).
#
# Sourced, never executed (only `source` keeps its environment changes).
# setup_cocoa.sh runs it during a full installation, and sourcing it on its
# own redoes just this package (see the README FAQ on installing a single
# package). Setup scripts MAY use the internet (git clone, wget, pip index);
# the matching compile step never does, so anything that must be downloaded
# happens here.
#
# Skipped (return 99) when IGNORE_COBAYA_CODE is set. Reruns are safe:
# existing results are kept unless OVERWRITE_EXISTING_COBAYA_CODE is set. On
# failure error() names the failing step and unsets everything defined here;
# on success the script returns 55, the value the runner caches so a finished
# step is not repeated.
# ------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -n "${IGNORE_COBAYA_CODE:-}" ]; then
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
  unset -v COB CCCOB COBLIKE COBTH URL TFILE TFOLDER
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder error cpfile cppatch cppatchfolder
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

# cpfolder: recursive copy that reports source/target and cleans up on
# failure.
cpfolder() {
  cp -r "${1:?}" "${2:?}"  \
    2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
}

# cpfile: file copy that reports source/target and cleans up on failure.
cpfile() {
  cp "${1:?}" "${2:?}" \
    2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset_env_vars || return 1

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

CCCOB="${CCIL:?}/cobaya_changes"  # CC = CoCoA, COB = Cobaya (Cocoa Cobaya)

COBLIKE="cobaya/likelihoods"      # COB = Cobaya, LIKE = likelihoods

COBTH="cobaya/theories"           # COB = Cobaya, TH = theory

cppatch() {
  cp "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" || 
    { error "CP FILE ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; return 1; }
}

cppatchfolder() {
  cp -r "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" || 
    { error "CP FOLDER ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; \
    return 1; }
}

# ----------------------------------------------------------------------------
  
ptop "SETUP COBAYA" || { unset_all; return 1; }

#-----------------------------------------------------------------------------
# Remove any previous installed cobaya folder --------------------------------
#-----------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_COBAYA_CODE:-}" ]; then

  rm -rf "${ROOTDIR:?}/cobaya"

fi

if [ ! -d "${ROOTDIR:?}/cobaya" ]; then

  #-----------------------------------------------------------------------------
  # Clone Cobaya from original git repo ----------------------------------------
  #-----------------------------------------------------------------------------
  cdfolder "${ROOTDIR:?}" || return 1;
  
  URL="${COBAYA_URL:-"https://github.com/CobayaSampler/cobaya.git"}"

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch cobaya \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  unset URL
  
  cdfolder "${COB:?}" || return 1;

  if [[ -n "${COBAYA_GIT_COMMIT:-}" ||
        -n "${COBAYA_GIT_BRANCH:-}" ||
        -n "${COBAYA_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${COBAYA_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${COBAYA_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${COBAYA_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${COBAYA_GIT_BRANCH:?}" "origin/${COBAYA_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${COBAYA_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${COBAYA_GIT_TAG:?}" -b "${COBAYA_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

  cdfolder "${ROOTDIR:?}" || return 1;

  #-----------------------------------------------------------------------------
  # Adjust Cobaya Files --------------------------------------------------------
  #-----------------------------------------------------------------------------
  
  TFILE="change_python_files.sh" # T = TMP
  
  cppatch "cobaya" "${TFILE:?}" || return 1;

  cdfolder "${COB:?}/cobaya/" || return 1;
  
  # parenthesis = run in a subshell
  ( sh change_python_files.sh ) || { error "${EC22:?} (CPF)"; return 1; }

  # --------------------------------------------------------------------------
  # PATCH FILE --------------------------------------------- -----------------
  # --------------------------------------------------------------------------  
  # instead of using JAX-CLIK - uses the original CLIK in FORTRAN
  declare -a TFOLDER=("cobaya/likelihoods/base_classes/") # Must include  
  declare -a TFILE=("planck_clik.py")
  declare -a TFILEP=("planck_clik.patch")

  for (( i=0; i<${#TFOLDER[@]}; i++ ));
  do
    cdfolder "${COB:?}/${TFOLDER[$i]}" || return 1;
  
    cpfolder "${CCCOB:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
      2>>${OUT2:?} || return 1;
  
   patch --quiet --batch --verbose -u -R "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done

  # --------------------------------------------------------------------------
  # PATCH FILE2 --------------------------------------------- -----------------
  # --------------------------------------------------------------------------  
  declare -a TFOLDER=("cobaya/") # Must include  
  declare -a TFILE=("model.py")
  declare -a TFILEP=("model.patch")

  for (( i=0; i<${#TFOLDER[@]}; i++ ));
  do
    cdfolder "${COB:?}/${TFOLDER[$i]}" || return 1;
  
    cpfolder "${CCCOB:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
      2>>${OUT2:?} || return 1;
  
    # HERE I CAN't USE THE -R
    patch --quiet --batch --verbose -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done
  
  cdfolder "${ROOTDIR:?}" || return 1;
  
  unset -v TFILE

  #-----------------------------------------------------------------------------
  # Adjust SN / STRONG LENSING -------------------------------------------------
  #-----------------------------------------------------------------------------
  
  TFOLDER="${COBLIKE:?}/sn" # T = TMP

  cp "${CCCOB:?}/${TFOLDER:?}/"roman_* "${COB:?}/${TFOLDER:?}/" \
    2>>${OUT2:?} || { error "CP ROMAN FILES"; return 1; }

  cppatchfolder "${COBLIKE:?}" "h0licow" || return 1;

  unset -v TFOLDER

  #-----------------------------------------------------------------------------
  # Remove native DES-Y1 cobaya likelihood -------------------------------------
  #-----------------------------------------------------------------------------
  
  rm -rf "${COB:?}/${COBLIKE:?}/des_y1"

  #-----------------------------------------------------------------------------
  # Remove Planck 2015 likelihood files ----------------------------------------
  #-----------------------------------------------------------------------------
  
  rm -rf "${COB:?}/${COBLIKE:?}"/planck_2015_*

  #-----------------------------------------------------------------------------
  # Adjust Planck 2018 CLIK ----------------------------------------------------
  #-----------------------------------------------------------------------------
  
  TFOLDER="${COBLIKE:?}/base_classes" # T = TMP
  
  cppatch "${TFOLDER:?}" "change_planck_clik.sh" || return 1

  cdfolder "${COB:?}/${TFOLDER}" || return 1;
  
  # parenthesis = run in a subshell
  ( sh change_planck_clik.sh ) || { error "${EC22:?} (CPC)"; return 1; }
  
  unset -v TFOLDER

  cdfolder "${ROOTDIR:?}" || return 1;

  #---------------------------------------------------------------------------
  # Adjust Planck 2018 low-ell -----------------------------------------------
  #---------------------------------------------------------------------------
  
  # note 1: we remove Cobaya native python CLIK implementation
  # note 2: In the original Cobaya, TT_clik.py/EE_clik.py are the files 
  # note 2: associated with CLIK. We rename them to simply TT/EE.py
  
  TFOLDER="${COBLIKE:?}/planck_2018_lowl"
  
  rm -f "${COB:?}/${TFOLDER:?}/EE_clik.py"
  rm -f "${COB:?}/${TFOLDER:?}/EE_clik.yaml"
  rm -f "${COB:?}/${TFOLDER:?}/EE_sroll2.py"
  #rm -f "${COB:?}/${TFOLDER:?}/EE_sroll2.bibtex"

  cppatch "${TFOLDER:?}" "EE.py" || return 1;

  cppatch "${TFOLDER:?}" "EE.yaml" || return 1;

  cppatch "${TFOLDER:?}" "TT.py" || return 1;

  cppatch "${TFOLDER:?}" "TT.yaml" || return 1;

  cppatch "${TFOLDER:?}" "EE_sroll2.py" || return 1;

  unset -v TFOLDER
  
  #-----------------------------------------------------------------------------
  # Adjust Planck 2018 lensing -------------------------------------------------
  #-----------------------------------------------------------------------------
  
  # note: we remove Cobaya native python CLIK implementation

  TFOLDER="${COBLIKE:?}/planck_2018_lensing"
  
  rm -f "${COB:?}/${TFOLDER:?}/CMBMarged.yaml"
  rm -f "${COB:?}/${TFOLDER:?}/native.yaml"

  cppatch "${TFOLDER:?}" "clik.yaml" || return 1;

  unset -v TFOLDER

  #-----------------------------------------------------------------------------
  # Adjust Planck 2018 high-ell (Plik) -----------------------------------------
  #-----------------------------------------------------------------------------
  
  # note 1: we remove Cobaya native python CLIK implementation
  # note 2: we remove unbinned likelihood for now 

  TFOLDER="${COBLIKE:?}/planck_2018_highl_plik"  # T = TMP
  
  rm -f "${COB:?}/${TFOLDER:?}/TTTEEE_unbinned.py"
  rm -f "${COB:?}/${TFOLDER:?}/TTTEEE_unbinned.yaml"
  rm -f "${COB:?}/${TFOLDER:?}/TT_unbinned.py"
  rm -f "${COB:?}/${TFOLDER:?}/TT_unbinned.yaml"
  rm -f "${COB:?}/${TFOLDER:?}/TTTEEE_lite_native.py"
  rm -f "${COB:?}/${TFOLDER:?}/TTTEEE_lite_native.yaml"
  rm -f "${COB:?}/${TFOLDER:?}/TT_lite_native.py"
  rm -f "${COB:?}/${TFOLDER:?}/TT_lite_native.yaml"

  cppatch "${TFOLDER:?}" "EE.yaml" || return 1;
  
  cppatch "${TFOLDER:?}" "TE.yaml" || return 1;
  
  cppatch "${TFOLDER:?}" "TT.yaml" || return 1;
  
  cppatch "${TFOLDER:?}" "TTTEEE.yaml" || return 1;
  
  cppatch "${TFOLDER:?}" "TT_lite.yaml" || return 1;
  
  cppatch "${TFOLDER:?}" "TTTEEE_lite.yaml" || return 1;
  
  unset -v TFOLDER

  #-----------------------------------------------------------------------------
  # SPT3G Y1 LIKELIHOOD --------------------------------------------------------
  #-----------------------------------------------------------------------------

  cppatchfolder "${COBLIKE:?}" "SPT3G_Y1" || return 1;
  
  #-----------------------------------------------------------------------------
  # ACT DR6 LENSLIKE LIKELIHOOD ------------------------------------------------
  #-----------------------------------------------------------------------------

  cppatchfolder "${COBLIKE:?}" "act_dr6_lenslike" || return 1;
  
  #-----------------------------------------------------------------------------
  # Fix renaming parameters in CAMB --------------------------------------------
  #-----------------------------------------------------------------------------
  
  cppatch "${COBTH:?}/camb" "camb.yaml" || return 1;

  #---------------------------------------------------------------------------
  # CAMSPEC ------------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  cppatch "${COBLIKE}/base_classes" "InstallableLikelihood.patch" || return 1;
  
  cdfolder "${COB:?}/${COBLIKE}/base_classes/" || return 1;
  
  patch --quiet --batch --verbose -u "InstallableLikelihood.py" -i "InstallableLikelihood.patch" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?}"; return 1; }
  
  cdfolder "${ROOTDIR:?}" || return 1;

  # ----------------------------------------------------------------------------
  # PIP COBAYA -----------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ${PIP3:?} install --editable cobaya --prefix="${ROOTDIR:?}/.local" \
    >>"${OUT1:?}" 2>>"${OUT2:?}" || { error "${EC13:?}"; return 1; }

  #-----------------------------------------------------------------------------

fi

pbottom "SETUP COBAYA" || { unset_all; return 1; }

cdfolder "${ROOTDIR:?}" || return 1;

#-----------------------------------------------------------------------------

unset_all || return 1;

#-----------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted
           
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
