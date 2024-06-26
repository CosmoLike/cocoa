#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_COBAYA_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR:?}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;

  unset_env_vars () {
    unset -v COB CCCOB COBLIKE COBTH URL TFILE TFOLDER TFOLDER2
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile flipop cppatch cppatchfolder
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
    fail_script_msg "$(basename ${BASH_SOURCE[0]})" "${1}"
    unset_all || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  cpfolder() {
    cp -r "${1:?}" "${2:?}"  \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  cpfile() {
    cp "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
#  ptop2 "UPDATING COBAYA PACKAGE" || return 1

  unset_env_vars || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

  COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

  CCCOB="${CCIL:?}/cobaya_changes"  # CC = CoCoA, COB = Cobaya (Cocoa Cobaya)
  
  COBLIKE="cobaya/likelihoods"      # COB = Cobaya, LIKE = likelihoods
  
  COBTH="cobaya/theories"           # COB = Cobaya, TH = theory

  cppatch() {
    cp "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" \
      2>"/dev/null" || 
      { error "CP FILE ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; return 1; }
  }

  cppatchfolder() {
    cp -r "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" \
      2>"/dev/null" || 
      { error "CP FOLDER ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; \
      return 1; }
  }

  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_COBAYA_INSTALLATION}" ]; then
    
    ptop "GETTING AND PATCHING COBAYA" || return 1;

    #---------------------------------------------------------------------------
    # Remove any previous installed cobaya folder ------------------------------
    #---------------------------------------------------------------------------
    rm -rf "${ROOTDIR:?}/cobaya"

    #---------------------------------------------------------------------------
    # Clone Cobaya from original git repo --------------------------------------
    #---------------------------------------------------------------------------
    cdfolder "${ROOTDIR:?}" || return 1;
    
    URL="${COBAYA_URL:-"https://github.com/CobayaSampler/cobaya.git"}"

    "${GIT:?}" clone "${URL:?}" cobaya \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    unset URL
    
    cdfolder "${COB:?}" || return 1;

    if [ -n "${COBAYA_GIT_COMMIT}" ]; then
      ${GIT:?} reset --hard "${COBAYA_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

    cdfolder "${ROOTDIR:?}" || return 1;

    #---------------------------------------------------------------------------
    # Adjust Cobaya Files ------------------------------------------------------
    #---------------------------------------------------------------------------
    
    TFILE="change_python_files.sh" # T = TMP
    
    cppatch "cobaya" "${TFILE:?}" || return 1

    cdfolder "${COB:?}/cobaya/" || return 1;
    
    # parenthesis = run in a subshell
    ( sh change_python_files.sh ) || { error "${EC22:?} (CPF)"; return 1; }

    cdfolder "${ROOTDIR:?}" || return 1;

    unset -v TFILE

    #---------------------------------------------------------------------------
    # Adjust SN / STRONG LENSING -----------------------------------------------
    #---------------------------------------------------------------------------
    
    TFOLDER="${COBLIKE:?}/sn" # T = TMP

    cp "${CCCOB:?}/${TFOLDER:?}/"roman_* "${COB:?}/${TFOLDER:?}/" \
      2>${OUT2:?} || { error "CP ROMAN FILES"; return 1; }
  
    cppatchfolder "${COBLIKE:?}" "h0licow" || return 1
  
    unset -v TFOLDER

    #---------------------------------------------------------------------------
    # Remove native DES-Y1 cobaya likelihood -----------------------------------
    #---------------------------------------------------------------------------
    
    rm -rf "${COB:?}/${COBLIKE:?}/des_y1"

    #---------------------------------------------------------------------------
    # Remove Planck 2015 likelihood files --------------------------------------
    #---------------------------------------------------------------------------
    
    rm -rf "${COB:?}/${COBLIKE:?}"/planck_2015_*

    #---------------------------------------------------------------------------
    # Adjust Planck 2018 CLIK --------------------------------------------------
    #---------------------------------------------------------------------------
    
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
    rm -f "${COB:?}/${TFOLDER:?}/EE_sroll2.bibtex"

    cppatch "${TFOLDER:?}" "EE.py" || return 1

    cppatch "${TFOLDER:?}" "EE.yaml" || return 1

    cppatch "${TFOLDER:?}" "TT.py" || return 1

    cppatch "${TFOLDER:?}" "TT.yaml" || return 1

    unset -v TFOLDER
    
    #---------------------------------------------------------------------------
    # Adjust Planck 2018 lensing -----------------------------------------------
    #---------------------------------------------------------------------------
    
    # note: we remove Cobaya native python CLIK implementation

    TFOLDER="${COBLIKE:?}/planck_2018_lensing"
    
    rm -f "${COB:?}/${TFOLDER:?}/CMBMarged.yaml"
    rm -f "${COB:?}/${TFOLDER:?}/native.yaml"

    cppatch "${TFOLDER:?}" "clik.yaml" || return 1

    unset -v TFOLDER

    #---------------------------------------------------------------------------
    # Adjust Planck 2018 high-ell (Plik) ---------------------------------------
    #---------------------------------------------------------------------------
    
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

    cppatch "${TFOLDER:?}" "EE.yaml" || return 1
    
    cppatch "${TFOLDER:?}" "TE.yaml" || return 1
    
    cppatch "${TFOLDER:?}" "TT.yaml" || return 1
    
    cppatch "${TFOLDER:?}" "TTTEEE.yaml" || return 1
    
    cppatch "${TFOLDER:?}" "TT_lite.yaml" || return 1
    
    cppatch "${TFOLDER:?}" "TTTEEE_lite.yaml" || return 1
    
    unset -v TFOLDER

    #---------------------------------------------------------------------------
    # SPT3G Y1 LIKELIHOOD ------------------------------------------------------
    #---------------------------------------------------------------------------

    cppatchfolder "${COBLIKE:?}" "SPT3G_Y1" || return 1
    
    #---------------------------------------------------------------------------
    # ACT DR6 LENSLIKE LIKELIHOOD ----------------------------------------------
    #---------------------------------------------------------------------------

    cppatchfolder "${COBLIKE:?}" "act_dr6_lenslike" || return 1
    
    #---------------------------------------------------------------------------
    # Fix renaming parameters in CAMB ------------------------------------------
    #---------------------------------------------------------------------------
    
    cppatch "${COBTH:?}/camb" "camb.yaml" || return 1

    # --------------------------------------------------------------------------
    # PIP COBAYA ---------------------------------------------------------------
    # --------------------------------------------------------------------------
    
    ${PIP3:?} install --editable cobaya --prefix="${ROOTDIR:?}/.local" \
      >"${OUT1:?}" 2>"${OUT2:?}" || { error "${EC13:?}"; return 1; }

    #---------------------------------------------------------------------------

    cdfolder "${ROOTDIR:?}" || return 1;

    pbottom "GETTING AND PATCHING COBAYA" || return 1;

  fi

  #-----------------------------------------------------------------------------
  # INSTALL SIMONS OBSERVATORY LIKELIHOOD --------------------------------------
  #-----------------------------------------------------------------------------
  # TODO - download from original git repo
   
  if [ -z "${IGNORE_SO_INSTALLATION}" ]; then

    ptop "GETTING AND PATCHING SIMONS OBSERVATORY LIKELIHOOD" || return 1;

    cppatchfolder "${COBLIKE:?}" "mflike" || return 1
    
    cdfolder "${ROOTDIR:?}" || return 1;

    pbottom "GETTING AND PATCHING SIMONS OBSERVATORY LIKELIHOOD" || return 1;

  fi

  #-----------------------------------------------------------------------------
  # INSTALL CAMSPEC 2021 LIKELIHOOD --------------------------------------------
  #-----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CAMSPEC_INSTALLATION}" ]; then

    # note: we always delete CAMSPEC2018 likelihood
    rm -rf "${COB:?}/${COBLIKE:?}"/planck_2018_highl_CamSpec

    ptop "PATCHING CAMSPEC 2021 LIKELIHOOD" || return 1;

    TFOLDER="${COBLIKE}/base_classes"
    
    TFILE="InstallableLikelihood"
    
    cppatch "${TFOLDER:?}" "${TFILE:?}.patch"
    
    cdfolder "${COB:?}/${TFOLDER}/" || return 1;

    patch -u "${TFILE:?}.py" -i "${TFILE:?}.patch" >${OUT1:?} 2>${OUT2:?} || 
      { error "${EC17:?}"; return 1; }
    
    unset -v TFOLDER TFILE
    
    cdfolder "${ROOTDIR:?}" || return 1;

    pbottom "PATCHING CAMSPEC 2021 LIKELIHOOD" || return 1;

  fi
  
  #-----------------------------------------------------------------------------
  # INSTALL LIPOP LIKELIHOOD ---------------------------------------------------
  #-----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_LIPOP_INSTALLATION}" ]; then

    flipop() {
      local TFOLDER="${COBLIKE:?}/${1:?}"
      local URL="${2:?}"

      # in case you run this script more than once
      rm -rf "${COB:?}/${TFOLDER:?}"
      rm -rf "${COB:?}/${COBLIKE:?}/tmp"

      cdfolder "${COB:?}/${COBLIKE:?}" || return 1;

      ${GIT:?} clone "${URL:?}" "tmp" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
      cdfolder "${COB:?}/${COBLIKE}/tmp" || return 1;

      "${GIT:?}" reset --hard "${3:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC23:?}"; return 1; }

      cpfolder ${1:?}/ "${COB:?}/${COBLIKE:?}" || return 1;
      
      cdfolder "${COB:?}/${TFOLDER:?}" || return 1;
      
      rm -rf "${COB:?}/${COBLIKE:?}/tmp"

      cppatch "${TFOLDER:?}" "init.patch" || return 1

      patch -u '__init__.py' -i 'init.patch' >${OUT1:?} 2>${OUT2:?} || 
        { error "${EC17:?}"; return 1; }

      cdfolder "${ROOTDIR:?}" || return 1; 
    }

    #---------------------------------------------------------------------------
    # HILLIPOP_URL LIKELIHOOD --------------------------------------------------
    #---------------------------------------------------------------------------

    ptop "GETTING AND PATCHING HILLIPOP LIKELIHOOD" || return 1

    TFOLDER="planck_2020_hillipop"
    
    URL="${HILLIPOP_URL:-"https://github.com/planck-npipe/hillipop.git"}"
    
    flipop "${TFOLDER:?}" "${URL:?}" "${HILLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;
    
    unset -v TFOLDER URL

    pbottom "GETTING AND PATCHING HILLIPOP LIKELIHOOD" || return 1

    #---------------------------------------------------------------------------
    # LOLLIPOP LIKELIHOOD ------------------------------------------------------
    #---------------------------------------------------------------------------

    ptop "GETTING AND PATCHING LOLLIPOP LIKELIHOOD" || return 1

    TFOLDER="planck_2020_lollipop"
    
    URL="${LOLLIPOP_URL:-"https://github.com/planck-npipe/lollipop.git"}"
    
    flipop "${TFOLDER:?}" "${URL:?}" "${LOLLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;
    
    unset -v TFOLDER URL

    pbottom "GETTING AND PATCHING LOLLIPOP LIKELIHOOD" || return 1

  fi

  #-----------------------------------------------------------------------------
  
  unset_all || return 1;
  
#  pbottom2 "UPDATING COBAYA PACKAGE" || return 1;

fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------