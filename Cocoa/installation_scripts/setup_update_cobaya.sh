#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_COBAYA_INSTALLATION}" ]; then
  pfail() {
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi

  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }

  if [ -z "${GIT}" ]; then
    pfail 'GIT'; cdroot; return 1
  fi

  unset_env_vars_ucob () {
    unset COB
    unset COCOACOB
    unset CBLIKE
    unset CBTH
    unset OUT1
    unset OUT2
    unset CBURL
    unset PL2020
    unset PL_LL
    unset BASECL
    unset NPIPE_URL
    unset HGL
    unset ftmp
    unset ACTDR6_LL
    unset unset_env_vars_ucob
    cdroot || return 1;
  }

  fail_ucb () {
    local MSG="\033[0;31m\t\t (setup_update_cobaya.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_ucb
    unset_env_vars_ucob
  }

  if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_ucb "CD FOLDER: ${1}"; return 1; }
  }
  
  export COB="${ROOTDIR}/cobaya/"
  export COCOACOB="${ROOTDIR}/../cocoa_installation_libraries/cobaya_changes"
  export CBLIKE=cobaya/likelihoods
  export CBTH=cobaya/theories

  ptop2 "UPDATING COBAYA PACKAGE"
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_COBAYA_INSTALLATION}" ]; then
    
    ptop "INSTALLING COBAYA"

    #---------------------------------------------------------------------------
    # Remove any previous installed cobaya folder ------------------------------
    #---------------------------------------------------------------------------
    rm -rf "${ROOTDIR}/cobaya"

    #---------------------------------------------------------------------------
    # Clone Cobaya from original git repo --------------------------------------
    #---------------------------------------------------------------------------
    cdfolder "${ROOTDIR}" || return 1;
    
    export CBURL="https://github.com/CobayaSampler/cobaya.git"

    $GIT clone $CBURL cobaya >${OUT1} 2>${OUT2} ||
      { fail_ucb "GIT CLONE"; return 1 }
    
    unset CBURL
    
    cdfolder "${COB}" || return 1;

    $GIT reset --hard ${COBAYA_GIT_COMMIT} >${OUT1} 2>${OUT2} ||
      { fail_ucb "GIT CHECKOUT"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    #---------------------------------------------------------------------------
    # Adjust Cobaya Files ------------------------------------------------------
    #---------------------------------------------------------------------------
    cp ${COCOACOB}/cobaya/change_python_files.sh ${COB}/cobaya/
    if [ $? -ne 0 ]; then
      fail_ucb "CP CHANGE_PYTHON_FILES SCRIPT (COBAYA)"
      return 1
    fi

    cdfolder "${COB}/cobaya/" || return 1;
    
    sh change_python_files.sh ||
      { fail_ucb "SCRIPT CHANGE_PYTHON_FILES (COBAYA)"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    #---------------------------------------------------------------------------
    # Adjust SN / STRONG LENSING -----------------------------------------------
    #---------------------------------------------------------------------------
    cp "${COCOACOB}/${CBLIKE}"/sn/roman_* "${COB}/${CBLIKE}"/sn/ \
      2>${OUT2} || { fail_ucb "CP ROMAN FILES (COBAYA)"; return 1; }
  
    cp -r "${COCOACOB}/${CBLIKE}"/h0licow/ "${COB}/${CBLIKE}/" \
      2>${OUT2} || { fail_ucb "CP HOLICOW LIKELIHOOD (COBAYA)"; return 1; }
  
    #---------------------------------------------------------------------------
    # Remove native DES-Y1 cobaya likelihood -----------------------------------
    #---------------------------------------------------------------------------
    rm -rf "${COB}/${CBLIKE}/des_y1"

    #---------------------------------------------------------------------------
    # Remove Planck 2015 likelihood files --------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE PLANCK 2015 DATA
    rm -rf "${COB}/${CBLIKE}"/planck_2015_*

    #---------------------------------------------------------------------------
    # Adjust Planck 2018 CLIK --------------------------------------------------
    #---------------------------------------------------------------------------
    export BASECL="${CBLIKE}/base_classes"
    
    cp "${COCOACOB}/${BASECL}"/change_planck_clik.sh "${COB}/${BASECL}" \
      2>${OUT2} || { fail_ucb "CP CHANGE_PLANCK_CLIK SCRIPT"; return 1; }
  
    cdfolder "${COB}/${BASECL}/" || return 1;
    
    sh change_planck_clik.sh || 
      { fail_ucb "SCRIPT CHANGE_PLANCK_CLIK"; return 1; }
    
    unset BASECL

    cdfolder "${ROOTDIR}" || return 1;

    #---------------------------------------------------------------------------
    # Adjust Planck 2018 low-ell -----------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD 
    export PLLWL="${CBLIKE}/planck_2018_lowl"
    
    rm -f "${COB}/${PLLWL}"/EE_clik.py
    rm -f "${COB}/${PLLWL}"/EE_clik.yaml
    rm -f "${COB}/${PLLWL}"/EE_sroll2.py
    rm -f "${COB}/${PLLWL}"/EE_sroll2.bibtex

    #We rename TT/EE.py to also mean TT/EE_clik.py (in the original cobaya)
    cp "${COCOACOB}/${PLLWL}"/EE.py "${COB}/${PLLWL}/EE.py" 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE LOW-ELL EE.PY"; return 1; }
    
    cp "${COCOACOB}/${PLLWL}"/EE.yaml "${COB}/${PLLWL}/EE.yaml" 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE LOW-ELL EE.YAML"; return 1; }

    cp "${COCOACOB}/${PLLWL}"/TT.py "${COB}/${PLLWL}/TT.py" 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE LOW-ELL TT.PY"; return 1; }
    
    cp "${COCOACOB}/${PLLWL}"/TT.yaml "${COB}/${PLLWL}/TT.yaml" 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE LOW-ELL TT.YAML"; return 1; }

    unset PLLWL
    
    #---------------------------------------------------------------------------
    # Adjust Planck 2018 lensing -----------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD
    export PL_LL="${CBLIKE}/planck_2018_lensing"
    
    rm -f "${COB}/${PL_LL}"/CMBMarged.yaml
    rm -f "${COB}/${PL_LL}"/native.yaml

    # WE COPY FIXED YAML FILE
    cp "${COCOACOB}/${PL_LL}"/clik.yaml "${COB}/${PL_LL}"/clik.yaml 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE LENSING CLIK.YAML"; return 1; }

    unset PL_LL
    
    #---------------------------------------------------------------------------
    # Adjust Planck 2018 high-ell (Plik) ---------------------------------------
    #---------------------------------------------------------------------------

    # WE REMOVE UNBINNED LIKELIHOOD FOR NOW
    export HGL="${CBLIKE}/planck_2018_highl_plik"
    
    rm -f "${COB}/${HGL}"/TTTEEE_unbinned.py
    rm -f "${COB}/${HGL}"/TTTEEE_unbinned.yaml
    rm -f "${COB}/${HGL}"/TT_unbinned.py
    rm -f "${COB}/${HGL}"/TT_unbinned.yaml

    # REMOVE LEWIS NATIVE REIMPLEMENTATION
    rm -f "${COB}/${HGL}"/TTTEEE_lite_native.py
    rm -f "${COB}/${HGL}"/TTTEEE_lite_native.yaml
    rm -f "${COB}/${HGL}"/TT_lite_native.py
    rm -f "${COB}/${HGL}"/TT_lite_native.yaml

    # FIX YAML
    cp "${COCOACOB}/${HGL}"/EE.yaml "${COB}/${HGL}"/EE.yaml 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE HIGH-ELL EE.YAML"; return 1; }
    
    cp "${COCOACOB}/${HGL}"/TE.yaml "${COB}/${HGL}"/TE.yaml 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE HIGH-ELL TE.YAML"; return 1; }
    
    cp "${COCOACOB}/${HGL}"/TT.yaml "${COB}/${HGL}"/TT.yaml 2>${OUT2} ||
      { fail_ucb "CP PLANCK2018 LIKELIHOOD FILE HIGH-ELL TT.YAML"; return 1; }
  
    cp "${COCOACOB}/${HGL}"/TTTEEE.yaml "${COB}/${HGL}"/TTTEEE.yaml \
      2>${OUT2} || { fail_ucb \
      "CP PLANCK2018 LIKELIHOOD FILE HIGH-ELL TTTEEE.YAML"; return 1; }
  
    cp "${COCOACOB}/${HGL}"/TT_lite.yaml "${COB}/${HGL}"/TT_lite.yaml \
      2>${OUT2} || { fail_ucb \
      "CP PLANCK2018 LIKELIHOOD FILE HIGH-ELL TT_LITE.YAML"; return 1; }
    
    cp "${COCOACOB}/${HGL}"/TTTEEE_lite.yaml "${COB}/${HGL}"/TTTEEE_lite.yaml \
      2>${OUT2} || { fail_ucb \
      "CP PLANCK2018 LIKELIHOOD FILE HIGH-ELL TTTEEE_LITE.YAML"; return 1; }
    
    unset HGL

    #---------------------------------------------------------------------------
    # SPT3G Y1 -----------------------------------------------------------------
    #---------------------------------------------------------------------------

    cp -r "${COCOACOB}/${CBLIKE}"/SPT3G_Y1 "${COB}/${CBLIKE}"/SPT3G_Y1 \
      2>${OUT2} || { fail_ucb "CP SPT-3G Y1 LIKELIHOOD FILES"; return 1; }
    

    #---------------------------------------------------------------------------
    # ACT DR6 LENSLIKE ---------------------------------------------------------
    #---------------------------------------------------------------------------

    export ACTDR6_LL="${CBLIKE}/act_dr6_lenslike"
    
    cp -r "${COCOACOB}/$ACTDR6_LL" "${COB}/$ACTDR6_LL" 2>${OUT2} ||
      { fail_ucb "CP ACT-DR6 LENSING LIKELIHOOD FILES"; return 1; }
    
    unset ACTDR6_LL

    #---------------------------------------------------------------------------
    # Fix renaming parameters in CAMB ------------------------------------------
    #---------------------------------------------------------------------------

    cp "${COCOACOB}/${CBTH}"/camb/camb.yaml "${COB}/${CBTH}"/camb/camb.yaml \
      2>${OUT2} || { fail_ucb "CP CAMB COBAYA THEORY FILES (CAMB.YAML)"; \
      return 1; }
  
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------    

    cdfolder "${ROOTDIR}" || return 1;

    pbottom "INSTALLING COBAYA"

  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # SO -------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # TODO - download from original git repo
  if [ -z "${SKIP_DECOMM_SIMONS_OBSERVATORY}" ]; then
    
    ptop "INSTALLING SIMONS OBSERVATORY LIKELIHOOD"

    cp -r "${COCOACOB}/${CBLIKE}"/mflike "${COB}/${CBLIKE}"/mflike 2>${OUT2} ||
      { fail_ucb "CP SIMONS OBSERVATORY LIKELIHOOD FILES"; return 1; }
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom "INSTALLING SIMONS OBSERVATORY LIKELIHOOD"

  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # CAMSPEC --------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_CAMSPEC}" ]; then

    ptop "INSTALLING CAMSPEC LIKELIHOOD"

    rm -rf "${COB}/${CBLIKE}"/planck_2018_highl_CamSpec
    
    export BASECL="${CBLIKE}/base_classes"
    
    cp "${COCOACOB}/${BASECL}"/InstallableLikelihood.patch \
      "${COB}/${BASECL}" 2>${OUT2} || { fail_ucb \
      "CP CAMSPEC BASE LIKELIHOOD INSTALLABLELIKELIHOOD PATCH"; return 1; }
    
    cdfolder "${COB}/${BASECL}/" || return 1;

    patch -u InstallableLikelihood.py \
      -i InstallableLikelihood.patch >${OUT1} 2>${OUT2} || 
      { fail_ucb "PATCH CAMSPEC BASE LIKELIHOOD INSTALLABLELIKELIHOOD.PY FILE"; \
        return 1; }
    
    unset BASECL
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom "INSTALLING CAMSPEC LIKELIHOOD"

  else

    rm -rf "${COB}/${CBLIKE}"/planck_2018_highl_CamSpec
    rm -rf "${COB}/${CBLIKE}"/planck_2018_highl_CamSpec2021

  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # INSTALL HILLIPOP -----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
    
    ptop "INSTALLING HILLIPOP LIKELIHOOD"

    export PL2020="${CBLIKE}/planck_2020_hillipop"
    export ftmp=hipoptmp

    rm -rf "${COB}/${PL2020}"
    rm -rf "${COB}/${CBLIKE}/${ftmp}"

    cdfolder "${COB}/${CBLIKE}" || return 1;

    export NPIPE_URL="https://github.com/planck-npipe"

    $GIT clone "${NPIPE_URL}/hillipop.git" "${ftmp}" >${OUT1} 2>${OUT2} || 
      { fail_ucb "GIT CLONE (Planck 2020HILLIPOP)"; return 1; }
  
    cdfolder "${COB}/${CBLIKE}/${ftmp}" || return 1;

    $GIT reset --hard $HILLIPOP_GIT_COMMIT >${OUT1} 2>${OUT2} ||
      { fail_ucb "GIT RESET (Planck2020 HILLIPOP)"; return 1; |
  
    mv planck_2020_hillipop/ "${COB}/${CBLIKE}" 2>${OUT2} ||
      { fail_ucb "MV LIKELIHOOD FOLDER (Planck2020 HILLIPOP)"; return 1; }
    
    #---------------------------------------------------------------------------
    # now patch the likelihood __init__ file
    #---------------------------------------------------------------------------
    cp "${COCOACOB}/${PL2020}/init.patch" "${COB}/${PL2020}" 2>${OUT2} ||
      { fail_ucb "CP INIT.PATCH (Planck2020 HILLIPOP)"; return 1; }
  
    cdfolder "${COB}/${PL2020}" || return 1;

    patch -u __init__.py -i init.patch >${OUT1} 2>${OUT2} || 
      { fail_ucb "PATCH LIKELIHOOD FILES (Planck2020 HILLIPOP)"; return 1; }
    
    rm -rf "${COB}/${CBLIKE}/${ftmp}"
    
    unset PL2020
    unset ftmp

    cdfolder "${ROOTDIR}" || return 1;
    
    pbottom "INSTALLING HILLIPOP LIKELIHOOD"

  fi

  #-----------------------------------------------------------------------------
  # INSTALL LOWLLIPOP ----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then

    ptop "INSTALLING LOLLIPOP LIKELIHOOD"

    export PL2020="${CBLIKE}/planck_2020_lollipop"
    export ftmp=lipoptmp

    rm -rf "${COB}/${PL2020}"
    rm -rf "${COB}/${CBLIKE}/${ftmp}"

    cdfolder "${COB}/${CBLIKE}" || return 1;

    export NPIPE_URL="https://github.com/planck-npipe" 

    $GIT clone "${NPIPE_URL}/lollipop.git" ${ftmp} >${OUT1} 2>${OUT2} || 
      { fail_ucb "GIT CLONE (Planck2020 LOLLIPOP)"; return 1; }
    
    cdfolder "${COB}/${CBLIKE}/${ftmp}" || return 1;

    $GIT reset --hard $LOLLIPOP_GIT_COMMIT >${OUT1} 2>${OUT2} || 
      { fail_ucb "GIT RESET (Planck2020 LOLLIPOP)"; return 1; }
    
    mv planck_2020_lollipop/ "${COB}/${CBLIKE}" 2>${OUT2} || 
      { fail_ucb "MV LIKELIHOOD FOLDER (Planck2020 LOLLIPOP)"; return 1; }
    
    #---------------------------------------------------------------------------
    # now patch the likelihood __init__ file
    #---------------------------------------------------------------------------
    cp "${COCOACOB}/${PL2020}/init.patch" "${COB}/${PL2020}" 2>${OUT2} || 
      { fail_ucb "CP INIT.PATCH (Planck2020 LOLLIPOP)"; return 1; }

    cdfolder "${COB}/${PL2020}" || return 1;

    patch -u __init__.py -i init.patch >${OUT1} 2>${OUT2} ||
      { fail_ucb "PATCH LIKELIHOOD FILES (Planck2020 LOLLIPOP)"; return 1; }
    
    rm -rf "${COB}/${CBLIKE}/lipoptmp"
    
    unset PL2020
    unset ftmp
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom "INSTALLING LOLLIPOP LIKELIHOOD"

  fi
    
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  
  unset_env_vars_ucob | return 1;
  
  pbottom2 "UPDATING COBAYA PACKAGE"

fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------