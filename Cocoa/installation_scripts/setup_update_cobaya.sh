#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_COBAYA_INSTALLATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${GIT}" ]; then
    pfail 'GIT'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_ucob () {
    cd $ROOTDIR
    unset COBAYA
    unset COBAYA_COCOA
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
    unset ACTDR6_LL
    unset unset_env_vars_ucob
  }
  fail_ucob () {
    export MSG="\033[0;31m (setup_update_cobaya.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_ucob
    unset MSG
    unset MSG2
    unset fail_ucob
  }
  if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  export COBAYA=$ROOTDIR/cobaya/
  export COBAYA_COCOA=$ROOTDIR/../cocoa_installation_libraries/cobaya_changes
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
    rm -rf $ROOTDIR/cobaya

    #---------------------------------------------------------------------------
    # Clone Cobaya from original git repo --------------------------------------
    #---------------------------------------------------------------------------
    cd $ROOTDIR
    
    export CBURL="https://github.com/CobayaSampler/cobaya.git"

    $GIT clone $CBURL cobaya > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "GIT CLONE"
      return 1
    fi
    unset CBURL
    
    cd $COBAYA

    $GIT reset --hard $COBAYA_GIT_COMMIT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "GIT CHECKOUT"
      return 1
    fi

    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # Adjust Cobaya Files ------------------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/cobaya/change_python_files.sh $COBAYA/cobaya/
    if [ $? -ne 0 ]; then
      fail_ucob "CP CHANGE_PYTHON_FILES SCRIPT (COBAYA)"
      return 1
    fi

    cd $COBAYA/cobaya/
    sh change_python_files.sh
    if [ $? -ne 0 ]; then
      fail_ucob "SCRIPT CHANGE_PYTHON_FILES (COBAYA)"
      return 1
    fi

    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # Adjust SN / STRONG LENSING -----------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/$CBLIKE/sn/roman_* $COBAYA/$CBLIKE/sn/
    if [ $? -ne 0 ]; then
      fail_ucob "CP ROMAN FILES (COBAYA)"
      return 1
    fi

    cp -r $COBAYA_COCOA/$CBLIKE/h0licow/ $COBAYA/$CBLIKE/
    if [ $? -ne 0 ]; then
      fail_ucob "CP HOLICOW LIKELIHOOD (COBAYA)"
      return 1
    fi

    #---------------------------------------------------------------------------
    # Remove native DES-Y1 cobaya likelihood -----------------------------------
    #---------------------------------------------------------------------------
    rm -rf $COBAYA/$CBLIKE/des_y1

    #---------------------------------------------------------------------------
    # Remove Planck 2015 likelihood files --------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE PLANCK 2015 DATA
    rm -rf $COBAYA/$CBLIKE/planck_2015_*

    #---------------------------------------------------------------------------
    # Adjust Planck 2018 CLIK --------------------------------------------------
    #---------------------------------------------------------------------------
    export BASECL="${CBLIKE}/base_classes"
    
    cp $COBAYA_COCOA/$BASECL/change_planck_clik.sh $COBAYA/$BASECL
    if [ $? -ne 0 ]; then
      fail_ucob "CP CHANGE_PLANCK_CLIK SCRIPT"
      return 1
    fi

    cd $COBAYA/$BASECL/
    
    sh change_planck_clik.sh
    if [ $? -ne 0 ]; then
      fail_ucob "SCRIPT CHANGE_PLANCK_CLIK"
      return 1
    fi
 
    unset BASECL   
    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # Adjust Planck 2018 low-ell -----------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD 
    export PL_LWL="${CBLIKE}/planck_2018_lowl"
    
    rm -f $COBAYA/$PL_LWL/EE_clik.py
    rm -f $COBAYA/$PL_LWL/EE_clik.yaml
    rm -f $COBAYA/$PL_LWL/EE_sroll2.py
    rm -f $COBAYA/$PL_LWL/EE_sroll2.bibtex

    #WE RENAME TT/EE.py to also mean TT/EE_clik.py
    cp $COBAYA_COCOA/$PL_LWL/EE.py $COBAYA/$PL_LWL/EE.py
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (LOW-ELL EE.PY)"
      return 1
    fi

    cp $COBAYA_COCOA/$PL_LWL/EE.yaml $COBAYA/$PL_LWL/EE.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (LOW-ELL EE.YAML)"
      return 1
    fi

    cp $COBAYA_COCOA/$PL_LWL/TT.py $COBAYA/$PL_LWL/TT.py
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (LOW-ELL TT.PY)"
      return 1
    fi

    cp $COBAYA_COCOA/$PL_LWL/TT.yaml $COBAYA/$PL_LWL/TT.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (LOW-ELL TT.YAML)"
      return 1
    fi

    unset PL_LWL
    
    #---------------------------------------------------------------------------
    # Adjust Planck 2018 lensing -----------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD
    export PL_LL="${CBLIKE}/planck_2018_lensing"
    
    rm -f $COBAYA/$PL_LL/CMBMarged.yaml
    rm -f $COBAYA/$PL_LL/native.yaml

    # WE COPY FIXED YAML FILE
    cp $COBAYA_COCOA/$PL_LL/clik.yaml $COBAYA/$PL_LL/clik.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (LENSING CLIK.YAML)"
      return 1
    fi

    unset PL_LL
    
    #---------------------------------------------------------------------------
    # Adjust Planck 2018 high-ell (Plik) ---------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE UNBINNED
    export HGL="${CBLIKE}/planck_2018_highl_plik"
    
    rm -f $COBAYA/$HGL/TTTEEE_unbinned.py
    rm -f $COBAYA/$HGL/TTTEEE_unbinned.yaml
    rm -f $COBAYA/$HGL/TT_unbinned.py
    rm -f $COBAYA/$HGL/TT_unbinned.yaml

    # REMOVE LEWIS NATIVE REIMPLEMENTATION
    rm -f $COBAYA/$HGL/TTTEEE_lite_native.py
    rm -f $COBAYA/$HGL/TTTEEE_lite_native.yaml
    rm -f $COBAYA/$HGL/TT_lite_native.py
    rm -f $COBAYA/$HGL/TT_lite_native.yaml

    # FIX YAML
    cp $COBAYA_COCOA/$HGL/EE.yaml $COBAYA/$HGL/EE.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (HIGH-ELL EE.YAML)"
      return 1
    fi

    cp $COBAYA_COCOA/$HGL/TE.yaml $COBAYA/$HGL/TE.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (HIGH-ELL TE.YAML)"
      return 1
    fi

    cp $COBAYA_COCOA/$HGL/TT.yaml $COBAYA/$HGL/TT.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (HIGH-ELL TT.YAML)"
      return 1
    fi

    cp $COBAYA_COCOA/$HGL/TTTEEE.yaml $COBAYA/$HGL/TTTEEE.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (HIGH-ELL TTTEEE.YAML)"
      return 1
    fi

    cp $COBAYA_COCOA/$HGL/TT_lite.yaml $COBAYA/$HGL/TT_lite.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (HIGH-ELL TT_LITE.YAML)"
      return 1
    fi

    cp $COBAYA_COCOA/$HGL/TTTEEE_lite.yaml $COBAYA/$HGL/TTTEEE_lite.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP PLANCK 2018 LIKELIHOOD FILES (HIGH-ELL TTTEEE_LITE.YAML)"
      return 1
    fi

    unset HGL

    #---------------------------------------------------------------------------
    # SPT3G Y1 -----------------------------------------------------------------
    #---------------------------------------------------------------------------
    cp -r $COBAYA_COCOA/$CBLIKE/SPT3G_Y1 $COBAYA/$CBLIKE/SPT3G_Y1
    if [ $? -ne 0 ]; then
      fail_ucob "CP SPT-3G Y1 LIKELIHOOD FILES"
      return 1
    fi

    #---------------------------------------------------------------------------
    # ACT DR6 LENSLIKE ---------------------------------------------------------
    #---------------------------------------------------------------------------
    export ACTDR6_LL="${CBLIKE}/act_dr6_lenslike"
    
    cp -r $COBAYA_COCOA/$ACTDR6_LL $COBAYA/$ACTDR6_LL
    if [ $? -ne 0 ]; then
      fail_ucob "CP ACT-DR6 LENSING LIKELIHOOD FILES"
      return 1
    fi

    unset ACTDR6_LL

    #---------------------------------------------------------------------------
    # Fix renaming parameters in CAMB ------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/$CBTH/camb/camb.yaml $COBAYA/$CBTH/camb/camb.yaml
    if [ $? -ne 0 ]; then
      fail_ucob "CP CAMB COBAYA THEORY FILES (CAMB.YAML)"
      return 1
    fi

    #---------------------------------------------------------------------------
    cd $ROOTDIR
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

    cp -r $COBAYA_COCOA/$CBLIKE/mflike $COBAYA/$CBLIKE/mflike
    
    cd $ROOTDIR
    pbottom "INSTALLING SIMONS OBSERVATORY LIKELIHOOD"
  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # CAMSPEC --------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_CAMSPEC}" ]; then
    ptop "INSTALLING CAMSPEC LIKELIHOOD"

    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec
    
    export BASECL="${CBLIKE}/base_classes"
    
    cp $COBAYA_COCOA/$BASECL/InstallableLikelihood.patch $COBAYA/$BASECL
    if [ $? -ne 0 ]; then
      fail_ucob "CP CAMSPEC BASE LIKELIHOOD INSTALLABLELIKELIHOOD PATCH"
      return 1
    fi 

    cd $COBAYA/$BASECL/
    if [ $? -ne 0 ]; then
      fail_ucob "CD BASE CLASSES (LIKELIHOOD) FOLDER"
      return 1
    fi

    patch -u InstallableLikelihood.py -i InstallableLikelihood.patch > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "PATCH CAMSPEC BASE LIKELIHOOD INSTALLABLELIKELIHOOD.PY FILE"
      return 1
    fi

    unset BASECL
    cd $ROOTDIR
    pbottom "INSTALLING CAMSPEC LIKELIHOOD"
  else
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec2021
    cd $ROOTDIR
  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # INSTALL HILLIPOP -----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
    ptop "INSTALLING HILLIPOP LIKELIHOOD"

    export PL2020="${CBLIKE}/planck_2020_hillipop"

    rm -rf $COBAYA/$PL2020
    rm -rf $COBAYA/$CBLIKE/hipoptmp

    cd $COBAYA/$CBLIKE
    if [ $? -ne 0 ]; then
      fail_ucob "CD BASE CLASSES (LIKELIHOOD) FOLDER"
      return 1
    fi

    export NPIPE_URL="https://github.com/planck-npipe"

    $GIT clone "${NPIPE_URL}/hillipop.git" hipoptmp > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "GIT CLONE (Planck 2020 HILLIPOP)"
      return 1
    fi

    cd $COBAYA/$CBLIKE/hipoptmp
    if [ $? -ne 0 ]; then
      fail_ucob "CD LIKELIHOOD TMP FOLDER (Planck 2020 HILLIPOP)"
      return 1
    fi

    $GIT reset --hard $HILLIPOP_GIT_COMMIT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "GIT RESET (Planck 2020 HILLIPOP)"
      return 1
    fi

    mv planck_2020_hillipop/ $COBAYA/$CBLIKE
    if [ $? -ne 0 ]; then
      fail_ucob "MV LIKELIHOOD FOLDER (Planck 2020 HILLIPOP)"
      return 1
    fi
    
    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$PL2020/init.patch $COBAYA/$PL2020
    if [ $? -ne 0 ]; then
      fail_ucob "CP INIT.PATCH (Planck 2020 HILLIPOP)"
      return 1
    fi

    cd $COBAYA/$PL2020
    if [ $? -ne 0 ]; then
      fail_ucob "CD LIKELIHOOD FOLDER (Planck 2020 HILLIPOP)"
      return 1
    fi

    patch -u __init__.py -i init.patch > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "PATCH LIKELIHOOD FILES (Planck 2020 HILLIPOP)"
      return 1
    fi

    rm -rf $COBAYA/$CBLIKE/hipoptmp
    unset PL2020
    cd $ROOTDIR
    pbottom "INSTALLING HILLIPOP LIKELIHOOD"
  fi

  #-----------------------------------------------------------------------------
  # INSTALL LOWLLIPOP ----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
    ptop "INSTALLING LOLLIPOP LIKELIHOOD"

    export PL2020="${CBLIKE}/planck_2020_lollipop"

    rm -rf $COBAYA/$PL2020
    rm -rf $COBAYA/$CBLIKE/lipoptmp

    cd $COBAYA/$CBLIKE
    if [ $? -ne 0 ]; then
      fail_ucob "CD BASE CLASSES (LIKELIHOOD) FOLDER"
      return 1
    fi

    export NPIPE_URL="https://github.com/planck-npipe" 

    $GIT clone "${NPIPE_URL}/lollipop.git" lipoptmp > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "GIT CLONE (Planck 2020 LOLLIPOP)"
      return 1
    fi

    cd $COBAYA/$CBLIKE/lipoptmp
    if [ $? -ne 0 ]; then
      fail_ucob "CD LIKELIHOOD TMP FOLDER (Planck 2020 LOLLIPOP)"
      return 1
    fi

    $GIT reset --hard $LOLLIPOP_GIT_COMMIT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "GIT RESET (Planck 2020 LOLLIPOP)"
      return 1
    fi

    mv planck_2020_lollipop/ $COBAYA/$CBLIKE
    if [ $? -ne 0 ]; then
      fail_ucob "MV LIKELIHOOD FOLDER (Planck 2020 LOLLIPOP)"
      return 1
    fi

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$PL2020/init.patch $COBAYA/$PL2020
    if [ $? -ne 0 ]; then
      fail_ucob "CP INIT.PATCH (Planck 2020 LOLLIPOP)"
      return 1
    fi

    cd $COBAYA/$PL2020
    if [ $? -ne 0 ]; then
      fail_ucob "CD LIKELIHOOD FOLDER (Planck 2020 LOLLIPOP)"
      return 1
    fi

    patch -u __init__.py -i init.patch > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_ucob "PATCH LIKELIHOOD FILES (Planck 2020 LOLLIPOP)"
      return 1
    fi

    rm -rf $COBAYA/$CBLIKE/lipoptmp
    unset PL2020
    cd $ROOTDIR

    pbottom "INSTALLING LOLLIPOP LIKELIHOOD"
  fi
    
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  unset_env_vars_ucob
  pbottom2 "UPDATING COBAYA PACKAGE"
fi
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------