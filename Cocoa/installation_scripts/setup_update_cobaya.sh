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
  fi
  if [ -z "${GIT}" ]; then
    pfail 'GIT'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars () {
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
    unset unset_env_vars
  }
  fail () {
    export FAILMSG="\033[0;31m WE CANNOT RUN \e[3m"
    export FAILMSG2="\033[0m"
    echo -e "${FAILMSG} ${1} ${FAILMSG2}"
    unset_env_vars
    unset FAILMSG
    unset FAILMSG2
    unset fail
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
  
  if [ -z "${IGNORE_COBAYA_INSTALLATION}" ]; then
    echo -e '\033[1;34m''\tINSTALLING COBAYA''\033[0m'

    #---------------------------------------------------------------------------
    # Remove any previous installed cobaya folder
    rm -rf $ROOTDIR/cobaya

    cd $ROOTDIR
    
    export CBURL="https://github.com/CobayaSampler/cobaya.git"

    $GIT clone $CBURL cobaya > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "GIT CLONE"
      return 1
    fi
    unset CBURL
    
    cd $COBAYA

    $GIT reset --hard $COBAYA_GIT_COMMIT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "GIT CHECKOUT"
      return 1
    fi

    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # Adjust Cobaya Files ------------------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/cobaya/change_python_files.sh $COBAYA/cobaya/
    cd $COBAYA/cobaya/
    sh change_python_files.sh
    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # Adjust SN / STRONG LENSING -----------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/$CBLIKE/sn/roman_* $COBAYA/$CBLIKE/sn/
    cp -r $COBAYA_COCOA/$CBLIKE/h0licow/ $COBAYA/$CBLIKE/

    #---------------------------------------------------------------------------
    # NATIVE DES-Y1 ------------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE COBAYA NATIVE DES-Y1
    rm -rf $COBAYA/$CBLIKE/des_y1

    #---------------------------------------------------------------------------
    # 2015 ---------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE PLANCK 2015 DATA
    rm -rf $COBAYA/$CBLIKE/planck_2015_*

    #---------------------------------------------------------------------------
    # CHANGE PLANCK CLIK -------------------------------------------------------
    #---------------------------------------------------------------------------
    export BASECL="${CBLIKE}/base_classes"
    
    cp -r $COBAYA_COCOA/$BASECL/change_planck_clik.sh $COBAYA/$BASECL
    cd $COBAYA/$BASECL/
    
    sh change_planck_clik.sh
    if [ $? -ne 0 ]; then
      fail "SCRIPT CHANGE_PLANCK_CLIK"
      return 1
    fi
    unset BASECL
    
    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # PLANCK 2018 LOWELL -------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD 
    export PL_LWL="${CBLIKE}/planck_2018_lowl"
    
    rm -f $COBAYA/$PL_LWL/EE_clik.py
    rm -f $COBAYA/$PL_LWL/EE_clik.yaml
    rm -f $COBAYA/$PL_LWL/EE_sroll2.py
    rm -f $COBAYA/$PL_LWL/EE_sroll2.bibtex

    #WE RENAME TT/EE.py to also mean TT/EE_clik.py
    cp $COBAYA_COCOA/$PL_LWL/EE.py $COBAYA/$PL_LWL/EE.py
    cp $COBAYA_COCOA/$PL_LWL/EE.yaml $COBAYA/$PL_LWL/EE.yaml

    cp $COBAYA_COCOA/$PL_LWL/TT.py $COBAYA/$PL_LWL/TT.py
    cp $COBAYA_COCOA/$PL_LWL/TT.yaml $COBAYA/$PL_LWL/TT.yaml

    unset PL_LWL
    
    #---------------------------------------------------------------------------
    # PLANCK 2018 LENSING ------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD
    export PL_LL="${CBLIKE}/planck_2018_lensing"
    
    rm -f $COBAYA/$PL_LL/CMBMarged.yaml
    rm -f $COBAYA/$PL_LL/native.yaml

    # WE COPY FIXED YAML FILE
    cp $COBAYA_COCOA/$PL_LL/clik.yaml $COBAYA/$PL_LL/clik.yaml

    unset PL_LL
    
    #---------------------------------------------------------------------------
    # HIGHELL_PLIK -------------------------------------------------------------
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
    cp $COBAYA_COCOA/$HGL/TE.yaml $COBAYA/$HGL/TE.yaml

    cp $COBAYA_COCOA/$HGL/TT.yaml $COBAYA/$HGL/TT.yaml
    cp $COBAYA_COCOA/$HGL/TTTEEE.yaml $COBAYA/$HGL/TTTEEE.yaml

    cp $COBAYA_COCOA/$HGL/TT_lite.yaml $COBAYA/$HGL/TT_lite.yaml
    cp $COBAYA_COCOA/$HGL/TTTEEE_lite.yaml $COBAYA/$HGL/TTTEEE_lite.yaml

    unset HGL

    #---------------------------------------------------------------------------
    # SPT3G Y1 -----------------------------------------------------------------
    #---------------------------------------------------------------------------
    cp -r $COBAYA_COCOA/$CBLIKE/SPT3G_Y1 $COBAYA/$CBLIKE/SPT3G_Y1

    #---------------------------------------------------------------------------
    # ACT DR6 LENSLIKE ---------------------------------------------------------
    #---------------------------------------------------------------------------
    export ACTDR6_LL="${CBLIKE}/act_dr6_lenslike"
    cp -r $COBAYA_COCOA/$ACTDR6_LL $COBAYA/$ACTDR6_LL

    unset ACTDR6_LL

    #---------------------------------------------------------------------------
    # FIX RENAMING IN CAMB -----------------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/$CBTH/camb/camb.yaml $COBAYA/$CBTH/camb/camb.yaml

    #---------------------------------------------------------------------------
    cd $ROOTDIR

    echo -e '\033[1;34m''\t\e[4mINSTALLING COBAYA DONE''\033[0m'
  fi

  #-----------------------------------------------------------------------------
  # SO -------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_SIMONS_OBSERVATORY}" ]; then
    cp -r $COBAYA_COCOA/$CBLIKE/mflike $COBAYA/$CBLIKE/mflike
    cd $ROOTDIR
  fi

  #-------------------------------------------------------------------------------
  # CAMSPEC ----------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_CAMSPEC}" ]; then
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec
    
    export BASECL="${CBLIKE}/base_classes"
    
    cp $COBAYA_COCOA/$BASECL/InstallableLikelihood.patch $COBAYA/$BASECL
    cd $COBAYA/$BASECL/

    patch -u InstallableLikelihood.py -i InstallableLikelihood.patch > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "PATCH LIKELIHOOD FILES (CAMSPEC)"
      return 1
    fi

    unset BASECL
    cd $ROOTDIR
  else
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec2021
    cd $ROOTDIR
  fi

  #-----------------------------------------------------------------------------
  # INSTALL HILLIPOP -----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
    
    export PL2020="${CBLIKE}/planck_2020_hillipop"

    rm -rf $COBAYA/$PL2020
    rm -rf $COBAYA/$CBLIKE/hipoptmp

    cd $COBAYA/$CBLIKE

    export NPIPE_URL="https://github.com/planck-npipe"

    $GIT clone "${NPIPE_URL}/hillipop.git" hipoptmp > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "GIT CLONE (HILLIPOP)"
      return 1
    fi

    cd $COBAYA/$CBLIKE/hipoptmp

    $GIT reset --hard $HILLIPOP_GIT_COMMIT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "GIT RESET (HILLIPOP)"
      return 1
    fi
    mv planck_2020_hillipop/ $COBAYA/$CBLIKE

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$PL2020/init.patch $COBAYA/$PL2020
    cd $COBAYA/$PL2020
    patch -u __init__.py -i init.patch > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "PATCH LIKELIHOOD FILES (HILLIPOP)"
      return 1
    fi

    rm -rf $COBAYA/$CBLIKE/hipoptmp
    unset PL2020
    
    cd $ROOTDIR
  fi

  #-----------------------------------------------------------------------------
  # INSTALL LOWLLIPOP ----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
    
    export PL2020="${CBLIKE}/planck_2020_lollipop"

    rm -rf $COBAYA/$PL2020
    rm -rf $COBAYA/$CBLIKE/lipoptmp

    cd $COBAYA/$CBLIKE

    export NPIPE_URL="https://github.com/planck-npipe" 

    $GIT clone "${NPIPE_URL}/lollipop.git" lipoptmp > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "GIT CLONE (LOLLIPOP)"
      return 1
    fi

    cd $COBAYA/$CBLIKE/lipoptmp
    $GIT reset --hard $LOLLIPOP_GIT_COMMIT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "GIT RESET (LOLLIPOP)"
      return 1
    fi
    mv planck_2020_lollipop/ $COBAYA/$CBLIKE

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$PL2020/init.patch $COBAYA/$PL2020
    cd $COBAYA/$PL2020
    patch -u __init__.py -i init.patch > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "PATCH LIKELIHOOD FILES (LOLLIPOP)"
      return 1
    fi

    rm -rf $COBAYA/$CBLIKE/lipoptmp

    cd $ROOTDIR
  fi
    
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  unset_env_vars
  pbottom2 "UPDATING COBAYA PACKAGE"
fi
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------