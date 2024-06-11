#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_COBAYA_INSTALLATION}" ]; then
  echo -e '\033[1;44m''UPDATING COBAYA PACKAGE''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${GIT}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE GIT IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
    export OUT_UCB_1="/dev/null"
    export OUT_UCB_2="/dev/null"
  else
    export OUT_UCB_1="/dev/tty"
    export OUT_UCB_2="/dev/tty"
  fi

  export COBAYA=$ROOTDIR/cobaya/
  export COBAYA_COCOA=$ROOTDIR/../cocoa_installation_libraries/cobaya_changes
  export CBLIKE=cobaya/likelihoods
  export CBTH=cobaya/theories
  
  if [ -z "${IGNORE_COBAYA_INSTALLATION}" ]; then
    echo -e '\033[1;34m''\tINSTALLING COBAYA''\033[0m'

    #-------------------------------------------------------------------------------
    # Remove any previous installed cobaya folder
    rm -rf $ROOTDIR/cobaya

    cd $ROOTDIR
    
    export CBURL="https://github.com/CobayaSampler/cobaya.git"

    $GIT clone $CBURL cobaya > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset CBURL
      return 1
    fi
    unset CBURL
    
    cd $COBAYA

    $GIT reset --hard $COBAYA_GIT_COMMIT > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset OUT_UCB_1
      unset OUT_UCB_2
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
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (CHANGE PLANCK) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset BASECL
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

    patch -u InstallableLikelihood.py -i InstallableLikelihood.patch > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (PATCH CAMSPEC) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset BASECL
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

    $GIT clone "${NPIPE_URL}/hillipop.git" hipoptmp > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (HILLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset NPIPE_URL
      unset PL2020
      return 1
    fi

    cd $COBAYA/$CBLIKE/hipoptmp

    $GIT reset --hard $HILLIPOP_GIT_COMMIT > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (HILLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset NPIPE_URL
      unset PL2020
      return 1
    fi
    mv planck_2020_hillipop/ $COBAYA/$CBLIKE

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$PL2020/init.patch $COBAYA/$PL2020
    cd $COBAYA/$PL2020
    patch -u __init__.py -i init.patch > ${OUT_UCB_1} 2> ${OUT_UCB_2}

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

    $GIT clone "${NPIPE_URL}/lollipop.git" lipoptmp > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (LOLLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset NPIPE_URL
      unset PL2020
      return 1
    fi

    cd $COBAYA/$CBLIKE/lipoptmp
    $GIT reset --hard $LOLLIPOP_GIT_COMMIT > ${OUT_UCB_1} 2> ${OUT_UCB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (LOLLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset OUT_UCB_1
      unset OUT_UCB_2
      unset NPIPE_URL
      unset PL2020
      return 1
    fi
    mv planck_2020_lollipop/ $COBAYA/$CBLIKE

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$PL2020/init.patch $COBAYA/$PL2020
    cd $COBAYA/$PL2020
    patch -u __init__.py -i init.patch > ${OUT_UCB_1} 2> ${OUT_UCB_2}

    rm -rf $COBAYA/$CBLIKE/lipoptmp

    cd $ROOTDIR
  fi
    
  #-------------------------------------------------------------------------------
  # UNSETKEYS --------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  cd $ROOTDIR
  unset COBAYA
  unset COBAYA_COCOA
  unset CBLIKE
  unset CBTH
  unset OUT_UCB_1
  unset OUT_UCB_2
  echo -e '\033[1;44m''\e[4mUPDATING COBAYA PACKAGE DONE''\033[0m'
fi
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------