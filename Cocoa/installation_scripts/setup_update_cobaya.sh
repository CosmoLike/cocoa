#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_COBAYA_INSTALLATION}" ]; then
  echo -e '\033[1;44m''Updating Cobaya Package''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE CXX_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE C_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PIP3}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE PIP3 IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
    export OUT_UPT_CB_1="/dev/null"
    export OUT_UPT_CB_2="/dev/null"
  else
    export OUT_UPT_CB_1="/dev/tty"
    export OUT_UPT_CB_2="/dev/tty"
  fi

  export COBAYA=$ROOTDIR/cobaya/
  export COBAYA_COCOA=$ROOTDIR/../cocoa_installation_libraries/cobaya_changes
  export CBLIKE=cobaya/likelihoods
  export CBTH=cobaya/theories
  export HGL=planck_2018_highl_plik

  cd $ROOTDIR

  if [ -z "${IGNORE_COBAYA_INSTALLATION}" ]; then
    #---------------------------------------------------------------------------
    # Remove any previous installed cobaya folder ------------------------------
    #---------------------------------------------------------------------------
    rm -rf $ROOTDIR/cobaya
    
    #---------------------------------------------------------------------------
    # CLONE COBAYA FROM ORIGINAL REPO ------------------------------------------
    #---------------------------------------------------------------------------
    git clone https://github.com/CobayaSampler/cobaya.git cobaya > ${OUT_UPT_CB_1} 2> ${OUT_UPT_CB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi

    cd $COBAYA

    git reset --hard $COBAYA_GIT_COMMIT
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
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
    # Native DES-Y1 ------------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE COBAYA NATIVE DES-Y1
    rm -rf $COBAYA/$CBLIKE/des_y1
    
    #---------------------------------------------------------------------------
    # Planck 2015 --------------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE PLANCK 2015 DATA
    rm -rf $COBAYA/$CBLIKE/planck_2015_*

    #---------------------------------------------------------------------------
    # CHANGE PLANCK CLIK -------------------------------------------------------
    #---------------------------------------------------------------------------
    cp -r $COBAYA_COCOA/$CBLIKE/base_classes/change_planck_clik.sh $COBAYA/$CBLIKE/base_classes/

    cd $COBAYA/$CBLIKE/base_classes/
    sh change_planck_clik.sh
    if [ $? -ne 0 ]; then
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi
    cd $ROOTDIR

    #---------------------------------------------------------------------------
    # PLANCK 2018 LOWELL -------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD 
    rm -f $COBAYA/$CBLIKE/planck_2018_lowl/EE_clik.py
    rm -f $COBAYA/$CBLIKE/planck_2018_lowl/EE_clik.yaml
    rm -f $COBAYA/$CBLIKE/planck_2018_lowl/EE_sroll2.py
    rm -f $COBAYA/$CBLIKE/planck_2018_lowl/EE_sroll2.bibtex

    #WE RENAME TT/EE.py to also mean TT/EE_clik.py
    cp $COBAYA_COCOA/$CBLIKE/planck_2018_lowl/EE.py $COBAYA/$CBLIKE/planck_2018_lowl/EE.py
    cp $COBAYA_COCOA/$CBLIKE/planck_2018_lowl/EE.yaml $COBAYA/$CBLIKE/planck_2018_lowl/EE.yaml

    cp $COBAYA_COCOA/$CBLIKE/planck_2018_lowl/TT.py $COBAYA/$CBLIKE/planck_2018_lowl/TT.py
    cp $COBAYA_COCOA/$CBLIKE/planck_2018_lowl/TT.yaml $COBAYA/$CBLIKE/planck_2018_lowl/TT.yaml

    #---------------------------------------------------------------------------
    # PLANCK 2018 LENSING ------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE LEWIS NATIVE LIKELIHOOD
    rm -f $COBAYA/$CBLIKE/planck_2018_lensing/CMBMarged.yaml
    rm -f $COBAYA/$CBLIKE/planck_2018_lensing/native.yaml

    # WE COPY FIXED YAML FILE
    cp $COBAYA_COCOA/$CBLIKE/planck_2018_lensing/clik.yaml $COBAYA/$CBLIKE/planck_2018_lensing/clik.yaml

    #---------------------------------------------------------------------------
    # HIGHELL_PLIK -------------------------------------------------------------
    #---------------------------------------------------------------------------
    # WE REMOVE UNBINNED
    rm -f $COBAYA/$CBLIKE/$HGL/TTTEEE_unbinned.py
    rm -f $COBAYA/$CBLIKE/$HGL/TTTEEE_unbinned.yaml
    rm -f $COBAYA/$CBLIKE/$HGL/TT_unbinned.py
    rm -f $COBAYA/$CBLIKE/$HGL/TT_unbinned.yaml

    # REMOVE LEWIS NATIVE REIMPLEMENTATION
    rm -f $COBAYA/$CBLIKE/$HGL/TTTEEE_lite_native.py
    rm -f $COBAYA/$CBLIKE/$HGL/TTTEEE_lite_native.yaml

    rm -f $COBAYA/$CBLIKE/$HGL/TT_lite_native.py
    rm -f $COBAYA/$CBLIKE/$HGL/TT_lite_native.yaml

    # FIX YAML
    cp $COBAYA_COCOA/$CBLIKE/$HGL/EE.yaml $COBAYA/$CBLIKE/$HGL/EE.yaml
    cp $COBAYA_COCOA/$CBLIKE/$HGL/TE.yaml $COBAYA/$CBLIKE/$HGL/TE.yaml

    cp $COBAYA_COCOA/$CBLIKE/$HGL/TT.yaml $COBAYA/$CBLIKE/$HGL/TT.yaml
    cp $COBAYA_COCOA/$CBLIKE/$HGL/TTTEEE.yaml $COBAYA/$CBLIKE/$HGL/TTTEEE.yaml

    cp $COBAYA_COCOA/$CBLIKE/$HGL/TT_lite.yaml $COBAYA/$CBLIKE/$HGL/TT_lite.yaml
    cp $COBAYA_COCOA/$CBLIKE/$HGL/TTTEEE_lite.yaml $COBAYA/$CBLIKE/$HGL/TTTEEE_lite.yaml

    #---------------------------------------------------------------------------
    # SPT3G Y1 -----------------------------------------------------------------
    #---------------------------------------------------------------------------
    cp -r $COBAYA_COCOA/$CBLIKE/SPT3G_Y1 $COBAYA/$CBLIKE/SPT3G_Y1

    #---------------------------------------------------------------------------
    # ACT DR6 LENSLIKE ---------------------------------------------------------
    #---------------------------------------------------------------------------
    cp -r $COBAYA_COCOA/$CBLIKE/act_dr6_lenslike $COBAYA/$CBLIKE/act_dr6_lenslike

    #---------------------------------------------------------------------------
    # FIX RENAMING IN CAMB -----------------------------------------------------
    #---------------------------------------------------------------------------
    cp $COBAYA_COCOA/$CBTH/camb/camb.yaml $COBAYA/$CBTH/camb/camb.yaml
  fi

  #-------------------------------------------------------------------------------
  # CAMSPEC ----------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  if [ -z "${IGNORE_CAMSPEC_INSTALLATION}" ]; then
    # WE REMOVE UNBINNED PLANCK CAMSPEC (TODO IN THE FUTURE: FIX THEIR YAML)
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec

    cp $COBAYA_COCOA/$CBLIKE/base_classes/InstallableLikelihood.patch $COBAYA/$CBLIKE/base_classes/
    cd COBAYA/$CBLIKE/base_classes/

    patch -u InstallableLikelihood.py -i InstallableLikelihood.patch
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (PATCH CAMSPEC) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi

    cd $ROOTDIR
  else
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec
    rm -rf $COBAYA/$CBLIKE/planck_2018_highl_CamSpec2021
  fi

  #-----------------------------------------------------------------------------
  # SO -------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${IGNORE_SO_INSTALLATION}" ]; then
    cp -r $COBAYA_COCOA/$CBLIKE/mflike $COBAYA/$CBLIKE/mflike
  fi

  #-----------------------------------------------------------------------------
  # INSTALL HILLIPOP -----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${IGNORE_LIPOP_INSTALLATION}" ]; then
    rm -rf $COBAYA/$CBLIKE/planck_2020_hillipop
    rm -rf $COBAYA/$CBLIKE/hillipop_tmp

    cd $COBAYA/$CBLIKE

    git clone https://github.com/planck-npipe/hillipop.git hillipop_tmp > ${OUT_UPT_CB_1} 2> ${OUT_UPT_CB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (HILLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi

    cd $COBAYA/$CBLIKE/hillipop_tmp

    git reset --hard $HILLIPOP_GIT_COMMIT > ${OUT_UPT_CB_1} 2> ${OUT_UPT_CB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (HILLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi
    mv planck_2020_hillipop/ $COBAYA/$CBLIKE

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$CBLIKE/planck_2020_hillipop/init.patch $COBAYA/$CBLIKE/planck_2020_hillipop/
    cd $COBAYA/$CBLIKE/planck_2020_hillipop/
    patch -u __init__.py -i init.patch

    cd $ROOTDIR

    rm -rf $COBAYA/$CBLIKE/hillipop_tmp
  fi

  #-----------------------------------------------------------------------------
  # INSTALL LOWLLIPOP ----------------------------------------------------------
  #-----------------------------------------------------------------------------
  if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
    rm -rf $COBAYA/$CBLIKE/planck_2020_lollipop
    rm -rf $COBAYA/$CBLIKE/lollipop_tmp

    cd $COBAYA/$CBLIKE

    git clone https://github.com/planck-npipe/lollipop.git lollipop_tmp
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (LOLLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi

    cd $COBAYA/$CBLIKE/lollipop_tmp
    git reset --hard $LOLLIPOP_GIT_COMMIT > ${OUT_UPT_CB_1} 2> ${OUT_UPT_CB_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"\t\t SETUP UPDATE COBAYA (LOLLIPOP) FAILED"'\033[0m'
      cd $ROOTDIR
      unset COBAYA
      unset COBAYA_COCOA
      unset CBLIKE
      unset CBTH
      unset HGL
      unset OUT_UPT_CB_1
      unset OUT_UPT_CB_2
      return 1
    fi

    mv planck_2020_lollipop/ $COBAYA/$CBLIKE

    # now patch the likelihood __init__ file
    cp $COBAYA_COCOA/$CBLIKE/planck_2020_lollipop/init.patch $COBAYA/$CBLIKE/planck_2020_lollipop/
    cd $COBAYA/$CBLIKE/planck_2020_lollipop/
    patch -u __init__.py -i init.patch

    cd $ROOTDIR
    rm -rf $COBAYA/$CBLIKE/lollipop_tmp
  fi

  #-----------------------------------------------------------------------------
  # UNSETKEYS ------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  cd $ROOTDIR
  unset COBAYA
  unset COBAYA_COCOA
  unset CBLIKE
  unset CBTH
  unset HGL
  unset OUT_UPT_CB_1
  unset OUT_UPT_CB_2
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------