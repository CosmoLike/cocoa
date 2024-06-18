#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_clean_planck () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_clean_planck
  }
  fail_clean_planck () {
    export MSG="\033[0;31m (clean_planck.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_clean_planck
    unset MSG
    unset MSG2
    unset fail_clean_planck
  }
  if [ -z "${DEBUG_PLANCK_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  ptop 'CLEANING PLANCK LIKELIHOOD'

  if [ -z "${USE_SPT_CLIK_PLANCK}" ]; then
    cd $ROOTDIR/external_modules/code/planck/code/plc_3.0/plc-3.1/
  else
    cd $ROOTDIR/external_modules/code/planck/code/spt_clik/
  fi

  rm -f $ROOTDIR/.local/bin/clik*

  rm -f $ROOTDIR/.local/lib/libclik_f90.so

  rm -f $ROOTDIR/.local/lib/libclik.so

  rm -rf $ROOTDIR/.local/lib/python/site-packages/clik

  rm -rf $ROOTDIR/.local/share/clik

  rm -f $ROOTDIR/.local/include/clik*

  rm -f .lock-waf_*

  $PYTHON3 waf distclean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_clean_planck "PYTHON WAF DISTCLEAN"
    return 1
  fi

  unset_env_vars_clean_planck
  pbottom 'CLEANING PLANCK LIKELIHOOD'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------