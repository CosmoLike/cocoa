#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  # ----------------------------------------------------------------------------
  source $ROOTDIR/installation_scripts/clean_planck.sh
  # ----------------------------------------------------------------------------
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    pfail 'MAKE_NUM_THREADS'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset pfail
    unset CLIK_LAPALIBS
    unset CLIK_CFITSLIBS
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
  if [ -z "${DEBUG_PLANCK_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi
  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    export CLIK_CFITSLIBS=$ROOTDIR/.local/lib
  else
    export CLIK_CFITSLIBS=$GLOBAL_PACKAGES_LOCATION
  fi
  
  if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
    export CLIK_LAPALIBS=$ROOTDIR/.local
  else
    export CLIK_LAPALIBS=$GLOBAL_PACKAGES_LOCATION
  fi

  # ---------------------------------------------------------------------------
  ptop2 'COMPILING PLANCK'

  if [ -z "${USE_SPT_CLIK_PLANCK}" ]; then
    cd $ROOTDIR/external_modules/code/planck/code/plc_3.0/plc-3.1/
  else
    cd $ROOTDIR/external_modules/code/planck/code/spt_clik/
  fi
  
  FC=$FORTRAN_COMPILER CC=$C_COMPILER CXX=$CXX_COMPILER $PYTHON3 waf configure \
    --gcc --gfortran --cfitsio_islocal --prefix $ROOTDIR/.local \
    --lapack_prefix=${CLIK_LAPALIBS} --cfitsio_lib=${CLIK_CFITSLIBS} \
    --python=${PYTHON3} > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "WAF CONFIGURE"
    return 1
  fi
  
  $PYTHON3 waf install -v > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "WAF INSTALL"
    return 1
  fi

  unset_env_vars
  pbottom2 'COMPILING PLANCK'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------