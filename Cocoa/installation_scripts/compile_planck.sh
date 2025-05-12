#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_PLANCK_LIKELIHOOD_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi
  
  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v ECPCF CLIK_LAPACK_LIBS CLIK_CFITSIO_LIBS PRINTNAME
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error
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
    fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
    unset_all || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
    
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  # ---------------------------------------------------------------------------
  # Name to be printed on this shell script messages
  PRINTNAME="PLANCK LIKELIHOOD"

  ptop "COMPILING ${PRINTNAME:?}" || return 1

  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    CLIK_CFITSIO_LIBS="${ROOTDIR:?}/.local/lib"
  else
    if [ -z "${GLOBAL_PACKAGES_LOCATION}" ]; then
      pfail "GLOBAL_PACKAGES_LOCATION"; cdroot; return 1;
    fi
    CLIK_CFITSIO_LIBS="${GLOBAL_PACKAGES_LOCATION:?}"
  fi
  
  if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
    CLIK_LAPACK_LIBS="${ROOTDIR:?}/.local"
  else
    if [ -z "${GLOBAL_PACKAGES_LOCATION}" ]; then
      pfail "GLOBAL_PACKAGES_LOCATION"; cdroot; return 1;
    fi
    CLIK_LAPACK_LIBS="${GLOBAL_PACKAGES_LOCATION:?}"
  fi
  
  ECPCF="external_modules/code/planck/code"
  
  if [ -z "${USE_SPT_CLIK_PLANCK}" ]; then 
    PACKDIR="${ROOTDIR:?}/${ECPCF:?}/plc_3.0/plc-3.1"
  else
    PACKDIR="${ROOTDIR:?}/${ECPCF:?}/spt_clik"
  fi

  cdfolder "${PACKDIR:?}" || return 1;

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation
  rm -f  "${ROOTDIR:?}/.local"/bin/clik*
  rm -f  "${ROOTDIR:?}/.local"/lib/libclik_f90.so
  rm -f  "${ROOTDIR:?}/.local"/lib/libclik.so
  rm -rf "${ROOTDIR:?}/.local"/lib/python/site-packages/clik
  rm -rf "${ROOTDIR:?}/.local"/share/clik
  rm -f  "${ROOTDIR:?}/.local"/include/clik*
  rm -f  "${PACKDIR:?}/".lock-waf_*

  "${PYTHON3:?}" waf distclean >${OUT1:?} 2>${OUT2:?} || { error "${EC18:?}"; return 1; }

  # ---------------------------------------------------------------------------
  
  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" CXX="${CXX_COMPILER:?}" \
   "${PYTHON3:?}" waf configure --gcc --gfortran --cfitsio_islocal --prefix "${ROOTDIR:?}/.local" \
   --lapack_prefix="${CLIK_LAPACK_LIBS:?}" --cfitsio_lib="${CLIK_CFITSIO_LIBS:?}" \
   --python="${PYTHON3:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC5:?}"; return 1; })
  
  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   "${PYTHON3:?}" waf install -v >${OUT1:?} 2>${OUT2:?} || { error "${EC6:?}"; return 1; })
  
  pbottom "COMPILING ${PRINTNAME:?}" || return 1

  # ---------------------------------------------------------------------------

  unset_all || return 1

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------