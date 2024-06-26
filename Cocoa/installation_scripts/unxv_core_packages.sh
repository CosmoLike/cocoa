#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_CORE_PACKAGES}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;

  unset_env_vars () {
    unset -v TKEYS TFILES AL CCIL
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
    fail_script_msg "$(basename ${BASH_SOURCE[0]})" "${1}"
    unset_all || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  declare -a TKEYS=("IGNORE_CMAKE_INSTALLATION"
                    "IGNORE_DISTUTILS_INSTALLATION"
                    "IGNORE_DISTUTILS_INSTALLATION"
                    "IGNORE_OPENBLAS_INSTALLATION"
                    "IGNORE_FORTRAN_LAPACK_INSTALLATION"
                    "IGNORE_HDF5_INSTALLATION"
                    "IGNORE_C_CFITSIO_INSTALLATION"
                    "IGNORE_C_FFTW_INSTALLATION"
                    "IGNORE_C_GSL_INSTALLATION"
                    "IGNORE_CPP_SPDLOG_INSTALLATION"
                    "IGNORE_CPP_ARMA_INSTALLATION"
                    "IGNORE_CPP_BOOST_INSTALLATION"
                    "IGNORE_CPP_CARMA_INSTALLATION"
                    "IGNORE_ALL_PIP_INSTALLATION"
                    "IGNORE_ALL_PIP_INSTALLATION"
                    "IGNORE_ALL_PIP_INSTALLATION"
                   ) # T = TMP

  declare -a TFILES=("cmake"
                     "binutils"
                     "texinfo"
                     "OpenBLAS"
                     "lapack"
                     "hdf5"
                     "cfitsio"
                     "fftw"
                     "gsl"
                     "spdlog"
                     "armadillo"
                     "boost"
                     "carma"
                     "pip_cache"
                     "expat"
                     "ee2"
                    ) # T = TMP

  # ----------------------------------------------------------------------------

  ptop 'DECOMPRESSING FILES ON COCOA_INSTALLATION_LIBRARIES' || return 1

  if [ "${#TKEYS[@]}" -eq "${#TFILES[@]}" ]; then
    AL=${#TFILES[@]}
  else
    error ("logical error - array sizes do not match");
  fi
  
  for (( i=0; i<${AL:?}; i++ ));
  do

    cdfolder "${CCIL:?}" || return 1

    if [ -z "${TKEYS[$i]}" ]; then
      
      tar xf "${ROOTDIR:?}/installation_scripts/${TFILES[$i]}.xz" \ 
        >${OUT1:?} 2>${OUT2:?} || 
        { error "${EC25:?} (${TFILES[$i]}.xz)"; return 1; }
    
    fi
  
  done

  # ----------------------------------------------------------------------------

  unset_all || return 1;

  pbottom 'DECOMPRESSING FILES ON COCOA_INSTALLATION_LIBRARIES' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------