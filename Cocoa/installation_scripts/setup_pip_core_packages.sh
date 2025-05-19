#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_PIP_CORE_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

  unset_env_vars () {
    unset -v URL_BASE URL FOLDER VER XZF CCIL CNAME PACKDIR PACKAGE_VERSION 
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder cpfile error wgetact gitact gitact1 gitact2
    unset -f unset_env_funcs wgetact1 wgetact2
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
   
  unset_env_vars || return 1;

  # ----------------------------------------------------------------------------
  # --------------------------- PIP CORE PACKAGES ------------------------------
  # ----------------------------------------------------------------------------  
  ptop "INSTALLING A FEW PYTHON CORE LIBRARIES VIA PIP" || return 1

  #PS: --force-reinstall - this helps CARMA to see numpy files
  #PS2: Need to include numpy in the same command to avoid numpy 2.0
  # mpi4py has a weird bug when installing from conda on a few machines 
  # (e.g., midway) no-cache-dir is important to fix this bug
  # https://github.com/mpi4py/mpi4py/issues/335
  if [ -n "${COCOA_FORCE_NUMPY_1_23}" ]; then
    COCOA_NUMPY_VERSION='1.23.5'

    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        "numpy==${COCOA_NUMPY_VERSION:?}" \
        'mpi4py==4.0.3' \
        'pyfftw==0.13.1' \
        'setuptools==80.3.1' \
      --no-cache-dir --prefer-binary \
      --prefix="${ROOTDIR:?}/.local" \
      --force-reinstall \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }
  else
    COCOA_NUMPY_VERSION='1.26.3'

    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        "numpy==${COCOA_NUMPY_VERSION:?}" \
        'mpi4py==4.0.3' \
        'setuptools==80.3.1' \
      --no-cache-dir --prefer-binary \
      --prefix="${ROOTDIR:?}/.local" \
      --force-reinstall \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }
  fi

  env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
      "numpy==${COCOA_NUMPY_VERSION:?}" \
      'mpi4py==4.0.3' \
      'notebook==7.4.2' \
      'ipyparallel==9.0.1' \
      'emcee== 3.1.6' \
      'sacc==1.0.2' \
      'george==0.4.4' \
    --no-cache-dir --prefer-binary --use-pep517 \
    --prefix="${ROOTDIR:?}/.local" \
    >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

  #env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
  #--upgrade setuptools --no-cache-dir >${OUT1:?} 2>${OUT2:?} \
  #|| { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

  pbottom "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1
  
  # ----------------------------------------------------------------------------
  # ----------------------------- PIP ML PACKAGES ------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
    ptop "PIP INSTALL MACHINE LEARNING GPU PACKAGES (takes a while O(5-10min)...)"

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow==2.17.0' \
        'tensorflow_probability==0.24.0' \
        'keras==3.9.2' \
        'keras-preprocessing==1.1.2' \
        'torch==2.6.0' \
        'torchvision==0.21.0' \
        'torchaudio==2.6.0' \
        'tensiometer==1.0.2' \
        'scikit-learn==1.6.1' \
        'jupyter==1.0.0' \
        'typing-extensions==4.13.2' \
        'mkdocs_material==9.6.13' \
        'mkdocstrings==0.29.1' \
        'pytest==8.3.5' \
        'tf-keras==2.17.0' \
      --no-cache-dir --prefer-binary \
      --extra-index-url "https://download.pytorch.org/whl/cu118" \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; } 

    # Without this code, jupyter breaks notebook
    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'notebook==7.4.2' \
        'ipyparallel==9.0.1' \
      --no-cache-dir --prefer-binary \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; }

    pbottom "PIP INSTALL MACHINE LEARNING GPU PACKAGES" || return 1
  
  fi

  # ----------------------------------------------------------------------------

  unset_all || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------