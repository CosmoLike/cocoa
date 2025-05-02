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
 
  if [ -z "${IGNORE_PIP_CORE_PACKAGES}" ]; then

    ptop "INSTALLING ALL PYTHON CORE LIBRARIES VIA PIP" || return 1

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'alabaster==0.7.13' \
        'anytree==2.8.0' \
        'appdirs==1.4.4' \
        'astropy==5.2.2' \
        'babel==2.12.1' \
        'cachetools==5.3.1' \
        'certifi==2023.5.7' \
        'charset-normalizer==3.1.0' \
        'configparser==5.3.0' \
        'contourpy==1.1.0' \
        'corner==2.2.1' \
        'coverage==7.5.4' \
        'cycler==0.11.0' \
        'cython==3.0.10' \
        'decorator==5.1.1' \
        'deprecated==1.2.14' \
        'dill==0.3.6' \
        'distlib==0.3.8' \
        'docutils==0.20.1' \
        'filelock==3.13.4' \
        'fonttools==4.40.0' \
        'fuzzywuzzy==0.18.0' \
        'getdist==1.4.8' \
        'gpy==1.10.0' \
        'h5py==3.8.0' \
        'idna==3.4' \
        'imageio==2.31.1' \
        'imagesize==1.4.1' \
        'iminuit==2.25.2' \
        'importlib_metadata==6.6.0' \
        'importlib-resources==5.12.0' \
        'jax==0.4.12' \
        'Jinja2==3.1.2' \
        'joblib==1.4.0' \
        'johnnydep==1.20.2' \
        'kiwisolver==1.4.4' \
        'lazy_loader==0.2' \
        'lenstronomy==1.11.2' \
        'llvmlite==0.40.1' \
        'markupsafe==2.1.3' \
        'matplotlib==3.7.5' \
        'ml-dtypes==0.2.0' \
        'mpmath==1.3.0' \
        'multiprocess==0.70.14' \
        'networkx==3.1' \
        'numba==0.57.0' \
        'numpy==1.23.5'  \
        'numpydoc==1.5.0' \
        'opt-einsum==3.3.0' \
        'oauthlib==3.2.2' \
        'oyaml==1.0' \
        'packaging==23.1' \
        'pandas==2.0.3' \
        'paramz==0.9.5' \
        'pgen==0.2.1' \
        'Pillow==9.5.0' \
        'platformdirs==2.6.2' \
        'portalocker==2.7.0' \
        'protobuf==4.23.2' \
        'Py-bobyqa==1.4' \
        'pybind11==2.12.0' \
        'pyDOE2==1.3.0' \
        'pyerfa==2.0.0.3' \
        'pyfftw==0.13.1' \
        'pygments==2.17.2' \
        'pyparsing==3.0.9' \
        'python-dateutil==2.8.2' \
        'pytz==2023.3' \
        'PyWavelets==1.4.1' \
        'pyxdg==0.28' \
        'PyYAML==6.0' \
        'qp-prob==0.8.3' \
        'requests==2.31.0' \
        'sacc==0.10' \
        'schwimmbad==0.3.2' \
        'scikit-image==0.21.0' \
        'scikit-learn==1.2.2' \
        'scipy==1.10.1' \
        'setuptools==67.7.2' \
        'setuptools_scm==8.2.0' \
        'six==1.16.0' \
        'snowballstemmer==2.2.0' \
        'sphinx==7.1.2' \
        'sphinxcontrib-applehelp==1.0.4' \
        'sphinxcontrib-devhelp==1.0.2' \
        'sphinxcontrib-htmlhelp==2.0.1' \
        'sphinxcontrib-jsmath==1.0.1' \
        'sphinxcontrib-qthelp==1.0.3' \
        'sphinxcontrib-serializinghtml==1.1.5' \
        'structlog==23.1.0' \
        'sympy==1.12' \
        'syslibrary==0.1' \
        'tables-io==0.8.1' \
        'tabulate==0.9.0' \
        'threadpoolctl==3.1.0' \
        'tifffile==2023.4.12' \
        'tokenizers==0.13.3' \
        'toml==0.10.2' \
        'tomli==2.0.1' \
        'tqdm==4.65.0' \
        'typing_extensions==4.6.3' \
        'tzdata==2023.3' \
        'urllib3==1.26.16' \
        'virtualenv==20.17.1' \
        'wget==3.2' \
        'wheel==0.40.0' \
        'wimpy==0.6' \
        'wrapt==1.14.1' \
        'zipfile38==0.0.3' \
        'zipp==3.15.0' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(CORE-PACKAGES) ${EC13:?}"; return 1; }

    # mpi4py has a weird bug when installing from conda on a few machines 
    # (e.g., midway) no-cache-dir is important to fix this bug
    # https://github.com/mpi4py/mpi4py/issues/335
    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        'numpy==1.23.5' \
        'mpi4py==3.1.4' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      --force-reinstall \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        'notebook==7.1.1' \
        'ipyparallel==8.8.0' \
        'emcee==3.1.4' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    pbottom "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1

  else
    
    ptop "INSTALLING A FEW PYTHON CORE LIBRARIES VIA PIP" || return 1

    #PS: --force-reinstall - this helps CARMA to see numpy files
    #PS2: Need to include numpy in the same command to avoid numpy 2.0
    
    # mpi4py has a weird bug when installing from conda on a few machines 
    # (e.g., midway) no-cache-dir is important to fix this bug
    # https://github.com/mpi4py/mpi4py/issues/335
    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        'pip==25.0.1' \
        'numpy==1.23.5' \
        'mpi4py==3.1.4' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      --force-reinstall \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        'numpy==1.23.5' \
        'mpi4py==3.1.4' \
        'notebook==7.1.1' \
        'ipyparallel==8.8.0' \
        'emcee==3.1.4' \
        'sacc==0.10' \
        'setuptools_scm==8.2.0' \
      --no-cache-dir --use-pep517 \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
    --upgrade setuptools --no-cache-dir >${OUT1:?} 2>${OUT2:?} \
    || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    pbottom "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1
  
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------- PIP ML PACKAGES ------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_EMULATOR_CPU_PIP_PACKAGES}" ]; then
  
    if [ -z "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
      error "${EC28:?} (GPU AND CPU EMULATOR FLAGS)"; return 1;
    fi

    ptop "PIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES" || return 1

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow-cpu==2.13.1' \
        'tensorflow_probability==0.21' \
        'keras==2.13.1' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cpu' \
        'torchvision==0.14.1+cpu' \
        'torchaudio==0.13.1' \
        'tensiometer==0.1.2' \
        'regex==2024.11.6' \
        'mkdocs-material==9.5.15' \
        'mkdocs==1.5.3' \
        'mkdocstrings==0.24.2' \
        'mkdocs_autorefs==1.4.1' \
        'pytest==8.2.0' \
        'gdown==5.2.0' \
        'tqdm==4.66.4' \
        'scikit-learn==1.5.0' \
        'jupyter==1.0.0' \
      --extra-index-url "https://download.pytorch.org/whl/cpu" \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; }
  
    pbottom "PIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES"
  
  fi

  if [ -z "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
  
    if [ -z "${IGNORE_EMULATOR_CPU_PIP_PACKAGES}" ]; then
      error "${EC28:?} (GPU AND CPU EMULATOR FLAGS)"; return 1;
    fi

    ptop "PIP INSTALL MACHINE LEARNING GPU PACKAGES"

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow==2.13.1' \
        'tensorflow_probability==0.21' \
        'keras==2.13.1' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cu116' \
        'torchvision==0.14.1+cu116' \
        'torchaudio==0.13.1' \
        'tensiometer==0.1.2' \
        'regex==2024.11.6' \
        'mkdocs-material==9.5.15' \
        'mkdocs==1.5.3' \
        'mkdocstrings==0.24.2' \
        'mkdocs_autorefs==1.4.1' \
        'pytest==8.2.0' \
        'gdown==5.2.0' \
        'tqdm==4.66.4' \
        'scikit-learn==1.5.0' \
        'jupyter==1.0.0' \
      --extra-index-url "https://download.pytorch.org/whl/cu116" \
      --prefix="${ROOTDIR:?}/.local" --use-deprecated=legacy-resolver \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; } 

    pbottom "PIP INSTALL MACHINE LEARNING GPU PACKAGES" || return 1
  
  fi

  # ----------------------------------------------------------------------------

  unset_all || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------