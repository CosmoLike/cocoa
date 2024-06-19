#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
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
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'
    return 1
  fi
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    return 1
  fi
  unset_env_vars_spp () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset PIP_MNT
    unset unset_env_vars_spp
  }
  fail_spp () {
    export MSG="\033[0;31m (setup_pip_packages.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_spp
    unset MSG
    unset MSG2
  }
  if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export PIP_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export PIP_MNT=1
  fi

  ptop2 "INSTALLING PYTHON PACKAGES VIA PIP"

  # ----------------------------------------------------------------------------
  # --------------------------------- LIBEXPAT ---------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${MINICONDA_INSTALLATION}" ]; then
    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_EXPAT_DIR
    
    FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure --prefix=$ROOTDIR/.local \
      --enable-shared=yes --enable-static=yes > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "LIBEXPAT CONFIGURE"
      return 1
    fi

    make -j $PIP_MNT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "LIBEXPAT MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "LIBEXPAT MAKE INSTALL"
      return 1
    fi

    cp $ROOTDIR/.local/lib/libexpat.so.1 $ROOTDIR/.local/lib/libexpat.so.0
    if [ $? -ne 0 ]; then
      fail_spp "LIBEXPAT CP LIBEXPAT.SO.1"
      return 1
    fi

   cd $ROOTDIR
  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # --------------------------- PIP CORE PACKAGES ------------------------------
  # ----------------------------------------------------------------------------
  ptop "PIP INSTALL CORE PACKAGES"
  
  if [ -z "${MINICONDA_INSTALLATION}" ]; then
    
    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'alabaster==0.7.13' \
        'appdirs==1.4.4' \
        'anytree==2.8.0' \
        'astropy==5.2.2' \
        'babel==2.12.1' \
        'cachetools==5.3.1' \
        'certifi==2023.5.7' \
        'charset-normalizer==3.1.0' \
        'configparser==5.3.0' \
        'contourpy==1.1.0' \
        'corner==2.2.1' \
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
        'GetDist==1.4.3' \
        'GPy==1.10.0' \
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
        'MarkupSafe==2.1.3' \
        'matplotlib==3.7.5' \
        'mpi4py==3.1.4' \
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
        'PGen==0.2.1' \
        'Pillow==9.5.0' \
        'platformdirs==2.6.2' \
        'portalocker==2.7.0' \
        'protobuf==4.23.2' \
        'Py-BOBYQA==1.4' \
        'pybind11==2.12.0' \
        'pyDOE2==1.3.0' \
        'pyerfa==2.0.0.3' \
        'pygments==2.17.2' \
        'pyparsing==3.0.9' \
        'python-dateutil==2.8.2' \
        'pytz==2023.3' \
        'PyWavelets==1.4.1' \
        'pyxdg==0.28' \
        'PyYAML==6.0' \
        'qp-prob==0.8.3' \
        'requests==2.31.0' \
        'sacc==0.8.1' \
        'schwimmbad==0.3.2' \
        'scikit-image==0.21.0' \
        'scikit-learn==1.2.2' \
        'scipy==1.10.1' \
        'setuptools==67.7.2' \
        'setuptools-scm==7.1.0' \
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
      --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "PIP INSTALL (CORE PACKAGES)"
      return 1
    fi
  
  else
  
    #PS: --force-reinstall - this helps CARMA to see numpy files
    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'numpy==1.23.5' \
      --prefix=$ROOTDIR/.local \
      --force-reinstall > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "PIP INSTALL (CORE PACKAGES)"
      return 1
    fi
  
  fi

  env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
    $ROOTDIR/../cocoa_installation_libraries/pip_cache/fgspectra \
    --prefix=$ROOTDIR/.local \
    --no-index > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_spp "PIP INSTALL FGSPECTRA"
    return 1
  fi

  pbottom "PIP INSTALL CORE PACKAGES"
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # ----------------------------- PIP ML PACKAGES ------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_EMULATOR_CPU_PIP_PACKAGES}" ]; then
    ptop "PIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES"

    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'tensorflow-cpu==2.12.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cpu' \
        'torchvision==0.14.1+cpu' \
        'torchaudio==0.13.1' \
      --extra-index-url https://download.pytorch.org/whl/cpu
      --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "PIP INSTALL (MACHINE LEARNING CPU-ONLY PACKAGES)"
      return 1
    fi

    pbottom "PIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES"
  fi

  if [ -z "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
    ptop "PIP INSTALL MACHINE LEARNING GPU PACKAGES"

    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'tensorflow==2.12.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cu116' \
        'torchvision==0.14.1+cu116' \
        'torchaudio==0.13.1' \
      --extra-index-url https://download.pytorch.org/whl/cu116 \
      --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_spp "PIP INSTALL (MACHINE LEARNING GPU PACKAGES)"
      return 1
    fi

    pbottom "PIP INSTALL MACHINE LEARNING GPU PACKAGES"
  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_spp
  pbottom2 'INSTALLING PYTHON PACKAGES VIA PIP DONE'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------