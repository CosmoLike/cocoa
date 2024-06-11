#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
echo -e '\033[1;44m''INSTALLING PYTHON PACKAGES VIA PIP''\033[0m'

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
  export OUTPUT_PIP_1="/dev/null"
  export OUTPUT_PIP_2="/dev/null"
  export PIP_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
else
  export OUTPUT_PIP_1="/dev/tty"
  export OUTPUT_PIP_2="/dev/tty"
  export PIP_MAKE_NUM_THREADS=1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# --------------------------------- LIBEXPAT ---------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  if [ -z "${MINICONDA_INSTALLATION}" ]; then
    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_EXPAT_DIR
    
    FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure --prefix=$ROOTDIR/.local \
      --enable-shared=yes --enable-static=yes > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}

    make -j $PIP_MAKE_NUM_THREADS > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"LIBEXPAT RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"LIBEXPAT COULD NOT RUN \e[3mMAKE"'\033[0m'
      return 1
    fi

    make install > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"LIBEXPAT RUN \e[3m MAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"LIBEXPAT COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
      return 1
    fi

   cp $ROOTDIR/.local/lib/libexpat.so.1 $ROOTDIR/.local/lib/libexpat.so.0
   cd $ROOTDIR
  fi
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# -------------------------- PIP BASIC PACKAGES ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  echo -e '\033[1;34m''\tPIP INSTALL BASIC PACKAGES''\033[0m'

  if [ -z "${MINICONDA_INSTALLATION}" ]; then
    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'cython==3.0.10' \
        'mpi4py==3.1.4' \
        'scipy==1.10.1' \
        'pandas==2.0.3' \
        'numpy==1.23.5' \
        'matplotlib==3.7.5' \
      --prefix=$ROOTDIR/.local > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}  
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL (BASIC PACKAGES)"'\033[0m'
      cd $ROOTDIR
      return 1
    else
      echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL (BASIC PACKAGES)\e[0m\e\033[0;32m DONE"'\033[0m'
    fi
  else
    #PS: --force-reinstall - this helps CARMA to see numpy files
    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'numpy==1.23.5' \
      --prefix=$ROOTDIR/.local \
      --force-reinstall > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL (BASIC PACKAGES)"'\033[0m'
      cd $ROOTDIR
      return 1
    else
      echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL (BASIC PACKAGES)\e[0m\e\033[0;32m DONE"'\033[0m'
    fi
  fi

  echo -e '\033[1;34m''\t\e[4mPIP INSTALL (BASIC PACKAGES) DONE''\033[0m'
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# --------------------------- PIP CORE PACKAGES ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  echo -e '\033[1;34m''\tPIP INSTALL CORE PACKAGES''\033[0m'
  
  if [ -z "${MINICONDA_INSTALLATION}" ]; then
    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'appdirs==1.4.4' \
        'anytree==2.8.0' \
        'astropy==5.2.2' \
        'Babel==2.12.1' \
        'certifi==2023.5.7' \
        'contourpy==1.1.0' \
        'cachetools==5.3.1' \
        'cycler==0.11.0' \
        'charset-normalizer==3.1.0' \
        'configparser==5.3.0' \
        'corner==2.2.1' \
        'dill==0.3.6' \
        'decorator==5.1.1' \
        'docutils==0.20.1' \
        'fonttools==4.40.0' \
        'fuzzywuzzy==0.18.0' \
        'h5py==3.8.0' \
        'idna==3.4' \
        'imageio==2.31.1' \
        'imagesize==1.4.1' \
        'importlib-resources==5.12.0' \
        'importlib_metadata==6.6.0' \
        'iminuit==2.21.3' \
        'GetDist==1.4.3' \
        'GPy==1.10.0' \
        'Jinja2==3.1.2' \
        'jax==0.4.12' \
        'johnnydep==1.20.2' \
        'importlib-metadata==6.6.0' \
        'kiwisolver==1.4.4' \
        'lazy_loader==0.2' \
        'llvmlite==0.40.1' \
        'lenstronomy==1.11.2' \
        'MarkupSafe==2.1.3' \
        'mpmath==1.3.0' \
        'multiprocess==0.70.14' \
        'numba==0.57.0' \
        'numpydoc==1.5.0' \
        'networkx==3.1' \
        'oauthlib==3.2.2' \
        'oyaml==1.0' \
        'pyDOE2==1.3.0' \
        'PyWavelets==1.4.1' \
        'packaging==23.1' \
        'pyparsing==3.0.9' \
        'python-dateutil==2.8.2' \
        'pytz==2023.3' \
        'Pillow==9.5.0' \
        'PyYAML==6.0' \
        'pyerfa==2.0.0.3' \
        'protobuf==4.23.2' \
        'Py-BOBYQA==1.4' \
        'pybind11==2.10.4' \
        'PGen==0.2.1' \
        'portalocker==2.7.0' \
        'paramz==0.9.5' \
        'qp-prob==0.8.3' \
        'requests==2.31.0' \
        'sympy==1.12' \
        'sacc==0.8.1' \
        'six==1.16.0' \
        'setuptools==67.7.2' \
        'setuptools-scm==7.1.0' \
        'scikit-image==0.21.0' \
        'scikit-learn==1.2.2' \
        'syslibrary==0.1' \
        'schwimmbad==0.3.2' \
        'tables-io==0.8.1' \
        'tifffile==2023.4.12' \
        'tzdata==2023.3' \
        'threadpoolctl==3.1.0' \
        'typing_extensions==4.6.3' \
        'tokenizers==0.13.3' \
        'tqdm==4.65.0' \
        'tomli==2.0.1' \
        'toml==0.10.2' \
        'tabulate==0.9.0' \
        'urllib3==1.26.16' \
        'wrapt==1.14.1' \
        'wheel==0.40.0' \
        'wget==3.2' \
        'wimpy==0.6' \
        'zipp==3.15.0' \
        'zipfile38==0.0.3' \
        'structlog==23.1.0' \
      --prefix=$ROOTDIR/.local > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL (CORE PACKAGES)"'\033[0m'
      cd $ROOTDIR
      return 1
    else
      echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL (CORE PACKAGES)\e[0m\e\033[0;32m DONE"'\033[0m'
    fi
  fi

  env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
    $ROOTDIR/../cocoa_installation_libraries/pip_cache/fgspectra \
    --prefix=$ROOTDIR/.local \
    --no-index > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL FGSPECTRA"'\033[0m'
    cd $ROOTDIR
    return 1
  else
    echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL FGSPECTRA\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  echo -e '\033[1;34m''\t\e[4mPIP INSTALL (CORE PACKAGES) DONE''\033[0m'
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------- PIP ML PACKAGES ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  if [ -z "${IGNORE_EMULATOR_CPU_PIP_PACKAGES}" ]; then
    echo -e '\033[1;34m''\tPIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES''\033[0m'

    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'tensorflow-cpu==2.12.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cpu' \
        'torchvision==0.14.1+cpu' \
        'torchaudio==0.13.1' \
      --extra-index-url https://download.pytorch.org/whl/cpu
      --prefix=$ROOTDIR/.local > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL (MACHINE LEARNING CPU-ONLY PACKAGES)"'\033[0m'
      cd $ROOTDIR
      return 1
    else
      echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL (MACHINE LEARNING CPU-ONLY PACKAGES)\e[0m\e\033[0;32m DONE"'\033[0m'
    fi

    echo -e '\033[1;34m''\t\e[4mPIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES DONE''\033[0m'
  fi
  if [ -z "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
    echo -e '\033[1;34m''\tPIP INSTALL MACHINE LEARNING GPU PACKAGES''\033[0m'

    env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install \
        'tensorflow==2.12.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cu116' \
        'torchvision==0.14.1+cu116' \
        'torchaudio==0.13.1' \
      --extra-index-url https://download.pytorch.org/whl/cu116 \
      --prefix=$ROOTDIR/.local > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
    if [ $? -ne 0 ]; then
      echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL (MACHINE LEARNING GPU PACKAGES)"'\033[0m'
      cd $ROOTDIR
      return 1
    else
      echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL (MACHINE LEARNING GPU PACKAGES)\e[0m\e\033[0;32m DONE"'\033[0m'
    fi

    echo -e '\033[1;34m''\t\e[4mPIP INSTALL MACHINE LEARNING GPU PACKAGES DONE''\033[0m'
  fi

fi

echo -e '\033[1;44m''\e[4mINSTALLING PYTHON PACKAGES VIA PIP DONE''\033[0m'
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------