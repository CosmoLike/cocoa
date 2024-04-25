#!/bin/bash

if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
  export OUTPUT_UNXV_ALL_1="/dev/null"
  export OUTPUT_UNXV_ALL_2="/dev/null"
else
  export OUTPUT_UNXV_ALL_1="/dev/tty"
  export OUTPUT_UNXV_ALL_2="/dev/tty"
fi

if [ -z "${ROOTDIR}" ]; then
    echo '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
fi

sh $ROOTDIR/../cocoa_installation_libraries/clean_all
cd $ROOTDIR/../cocoa_installation_libraries/

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
    echo '\033[0;32m'"\t\t DECOMPRESSING OPENBLAS"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf OpenBLAS.xz
    else
        tar xf OpenBLAS.xz &
        proc1=$!
    fi
else
  proc1=1
fi

if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
    echo '\033[0;32m'"\t\tDECOMPRESSING PIP CACHE"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf pip_cache.xz
        proc2=$!
    else
        tar xf pip_cache.xz &
        proc2=$!
    fi
    if [ -z "${MINICONDA_INSTALLATION}" ]; then
      if [ -z "${THREAD_UNXZ}" ]; then
        tar xf expat.xz
        proc2A=$!
      else
        tar xf expat.xz &
        proc2A=$!
      fi
    else
      proc2A=1
    fi
else
  proc2=1
  proc2A=1
fi

if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    echo '\033[0;32m'"\t\tDECOMPRESSING CMAKE LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf cmake.xz
        proc3=$!
    else
        tar xf cmake.xz &
        proc3=$!
    fi
else
  proc3=1
fi

if [ -z "${IGNORE_CPP_INSTALLATION}" ]; then
    if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
        echo '\033[0;32m'"\t\t DECOMPRESSING CPP ARMA LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf armadillo.xz
            proc4=$!
        else
            tar xf armadillo.xz &
            proc4=$!
        fi
    else
      proc4=1
    fi

    if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
        echo '\033[0;32m'"\t\t DECOMPRESSING CPP BOOST LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf boost.xz
            proc5=$!
        else
            tar xf boost.xz &
            proc5=$!
        fi
    else
      proc5=1
    fi

    if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
        echo '\033[0;32m'"\t\t DECOMPRESSING CPP SPDLOG LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf spdlog.xz
            proc6=$!
        else
            tar xf spdlog.xz &
            proc6=$!
        fi
    else
      proc6=1
    fi

    if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
        echo '\033[0;32m'"\t\t DECOMPRESSING CPP CARMA LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf carma.xz
            proc7=$!
        else
            tar xf carma.xz &
            proc7=$!
        fi
    else
      proc7=1
    fi
else
  proc4=1
  proc5=1
  proc6=1
  proc7=1
fi

if [ -z "${IGNORE_C_INSTALLATION}" ]; then
    if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
        if [ -z "${THREAD_UNXZ}" ]; then
          echo '\033[0;32m'"\t\t DECOMPRESSING C FFTW LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
          tar xf fftw.xz
          proc8=$!
        else
          tar xf fftw.xz &
          proc8=$!
        fi
    else
      proc8=1
    fi
    if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
        echo '\033[0;32m'"\t\t DECOMPRESSING C CFITSIO LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf cfitsio.xz
            proc9=$!
        else
            tar xf cfitsio.xz &
            proc9=$!
        fi
    else
      proc9=1
    fi
    if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
        if [ -z "${THREAD_UNXZ}" ]; then
          echo '\033[0;32m'"\t\t DECOMPRESSING C GSL LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
          tar xf gsl.xz
          proc10=$!
        else
          tar xf gsl.xz &
          proc10=$!
        fi
    else
      proc10=1
    fi
else
  proc8=1
  proc9=1
  proc10=1
fi

if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    if [ -z "${THREAD_UNXZ}" ]; then
      echo '\033[0;32m'"\t\t DECOMPRESSING FORTRAN LAPACK LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
      tar xf lapack.xz
      proc11=$!
    else
      tar xf lapack.xz &
      proc11=$!
    fi
else
    proc11=1
fi

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
    if [ -z "${THREAD_UNXZ}" ]; then
        echo '\033[0;32m'"\t\t DECOMPRESSING BINUTILS"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
        tar xf texinfo.xz
        proc12=$!
        tar xf binutils.xz
        proc13=$!
    else
        tar xf texinfo.xz &
        proc12=$!
        tar xf binutils.xz &
        proc13=$!
    fi
else
  proc12=1
  proc13=1
fi

if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
    if [ -z "${THREAD_UNXZ}" ]; then
      echo '\033[0;32m'"\t\t DECOMPRESSING HDF5 LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
      tar xf hdf5.xz
      proc14=$!
    else
      tar xf hdf5.xz &
      proc14=$!
    fi
else
    proc14=1
fi

if [ -n "${THREAD_UNXZ}" ]; then
  echo '\033[0;32m'"\t\t DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
fi

wait "$proc1" "$proc2" "$proc2A" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9" "$proc10" "$proc11" "$proc12" "$proc13" "$proc14" "$proc15" 2>/dev/null > /dev/null

cd $ROOTDIR/
