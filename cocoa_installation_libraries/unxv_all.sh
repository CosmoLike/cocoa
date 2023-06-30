#!/bin/bash

if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED' >&2
    return 1
fi

source $ROOTDIR/../cocoa_installation_libraries/clean_all
cd $ROOTDIR/../cocoa_installation_libraries/

echo 'DECOMPRESSION OF INSTALLATION LIBRARIES - THAT MIGHT TAKE A WHILE'

if [ -n "${THREAD_UNXZ}" ]; then
     echo 'DECOMPRESSION WILL HAPPEN IN PARALLEL'
fi

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
    echo 'DECOMPRESSING OPENBLAS - THAT MIGHT TAKE A WHILE'
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
    echo 'DECOMPRESSING PIP CACHE - THAT MIGHT TAKE A WHILE'
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
    echo 'DECOMPRESSING CMAKE LIBRARY - THAT MIGHT TAKE A WHILE'
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
        echo 'DECOMPRESSING CPP ARMA LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf armadillo.xz
            proc4=$!
        else
            tar xf armadillo.xz &
            proc4=$!
        fi
    fi

    if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP BOOST LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf boost.xz
            proc5=$!
        else
            tar xf boost.xz &
            proc5=$!
        fi
    fi

    if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP SPDLOG LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf spdlog.xz
            proc6=$!
        else
            tar xf spdlog.xz &
            proc6=$!
        fi
    fi

    if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP CARMA LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf carma.xz
            proc7=$!
        else
            tar xf carma.xz &
            proc7=$!
        fi
    fi
else
  proc4=1
  proc5=1
  proc6=1
  proc7=1
fi

if [ -z "${IGNORE_C_INSTALLATION}" ]; then
    if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
        echo 'DECOMPRESSING C FFTW LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf fftw.xz
            proc8=$!
        else
            tar xf fftw.xz &
            proc8=$!
        fi
    fi
    if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
        echo 'DECOMPRESSING C CFITSIO LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf cfitsio.xz
            proc9=$!
        else
            tar xf cfitsio.xz &
            proc9=$!
        fi
    fi
    if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
        echo 'DECOMPRESSING C GSL LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            tar xf gsl.xz
            proc10=$!
        else
            tar xf gsl.xz &
            proc10=$!
        fi
    fi
else
  proc8=1
  proc9=1
  proc10=1
fi

if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    echo 'DECOMPRESSING FORTRAN LAPACK LIBRARY - THAT MIGHT TAKE A WHILE'
    if [ -z "${THREAD_UNXZ}" ]; then
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
    echo 'DECOMPRESSING BINUTILS - THAT MIGHT TAKE A WHILE'
    if [ -z "${THREAD_UNXZ}" ]; then
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
    echo 'DECOMPRESSING HDF5 LIBRARY - THAT MIGHT TAKE A WHILE'
    if [ -z "${THREAD_UNXZ}" ]; then
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
     echo 'DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH'
fi

wait "$proc1" "$proc2" "$proc2A" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9" "$proc10" "$proc11" "$proc12" "$proc13" "$proc14" 2>/dev/null > /dev/null

cd $ROOTDIR/

echo 'DECOMPRESSION OF INSTALLATION LIBRARIES - DONE'