#!/bin/bash

if [ -n "${THREAD_UNXZ}" ]; then
     echo 'DECOMPRESSION WILL HAPPEN IN PARALLEL'
fi

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
    echo 'DECOMPRESSING OPENBLAS'
    rm -rf ./OpenBLAS/
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
    echo 'DECOMPRESSING PIP CACHE'
    rm -rf ./pip_cache/
    
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf pip_cache.xz
        proc2=$!
    else
        tar xf pip_cache.xz &
        proc2=$!
    fi

    rm -rf ./expat241/
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf expat241.xz
        proc2X=$!
    else
        tar xf expat241.xz &
        proc2A=$!
    fi
else
  proc2=1
  proc2A=1
fi

if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    echo 'DECOMPRESSING CMAKE LIBRARY - THAT MIGHT TAKE A WHILE'
    rm -rf ./cmake-3.17.1/
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf cmake3171.xz
        proc3=$!
    else
        tar xf cmake3171.xz &
        proc3=$!
    fi
else
  proc3=1
fi

if [ -z "${IGNORE_CPP_INSTALLATION}" ]; then
    if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP ARMA LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./armadillo-9.860.2/
            tar xf armadillo.xz
            proc4=$!
        else
            rm -rf ./armadillo-9.860.2/
            tar xf armadillo.xz &
            proc4=$!
        fi
    fi

    if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP BOOST LIBRARY - THAT MIGHT TAKE A WHILE'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./boost_1_72_0/
            tar xf boost_1_72_0.xz
            proc5=$!
        else
            rm -rf ./boost_1_72_0/
            tar xf boost_1_72_0.xz &
            proc5=$!
        fi
    fi

    if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP SPDLOG LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./spdlog/
            tar xf spdlog.xz
            proc6=$!
        else
            rm -rf ./spdlog/
            tar xf spdlog.xz &
            proc6=$!
        fi
    fi

    if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
        echo 'DECOMPRESSING CPP CARMA LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./carma/
            tar xf carma.xz
            proc7=$!
        else
            rm -rf ./carma/
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
        echo 'DECOMPRESSING C FFTW LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./fftw-3.3.8/
            tar xf fftw338.xz
            proc8=$!
        else
            rm -rf ./fftw-3.3.8/
            tar xf fftw338.xz &
            proc=8$!
        fi
    fi
    if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
        echo 'DECOMPRESSING C CFITSIO LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./cfitsio-3.47/
            tar xf cfitsio347.xz
            proc9=$!
        else
            rm -rf ./cfitsio-3.47/
            tar xf cfitsio347.xz &
            proc9=$!
        fi
    fi
    if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
        echo 'DECOMPRESSING C GSL LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
            rm -rf ./gsl-2.7/
            tar xf gsl-2.7.xz
            proc10=$!
        else
            rm -rf ./gsl-2.7/
            tar xf gsl-2.7.xz &
            proc10=$!
        fi
    fi
else
  proc8=1
  proc9=1
  proc10=1
fi

if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
    if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
        echo 'DECOMPRESSING FORTRAN LAPACK LIBRARY'
        if [ -z "${THREAD_UNXZ}" ]; then
          rm -rf ./lapack-3.9.0/
          tar xf lapack390.xz
          proc11=$!
        else
          rm -rf ./lapack-3.9.0/
          tar xf lapack390.xz &
          proc11=$!
        fi
    else
        proc11=1
    fi
else
  proc11=1
fi

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
    echo 'DECOMPRESSING BINUTILS'
    rm -rf ./binutils-2.37/
    rm -rf ./texinfo-6.7/
    if [ -z "${THREAD_UNXZ}" ]; then
        tar xf texinfo-6.7.xz
        proc12=$!
        tar xf binutils-2.37.xz
        proc13=$!
    else
        tar xf texinfo-6.7.xz &
        proc12=$!
        tar xf binutils-2.37.xz &
        proc13=$!
    fi
else
  proc12=1
  proc13=1
fi

if [ -n "${THREAD_UNXZ}" ]; then
     echo 'DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH'
fi

wait "$proc1" "$proc2" "$proc2A" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9" "$proc10" "$proc11" "$proc12" "$proc13" 2>/dev/null > /dev/null