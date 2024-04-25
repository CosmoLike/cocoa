if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING CMAKE LIBRARY"'\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/
    
    git clone  https://github.com/Kitware/CMake.git cmake-3.26.4
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CMAKE: COULD NOT RUN \e[3mGIT CLONE"'\033[0m'
        cd $ROOTDIR
        return 1
    fi
    
    cd ./cmake-3.26.4
    git checkout v3.26.4
    rm -rf ./.git/
    cd ../
    
    tar -cf - cmake-3.26.4/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > cmake.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CMAKE: COULD NOT COMPRESS \e[3mCMAKE FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    echo -e '\033[0;32m'"\t\tGETTING CMAKE LIBRARY DONE"'\033[0m'
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

fi

if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING CFITSIO LIBRARY"'\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/

    wget -q http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.0.0.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CFITSIO: COULD NOT RUN \e[3mWGET"'\033[0m'
        cd $ROOTDIR
        return 1
    fi
    tar zxvf cfitsio-4.0.0.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CFITSIO: COULD NOT RUN \e[3mTAR"'\033[0m'
        cd $ROOTDIR
        return 1
    fi
    rm -f cfitsio-4.0.0.tar.gz

    tar -cf - cfitsio-4.0.0/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > cfitsio.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CFITSIO: COULD NOT COMPRESS \e[3mCFITSIO FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    echo -e '\033[0;32m'"\t\tGETTING CFITSIO LIBRARY DONE"'\033[0m'
fi

if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY"'\033[0m'
    
    wget -q http://www.fftw.org/fftw-3.3.10.tar.gz 
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"FFTW: COULD NOT RUN \e[3mWGET"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    tar zxvf fftw-3.3.10.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"FFTW: COULD NOT RUN \e[3mTAR"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    rm -f fftw-3.3.10.tar.gz

    tar -cf - fftw-3.3.10/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > fftw.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"FFTW: COULD NOT COMPRESS \e[3mFFTW FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY DONE"'\033[0m'
fi

if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING GSL LIBRARY"'\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/

    wget -q http://ftp.wayne.edu/gnu/gsl/gsl-2.7.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"GSL: COULD NOT RUN \e[3mWGET"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    tar zxvf gsl-2.7.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"GSL: COULD NOT RUN \e[3mTAR"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    rm -f ROOTDIR/../cocoa_installation_libraries/gsl-2.7.tar.gz

    tar -cf - gsl-2.7/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > gsl.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"GSL: COULD NOT COMPRESS \e[3mGSL FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY DONE"'\033[0m'
fi

if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING ARMA LIBRARY"'\033[0m'
    
    cd $ROOTDIR/../cocoa_installation_libraries/

    wget -q https://sourceforge.net/projects/arma/files/armadillo-12.8.2.tar.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"ARMA: COULD NOT RUN \e[3mWGET"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    tar xf armadillo-12.8.2.tar.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"ARMA: COULD NOT RUN \e[3mTAR"'\033[0m'
        cd $ROOTDIR
        return 1
    fi
    
    rm -f ROOTDIR/../cocoa_installation_libraries/armadillo-12.8.2.tar.xz

    tar -cf - armadillo-12.8.2/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > armadillo.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"ARMA: COULD NOT COMPRESS \e[3mGSL FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    echo -e '\033[0;32m'"\t\tGETTING ARMA LIBRARY DONE"'\033[0m'
fi

if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING BOOST LIBRARY"'\033[0m'
    
    cd $ROOTDIR/../cocoa_installation_libraries/

    wget -q https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"BOOST: COULD NOT RUN \e[3mWGET"'\033[0m'
        cd $ROOTDIR
        
    fi
    
    tar xvzf boost_1_81_0.tar.gz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"BOOST: COULD NOT RUN \e[3mTAR"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    rm -f ROOTDIR/../cocoa_installation_libraries/boost_1_81_0.tar.gz
    
    tar -cf - boost_1_81_0/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > boost.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"BOOST: COULD NOT COMPRESS \e[3mGSL FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi

    echo -e '\033[0;32m'"\t\tGETTING BOOST LIBRARY DONE"'\033[0m'  
fi

if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING CARMA LIBRARY"'\033[0m'
    
    cd $ROOTDIR/../cocoa_installation_libraries/

    git clone https://github.com/RUrlus/carma.git carma_tmp
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CARMA: COULD NOT RUN \e[3mWGET"'\033[0m'
        cd $ROOTDIR
        
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/carma_tmp

    mv ./include ../

    cd ../

    mv ./include carma

    rm -rf $ROOTDIR/../cocoa_installation_libraries/carma_tmp

    tar -cf - carma/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > carma.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CARMA: COULD NOT COMPRESS \e[3mGSL FOLDER"'\033[0m'
        cd $ROOTDIR
        return 1
    fi  

    echo -e '\033[0;32m'"\t\tGETTING CARMA LIBRARY DONE"'\033[0m'  
fi
