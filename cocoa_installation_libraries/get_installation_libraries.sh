
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\tGETTING CMAKE LIBRARY"'\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/
    
    git clone  https://github.com/Kitware/CMake.git cmake-3.26.4
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CMAKE COULD NOT RUN \e[3mGIT CLONE"'\033[0m'
        cd $ROOTDIR
        return 1
    fi
    
    cd ./cmake-3.26.4
    git checkout v3.26.4
    rm -rf ./.git/
    cd ../
    
    tar -cf - cmake-3.26.4/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > cmake.xz
    if [ $? -neq 0 ];then
        echo -e '\033[0;31m'"CMAKE COULD NOT COMPRESS \e[3mCMAKE FOLDER"'\033[0m'
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


if [ -z "${IGNORE_C_INSTALLATION}" ]; then
    if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
        echo -e '\033[0;32m'"\t\tGETTING CFITSIO LIBRARY"'\033[0m'

        cd $ROOTDIR/../cocoa_installation_libraries/

        wget -q http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.0.0.tar.gz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"CFITSIO COULD NOT RUN \e[3mWGET"'\033[0m'
            cd $ROOTDIR
            return 1
        fi
        tar zxvf cfitsio-4.0.0.tar.gz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"CFITSIO COULD NOT RUN \e[3mTAR"'\033[0m'
            cd $ROOTDIR
            return 1
        fi
        rm -f cfitsio-4.0.0.tar.gz

        tar -cf - cfitsio-4.0.0/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > cfitsio.xz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"CMAKE COULD NOT COMPRESS \e[3mCFITSIO FOLDER"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        echo -e '\033[0;32m'"\t\tGETTING CFITSIO LIBRARY DONE"'\033[0m'
    fi

    if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
        echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY"'\033[0m'
        
        wget -q http://www.fftw.org/fftw-3.3.10.tar.gz 
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"FFTW COULD NOT RUN \e[3mWGET"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        tar zxvf fftw-3.3.10.tar.gz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"FFTW COULD NOT RUN \e[3mTAR"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        rm -f fftw-3.3.10.tar.gz

        tar -cf - fftw-3.3.10/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > fftw.xz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"FFTW COULD NOT COMPRESS \e[3mFFTW FOLDER"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY DONE"'\033[0m'
    fi

    if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
        echo -e '\033[0;32m'"\t\tGETTING GSL LIBRARY"'\033[0m'

        wget -q http://ftp.wayne.edu/gnu/gsl/gsl-2.7.tar.gz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"GSL COULD NOT RUN \e[3mWGET"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        tar zxvf gsl-2.7.tar.gz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"GSL COULD NOT RUN \e[3mTAR"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        rm -f gsl-2.7.tar.gz

        tar -cf - gsl-2.7/ | xz -k -1 --threads=$MAKE_NUM_THREADS -c - > gsl.xz
        if [ $? -neq 0 ];then
            echo -e '\033[0;31m'"GSL COULD NOT COMPRESS \e[3mGSL FOLDER"'\033[0m'
            cd $ROOTDIR
            return 1
        fi

        echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY DONE"'\033[0m'
    fi
fi


