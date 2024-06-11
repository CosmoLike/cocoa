# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
echo -e '\033[1;34m''\tCOMPILING CLASS''\033[0m'

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
if [ -z "${FORTRAN_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE FORTRAN_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi
if [ -z "${PYTHON3}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE PYTHON3 IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi
if [ -z "${MAKE_NUM_THREADS}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE_NUM_THREADS IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi

if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then

  cd $ROOTDIR/external_modules/code/class_public/
  
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUTPUT_CLASS_1="/dev/null"
    export OUTPUT_CLASS_2="/dev/null"
  else
    export OUTPUT_CLASS_1="/dev/tty"
    export OUTPUT_CLASS_2="/dev/tty"
  fi

  #Workaround Cocoa .gitignore entry on /include
  rm -rf ./include > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  cp -r  ./include2 ./include > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}

  make clean > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mMAKE CLEAN\e[0m"'\033[0m'
  else
    echo -e '\033[0;32m'"\t\tCLASS RUN \e[3mMAKE CLEAN RUN\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  rm -f class

  cd ./python

  $PYTHON3 setup.py clean > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mPYTHON3 SETUP.PY CLEAN\e[0m"'\033[0m'
    cd $ROOTDIR
    return 1
  else
    echo -e '\033[0;32m'"\t\tCLASS RUN \e[3mPYTHON3 SETUP.PY CLEAN\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  rm -rf ./build/

  cd ../

  CC=$C_COMPILER PYTHON=$PYTHON3 make all > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mMAKE ALL\e[0m"'\033[0m'
    cd $ROOTDIR
    return 1
  else
    echo -e '\033[0;32m'"\t\tCLASS RUN \e[3mMAKE ALL\e[0m\e\033[0;32m DONE"'\033[0m'
  fi
   
  cd ./python

  CC=$C_COMPILER $PYTHON3 setup.py build > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mPYTHON3 SETUP.PY BUILD\e[0m"'\033[0m'
    cd $ROOTDIR
    return 1
  else
    echo -e '\033[0;32m'"\t\tCLASS RUN \e[3mPYTHON3 SETUP.PY BUILD\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
fi

echo -e '\033[1;34m''\t\e[4mCOMPILING CLASS DONE''\033[0m'
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------