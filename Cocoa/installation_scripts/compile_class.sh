#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
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

  source $ROOTDIR/installation_scripts/clean_class.sh

  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OU_CL_1="/dev/null"
    export OU_CL_2="/dev/null"
  else
    export OU_CL_1="/dev/tty"
    export OU_CL_2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/class_public/
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  #Workaround Cocoa .gitignore entry on /include
  cp -r  ./include2 ./include > ${OU_CL_1} 2> ${OU_CL_2}

  CC=$C_COMPILER PYTHON=$PYTHON3 make all > ${OU_CL_1} 2> ${OU_CL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mMAKE ALL\e[0m"'\033[0m'
    cd $ROOTDIR
    unset OU_CL_1
    unset OU_CL_2
    return 1
  else
    echo -e '\033[0;32m'"\t\tCLASS RUN \e[3mMAKE ALL\e[0m\e\033[0;32m DONE"'\033[0m'
  fi
   
  cd ./python

  CC=$C_COMPILER $PYTHON3 setup.py build > ${OU_CL_1} 2> ${OU_CL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mPYTHON3 SETUP.PY BUILD\e[0m"'\033[0m'
    cd $ROOTDIR
    unset OU_CL_1
    unset OU_CL_2
    return 1
  else
    echo -e '\033[0;32m'"\t\t CLASS RUN \e[3mPYTHON3 SETUP.PY BUILD\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OU_CL_1
  unset OU_CL_2
  echo -e '\033[1;34m''\t\e[4mCOMPILING CLASS DONE''\033[0m'
fi
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------