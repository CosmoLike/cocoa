#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${OLD_OMP_NUM_THREADS}" ]; then
  if [ "${OLD_OMP_NUM_THREADS}" != "x" ]; then
    export OMP_NUM_THREADS=$OLD_OMP_NUM_THREADS
  fi
  unset OLD_OMP_NUM_THREADS
fi

if [ -n "${OLD_OMP_PROC_BIND}" ]; then
  if [ "${OLD_OMP_PROC_BIND}" != "x" ]; then
    export OMP_PROC_BIND=$OLD_OMP_PROC_BIND
  fi
  unset OLD_OMP_PROC_BIND
fi

if [ -n "${OLD_LD_LIBRARY_PATH}" ]; then
  if [ "${OLD_LD_LIBRARY_PATH}" != "x" ]; then
    export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
  fi
  unset OLD_LD_LIBRARY_PATH
fi

if [ -n "${OLD_PYTHONPATH}" ]; then
  if [ "${OLD_PYTHONPATH}" != "x" ]; then
    export PYTHONPATH=$OLD_PYTHONPATH 
  else
    unset PYTHONPATH
  fi
  unset OLD_PYTHONPATH
fi

if [ -n "${OLD_PATH}" ]; then
  if [ "${OLD_PATH}" != "x" ]; then
    export PATH=$OLD_PATH
  else
    unset PATH
  fi
  unset OLD_PATH
fi

if [ -n "${OLD_C_INCLUDE_PATH}" ]; then
  if [ "${OLD_C_INCLUDE_PATH}" != "x" ]; then
    export C_INCLUDE_PATH=$OLD_C_INCLUDE_PATH 
  else
    unset C_INCLUDE_PATH
  fi
  unset OLD_C_INCLUDE_PATH
fi

if [ -n "${OLD_CPLUS_INCLUDE_PATH}" ]; then
  if [ "${OLD_CPLUS_INCLUDE_PATH}" != "x" ]; then
    export CPLUS_INCLUDE_PATH=$OLD_CPLUS_INCLUDE_PATH
  else
    unset CPLUS_INCLUDE_PATH
  fi
  unset OLD_CPLUS_INCLUDE_PATH
fi

if [ -n "${OLD_LDFLAGS}" ]; then
  if [ "${OLD_LDFLAGS}" != "x" ]; then
    export LDFLAGS=$OLD_LDFLAGS 
  else
    unset LDFLAGS
  fi
  unset OLD_LDFLAGS
fi

if [ -n "${OLD_CPATH}" ]; then
  if [ "${OLD_CPATH}" != "x" ]; then
    export CPATH=$OLD_CPATH
  else
    unset CPATH
  fi
  unset OLD_CPATH
fi

if [ -n "${OLD_LD_RUN_PATH}" ]; then
  if [ "${OLD_LD_RUN_PATH}" != "x" ]; then
    export LD_RUN_PATH=$OLD_LD_RUN_PATH 
  else
    unset LD_RUN_PATH
  fi
  unset OLD_LD_RUN_PATH
fi

if [ -n "${OLD_CMAKE_LIBRARY_PATH}" ]; then
  if [ "${OLD_CMAKE_LIBRARY_PATH}" != "x" ]; then
    export CMAKE_LIBRARY_PATH=$OLD_CMAKE_LIBRARY_PATH
  else
    unset CMAKE_LIBRARY_PATH
  fi
  unset OLD_CMAKE_LIBRARY_PATH
fi

if [ -n "${OLD_CMAKE_INCLUDE_PATH}" ]; then
  if [ "${OLD_CMAKE_INCLUDE_PATH}" != "x" ]; then
    export CMAKE_INCLUDE_PATH=$OLD_CMAKE_INCLUDE_PATH
  else
    unset CMAKE_INCLUDE_PATH
  fi
  unset OLD_CMAKE_INCLUDE_PATH
fi

if [ -n "${OLD_LIBRARY_PATH}" ]; then
  if [ "${OLD_LIBRARY_PATH}" != "x" ]; then
    export LIBRARY_PATH=$OLD_LIBRARY_PATH
  else
    unset LIBRARY_PATH
  fi
  unset OLD_LIBRARY_PATH
fi

if [ -n "${OLD_INCLUDEPATH}" ]; then
  if [ "${OLD_INCLUDEPATH}" != "x" ]; then
    export INCLUDEPATH=$OLD_INCLUDEPATH
  else
    unset INCLUDEPATH
  fi
  unset OLD_INCLUDEPATH
fi

if [ -n "${OLD_INCLUDE}" ]; then
  if [ "${OLD_INCLUDE}" != "x" ]; then
    export INCLUDE=$OLD_INCLUDE
  else
    unset INCLUDE
  fi
  unset OLD_INCLUDE
fi

if [ -n "${OLD_CPATH}" ]; then
  if [ "${OLD_CPATH}" != "x" ]; then
    export $CPATH=OLD_CPATH
  else
    unset CPATH
  fi
  unset CPATH
fi

if [ -n "${OLD_OBJC_INCLUDE_PATH}" ]; then
  if [ "${OLD_OBJC_INCLUDE_PATH}" != "x" ]; then
    export OBJC_INCLUDE_PATH=$OLD_OBJC_INCLUDE_PATH
  else
    unset OBJC_INCLUDE_PATH
  fi
  unset OLD_OBJC_INCLUDE_PATH
fi

if [ -n "${OLD_OBJC_PATH}" ]; then
  if [ "${OLD_OBJC_PATH}" != "x" ]; then
    export OBJC_PATH=$OLD_OBJC_PATH 
  else
    unset OBJC_PATH
  fi
  unset OLD_OBJC_PATH
fi

if [ -n "${OLD_NUMPY_HEADERS}" ]; then
  if [ "${OLD_NUMPY_HEADERS}" != "x" ]; then
    export NUMPY_HEADERS=$OLD_NUMPY_HEADERS
  else
    unset NUMPY_HEADERS
  fi
  unset OLD_NUMPY_HEADERS
fi

if [ -n "${OLD_CFLAGS}" ]; then
  if [ "${OLD_CFLAGS}" != "x" ]; then
    export CFLAGS=$OLD_CFLAGS
  else
    unset CFLAGS
  fi
  unset OLD_CFLAGS
fi

if [ -n "${OLD_INCLUDE_PATH}" ]; then
  if [ "${OLD_INCLUDE_PATH}" != "x" ]; then
    export INCLUDE_PATH=$OLD_INCLUDE_PATH
  else
    unset INCLUDE_PATH
  fi
  unset OLD_INCLUDE_PATH
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------