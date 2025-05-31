#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
function addvar () {
    local tmp="${!1}" ;
    tmp="${tmp//:${2}:/:}" ;
    tmp="${tmp/#${2}:/}" ;
    tmp="${tmp/%:${2}/}" ;
    export $1="${2}:${tmp}" ;
}

if [ -n "${OMP_NUM_THREADS}" ]; then
  export OLD_OMP_NUM_THREADS=$OMP_NUM_THREADS
else
  export OLD_OMP_NUM_THREADS="x"
fi

if [ -n "${OMP_PROC_BIND}" ]; then
  export OLD_OMP_PROC_BIND=$OMP_PROC_BIND
else
  export OLD_OMP_PROC_BIND="x"
fi

if [ -n "${LD_LIBRARY_PATH}" ]; then
  export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
else
  export OLD_LD_LIBRARY_PATH="x"
fi

if [ -n "${PYTHONPATH}" ]; then
  export OLD_PYTHONPATH=$PYTHONPATH
else
  export OLD_PYTHONPATH="x"
fi

if [ -n "${PATH}" ]; then
  export OLD_PATH=$PATH
else
  export OLD_PATH="x"
fi

if [ -n "${C_INCLUDE_PATH}" ]; then
  export OLD_C_INCLUDE_PATH=$C_INCLUDE_PATH
else
  export OLD_C_INCLUDE_PATH="x"
fi

if [ -n "${CPLUS_INCLUDE_PATH}" ]; then
  export OLD_CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH
else
  export OLD_CPLUS_INCLUDE_PATH="x"
fi

if [ -n "${LDFLAGS}" ]; then
  export OLD_LDFLAGS=$LDFLAGS
else
  export OLD_LDFLAGS="x"
fi

if [ -n "${CPATH}" ]; then
  export OLD_CPATH=$CPATH
else
  export OLD_CPATH="x"
fi

if [ -n "${LD_RUN_PATH}" ]; then
  export OLD_LD_RUN_PATH=$LD_RUN_PATH
else
  export OLD_LD_RUN_PATH="x"
fi

if [ -n "${CMAKE_LIBRARY_PATH}" ]; then
  export OLD_CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH
else
  export OLD_CMAKE_LIBRARY_PATH="x"
fi

if [ -n "${CMAKE_INCLUDE_PATH}" ]; then
  export OLD_CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH
else
  export OLD_CMAKE_INCLUDE_PATH="x"
fi

if [ -n "${LIBRARY_PATH}" ]; then
  export OLD_LIBRARY_PATH=$LIBRARY_PATH
else
  export OLD_LIBRARY_PATH="x"
fi

if [ -n "${INCLUDEPATH}" ]; then
  export OLD_INCLUDEPATH=$INCLUDEPATH
else
  export OLD_INCLUDEPATH="x"
fi

if [ -n "${INCLUDE}" ]; then
  export OLD_INCLUDE=$INCLUDE
else
  export OLD_INCLUDE="x"
fi

if [ -n "${CPATH}" ]; then
  export OLD_CPATH=$CPATH
else
  export OLD_CPATH="x"
fi

if [ -n "${NUMPY_HEADERS}" ]; then
  export OLD_NUMPY_HEADERS=$NUMPY_HEADERS
else
  export OLD_NUMPY_HEADERS="x"
fi

if [ -n "${OBJC_INCLUDE_PATH}" ]; then
  export OLD_OBJC_INCLUDE_PATH=$OBJC_INCLUDE_PATH
else
  export OLD_OBJC_INCLUDE_PATH="x"
fi

if [ -n "${OBJC_PATH}" ]; then
  export OLD_OBJC_PATH=$OBJC_PATH
else
  export OLD_OBJC_PATH="x"
fi

if [ -n "${CFLAGS}" ]; then
  export OLD_CFLAGS=$CFLAGS
else
  export OLD_CFLAGS="x"
fi

if [ -n "${INCLUDE_PATH}" ]; then
  export OLD_INCLUDE_PATH=$INCLUDE_PATH
else
  export OLD_INCLUDE_PATH="x"
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------