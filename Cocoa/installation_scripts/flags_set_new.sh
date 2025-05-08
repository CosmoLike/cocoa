#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
export CFLAGS="${CFLAGS} -I${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/numpy/core/include"

export CMAKE_INCLUDE_PATH=${ROOTDIR:?}/.local/include/python${PYTHON_VERSION}m/:$CMAKE_INCLUDE_PATH
export CMAKE_INCLUDE_PATH=${ROOTDIR:?}/.local/include/:$CMAKE_INCLUDE_PATH    

export CMAKE_LIBRARY_PATH=${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages:$CMAKE_LIBRARY_PATH
export CMAKE_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$CMAKE_LIBRARY_PATH

export C_INCLUDE_PATH=$C_INCLUDE_PATH:${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/numpy/core/include

export CPATH=${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/numpy/core/include:$CPATH
export CPATH=${ROOTDIR:?}/.local/include/:$CPATH

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/numpy/core/include

export INCLUDE_PATH=${ROOTDIR:?}/.local/include/:$INCLUDE_PATH

export INCLUDEPATH=${ROOTDIR:?}/.local/include/:$INCLUDEPATH

export INCLUDE=${ROOTDIR:?}/.local/include/:$INCLUDE

export LDFLAGS="${LDFLAGS} -L${ROOTDIR:?}/.local/lib"
export LDFLAGS="${LDFLAGS} -L${ROOTDIR:?}/external_modules/code/${POLY_NAME}/lib  -Wl,-rpath-link,${ROOTDIR:?}/external_modules/code/${POLY_NAME}/lib"

export LD_RUN_PATH=${ROOTDIR:?}/.local/lib:$LD_RUN_PATH

export LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LIBRARY_PATH

export NUMPY_HEADERS=${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/numpy/core/include

export OBJC_INCLUDE_PATH=${ROOTDIR:?}/.local/:OBJC_INCLUDE_PATH

export OBJC_PATH=${ROOTDIR:?}/.local/:$OBJC_PATH

export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages/numpy/core/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOTDIR:?}/external_modules/code/${POLY_NAME}/lib

export PATH=${ROOTDIR:?}/.local/bin:$PATH

export PYTHONPATH=${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages:$PYTHONPATH
export PYTHONPATH=${ROOTDIR:?}/.local/lib:$PYTHONPATH
export PYTHONPATH=${ROOTDIR:?}/.local/lib/python/site-packages:$PYTHONPATH

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------