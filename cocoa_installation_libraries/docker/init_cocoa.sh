#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc "${CONDA_PREFIX}"/bin/gcc
ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-g++ "${CONDA_PREFIX}"/bin/g++
ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gfortran "${CONDA_PREFIX}"/bin/gfortran
ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ar "${CONDA_PREFIX}"/bin/gcc-ar
ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ranlib "${CONDA_PREFIX}"/bin/gcc-ranlib
ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-ld "${CONDA_PREFIX}"/bin/ld

git-lfs install
