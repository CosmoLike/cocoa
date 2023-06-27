#! /bin/bash

rm -r egfs_maude

mkdir egfs_maude
cp ../../minipmc/*.c egfs_maude
cp ../../minipmc/*.h egfs_maude
cp ../*.f90 egfs_maude
cp ../../clik_egfs.c egfs_maude
cp ../../clik_egfs.h egfs_maude
cp ../../clik_parametric.c egfs_maude
cp ../../clik_parametric.h egfs_maude
cp ../../python/clik/egfs.pyx egfs_maude
cp ../../python/clik/parametric.pyx egfs_maude
ln -s ../../python/clik/clik.parametric.pxd egfs_maude/parametric.pxd
cp ../../clik_parametric.c egfs_maude
cp ../../clik_parametric.h egfs_maude

cp setup.py egfs_maude
cp Makefile egfs_maude
mkdir egfs_maude/data
cp -r ../egfs_data/*.dat egfs_maude/data

tar -jcvf egfs_maude.tar.bz2 egfs_maude/
