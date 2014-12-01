#!/bin/bash
#
# 
#
echo "Get FFTW ..."
wget http://www.fftw.org/fftw-3.3.4.tar.gz
echo "Get HDF5 ..."
wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
echo "Get BOOST ..."
wget http://sourceforge.net/projects/boost/files/boost/1.47.0/boost_1_47_0.tar.gz
echo "Get UG ..."
wget http://conan.iwr.uni-heidelberg.de/download/ug-3.11.0.tar.gz
echo "Get metis (for ALUGrid) ..."
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
echo "Get ALUGrid version 1.52 ..."
wget http://dune.mathematik.uni-freiburg.de/downloads/ALUGrid-1.52.tar.gz
