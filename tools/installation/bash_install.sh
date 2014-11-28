#!/bin/bash
#
# Script to install all extra packages needed for DUNE-GESIS
# Author: A. Ngo

#
# Usage: 
# ./bash_install <package-name-version>.tar.gz
#
# Example:
# ./bash_install fftw-3.3.4.tar.gz
#
# would unpack the contents into a new directory "fftw-3.3.4/",
# run the configuration and the build before installing the binaries into
# the destination directory specified by the variables
# $installdir and $softwaredir (see below if you would like to adjust these paths)
# The local directory "fftw-3.3.4/" will be removed after the job is done.
#
#
# Dependencies:
#
#   mpic++ should point mpic++.openmpi
#   mpicc should point mpicc.openmpi
#
#   
#
# Building software libraries on your own and installing them 
# into your home directory is generally not a bad idea 
# if you have enough harddisk space available.
# This makes you independent of the availibility of software packages or modules 
# on remote systems such as compute servers or HPC clusters 
# where you ususally do not have administrative privileges!
# This installation path can be customized according to your environment:
softwaredir=$HOME/Software
#
#
#
create_directory()
{
    if [ ! -d $1 ]
    then
	mkdir $1
    else
	echo "$1 exists."
    fi
}


make_install()
{
    packagename=$1
    installdir=$softwaredir/$packagename
    metis_installdir=$softwaredir/metis-4.0.3 # metis is required only for dune-alugrid (3D parallel adaptive code)

    echo "processing package '$packagename'"
    if [ ! -f "$packagename.tar.gz" ]
    then
	echo "error: $packagename.tar.gz not found!"
	exit
    fi
    #
    if [[ "$packagename" == *"boost"* ]]; # means: $packagename contains boost
    then
        echo "Special unpacking BOOST to destination directory + symbolic link to itself named include/."
        echo "This is required ONLY due to the dune-common m4 check for boost."
        echo "processing package '$packagename.tar.gz'"
        create_directory "$installdir"
        tar xzvf $packagename.tar.gz -C $installdir/.
        cd "$installdir" 
        ln -s . include
        echo "Special unpacking BOOST to destination directory + symbolic link to itself named include/ done. Finished. No build needed."
        echo "Your special Boost version is now in $installdir."
        exit
    else
        echo "Default unpacking."
        tar xzvf $packagename.tar.gz
        echo "Default unpacking done."
    fi
    #
    if [ "$packagename" = "superlu_4.2" ]
    then
        echo "Replacing SuperLU_4.2/make.inc"
        cp superlu_4.2_make.inc SuperLU_4.2/make.inc
	cd SuperLU_4.2
        echo "Build blaslib"
        make blaslib
        echo "Build SuperLU"
        make
        echo "Copy result to Software directory '$installdir'"
        create_directory "$installdir"
        create_directory "$installdir/SRC"
        cp SRC/*.h $installdir/SRC/
        create_directory "$installdir/lib"
        cp lib/*.a $installdir/lib/
        cd ..
        echo "# Add this line to the CONFIGURE_FLAGS of your DUNE .opts file" >> addToMyOptsFile.txt
        echo "--with-blas=$installdir/lib --with-superlu=$installdir --with-superlu-lib=libsuperlu_4.2.a --with-superlu-blaslib=libblas.a" >> addToMyOptsFile.txt
        echo "# " >> addToMyOptsFile.txt
        echo "Clean up"
	rm -rvf SuperLU_4.2/
        exit
    fi
    #
    if [ "$packagename" = "superlu_4.3" ]
    then
        echo "Replacing SuperLU_4.3/make.inc"
        cp superlu_4.3_make.inc SuperLU_4.3/make.inc
	cd SuperLU_4.3
        echo "Build blaslib"
        make blaslib
        echo "Build SuperLU"
        make
        echo "Copy result to Software directory '$installdir'"
        create_directory "$installdir"
        create_directory "$installdir/SRC"
        cp SRC/*.h $installdir/SRC/
        create_directory "$installdir/lib"
        cp lib/*.a $installdir/lib/
        cd ..
        echo "# Add this line to the CONFIGURE_FLAGS of your DUNE .opts file" >> addToMyOptsFile.txt
        echo "--with-blas=$installdir/lib --with-superlu=$installdir --with-superlu-lib=libsuperlu_4.3.a --with-superlu-blaslib=libblas.a" >> addToMyOptsFile.txt
        echo "# " >> addToMyOptsFile.txt
        echo "Clean up"
	rm -rf SuperLU_4.3/
        exit
    fi
    #
    if [ "$packagename" = "metis-4.0.3" ]
    then
	make
        create_directory "$installdir"
        create_directory "$installdir/Lib"
        cp libmetis.a $installdir
        cp Lib/*.h $installdir/Lib
        exit
    fi
    #
    if [ "$packagename" = "ParMetis-3.1.1" ]
    then
	make
        exit
    fi
    #
    echo "Default change to package directory..."
    cd $packagename
    #
    if [[ "$packagename" == *"fftw-3.3"* ]]; # means: $packagename contains fftw-3.3
    then
        echo "Special configure of parallel FFTW."
	./configure --prefix=$installdir MPICC=mpicc --enable-threads --enable-mpi
        echo "Special configure of parallel FFTW done."
    elif [[ "$packagename" == *"hdf5-1.8"* ]]; # means: $packagename contains hdf5-1.8
    then
        echo "Special configure of parallel HDF5."
	./configure --prefix=$installdir CC=mpicc --enable-parallel --enable-shared
        echo "Special configure of parallel HDF5 done."
    elif [ "$packagename" = "ALUGrid-1.52" ]
    then
        echo "Special configure of ALUGrid 1.52 with metis-4.0.3 for DUNE."
	./configure --prefix=$installdir CXX=mpic++ CXXFLAGS="-DNDEBUG" --with-metis="$metis_installdir"
    elif [[ "$packagename" = "ug-3.11.0" ]];
    then
        echo "Special configure of UG for DUNE."
	./configure --prefix=$installdir --enable-dune --enable-parallel CC=g++ MPICC=mpicc CXXFLAGS="-O3 -DNDEBUG" F77=gfortran
        echo "Special configure of UG for DUNE done."
    else
        echo "Default configure."
	./configure --prefix=$installdir
        echo "Default configure done."
    fi
    make
    create_directory "$installdir"
    make install
    cd ..
    #
    echo "clean up"
    rm -rf $packagename
}




main(){
    create_directory "$softwaredir"
    tar_gz_file=$1
    tar_file=${tar_gz_file%.*}
    file=${tar_file%.*}
    echo $file
    make_install $file
    exit
}

main $1

