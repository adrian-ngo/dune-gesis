#!/bin/bash
if [ "$2" == "dbg" ];
then
    echo "build $1 debug version"
    build_flags="-O0 -g -Wall"
    build_dir="../dbg"
else
    echo "build $1 release version"
    build_flags="-O3 -funroll-loops -DNDEBUG"
    build_dir="../bin"
fi
#
make CXXFLAGS="$build_flags" $1
returnValue="$?"
if [ -f "$1" ];
then
    mv -v $1 $build_dir
    returnValue="$?"
else
    echo "File not found: $1"
fi
echo "Exit with " $returnValue
exit $returnValue
