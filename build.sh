#!/bin/bash
if [ $# == 0 ]
	then
    	echo "Providing arguments"
	exit
fi

m=$2
r=$4
mode=$6

#
echo "Recompile for m = "$m" and r = "$r
dim=$((m-r))
echo "Dimension = "$dim

sed -i "s/DIM[ ][0-9]*/DIM $dim/" include/Utils.h
DEFS="-DMODE_$mode$"

# recompile
LD_PRELOAD=libtbbmalloc_proxy.so.2
export LD_PRELOAD

#remove compiled executables 
rm GenerateInputHybridAttackBinaryError
rm GenerateInputHybridAttackRandomError

cd GenerateBinaryError
make clean
make
cp GenerateInputHybridAttackBinaryError ../

cd ../GenerateRandomError
make clean
make
cp GenerateInputHybridAttackRandomError ../

cd ..
make clean
make -j8 DEFS=$DEFS



