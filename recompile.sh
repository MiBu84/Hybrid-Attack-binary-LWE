r=20
beta=23
m=160
c=5


#
dim=$((m-r))
echo $dim

sed -i "s/DIM[ ][0-9]*/DIM $dim/" include/Utils.h

# recompile
make clean
make -j8

