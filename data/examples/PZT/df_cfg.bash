#!/bin/bash

#go through all the combination
rm -rf out
mkdir -p out
#Create cell 4x2x1
supercell  -s 4x2x1 -i PZT-PbZr05Ti05O3.cif -m -o out/PZT421 > /dev/null

#Create table
for i in {0..9}
do
  echo $((i+1))
  for y in 0.25000 0.75000
  do
    for x in 0.12500 0.37500 0.62500 0.87500
    do
      vl=`grep -E "(Zr1|Ti1).+$x\s+$y\s+0.56490" out/PZT421_i0${i}_w*.cif | sed -r 's/\s+([ZT])[ri]1.*/\1/g'`
      echo -n $vl
    done
    echo ""
  done
  echo ""
done
