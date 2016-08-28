#!/bin/bash

cells=( 1x1x2 1x2x2 2x2x2 )
xm1=( 0.0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1.0 )

mkdir -p out

echo "The process is long, please wait"

for i in "${cells[@]}"
do
  echo "Cell $i"
  M=`echo $i | sed 's/x/ /g' | awk '{print 4*$1*$2*$3}'`
  for j in "${xm1[@]}"
  do
    P=`echo "$M $j" | awk '{print $1*$2}'`
    supercell -i Ti1-xMgxN.cif -s $i -m -p "Mg:p=$P" -p "Ti:p=$(($M-$P))" -o out/TiMgN-$i-x$j -n r1  > /dev/null  2>&1 &
  done
  wait
done

