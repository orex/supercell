#!/bin/bash

cells=( 1x1x1 1x1x2 1x1x3 1x1x4 1x2x2 2x2x2 2x2x3 2x3x3 3x3x3 )

mkdir -p out
for i in "${cells[@]}"
do
  supercell -d -i FeSbO4.cif -s $i -m > out/cell_$i.out 2>/dev/null &
done

wait

echo -e "cell\tN_sym\tN_tot\t\t\tN_unq"

for i in "${cells[@]}"
do
  N_sym=`grep -Po '[0-9]* symmetry operation found for supercell.' out/cell_$i.out | awk '{print $1}'`
  N_tot=`grep -Po '(?<=The total number of combinations is ).*' out/cell_$i.out`
  N_unq=`grep -Po '(?<=Combinations after merge: ).*' out/cell_$i.out`
  echo -e "$i\t${N_sym}\t${N_tot}\t\t\t${N_unq}"
done




