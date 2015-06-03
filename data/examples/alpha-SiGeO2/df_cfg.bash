#!/bin/bash

for i in {0..6}
do
  mkdir -p cell112/Si$i
  supercell -i alpha-SiGeO2.cif -s 1x1x2 -p "Si1:p=$i" -p "Ge1:p=$((6-$i))" -m -o cell112/Si$i/SiGeO2_112-Si$i > cell112/Si$i/log.out &
done

wait

echo -ne "x"
for i in {0..6}; do echo -ne "\t$i/6"; done
echo ""

echo -ne "N_tot"
for i in {0..6}
do
  N_tot=`grep -Po '(?<=The total number of combinations is ).*' cell112/Si$i/log.out`
  echo -ne "\t${N_tot}" 
done
echo ""

echo -ne "N_irr"
for i in {0..6}
do
  N_irr=`grep -Po '(?<=Combinations after merge: ).*' cell112/Si$i/log.out`
  echo -ne "\t${N_irr}" 
done
echo ""

echo -ne "M"
for i in {0..6}
do
  Q=`echo cell112/Si$i/*.cif | sed 's/ *cell[^ ]*w.\([0-9]*\).cif */\1 /g'`
  echo -ne "\t${Q}" 
done
echo ""

echo -ne "N_ops"
for i in {0..6}
do
  Q=`echo cell112/Si$i/*.cif | sed 's/ *cell[^ ]*w.\([0-9]*\).cif */\1 /g'`
  unset N_ops
  for j in $Q; do
    N_ops=${N_ops}" "`echo $j | awk '{print 12/$1}'`
  done
  echo -ne "\t${N_ops:1}" 
done
echo ""








