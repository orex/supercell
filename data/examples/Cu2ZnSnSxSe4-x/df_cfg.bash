#!/bin/bash

for i in {0..8}
do
  mkdir -p cell111/S$i
  supercell -i stannite.cif -s 1x1x1 -p "S:p=$i" -p "Se:p=$((8-$i))" -m -o cell111/S$i/stannite-S$i > cell111/S$i/log.out &
done

wait

echo -e "x\tg_i"
for i in {0..8}
do
  x=`echo $i | awk '{print $i/2}'`
  for j in cell111/S$i/stannite-S$i*.cif
  do 
    gi=`echo $j | sed 's/ *cell[^ ]*w.\([0-9]*\).cif */\1 /g'`
    echo -e "${x}\t${gi}"
  done
done

