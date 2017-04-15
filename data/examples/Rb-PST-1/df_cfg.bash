#!/bin/bash

rm -rf out
mkdir out

supercell -v 2 -i RB-PST-1-DEHY_1_new.cif -p "Rb2:p=0" -c yes -s 1x1x2 -p "r(Si1|Ga1):fixed" -m -o out/Rb-PST-1-cell1x1x2_stage1

cd out

for i in Rb-PST-1-cell1x1x2_stage1*.cif
do
  out_dir=`basename $i .cif`
  mkdir $out_dir
  supercell -v 2 -i $i -c try -p "Si*:c=4" -p "Ga*:c=3" -p "Rb*:c=1" -p "O*:c=-2" -o $out_dir/Rb-PST-1-cell1x1x2_stage2 -n r5000 > ${out_dir}.log &
done

echo "Please wait supercell program finish. It can take a long time depend on your system."

wait
