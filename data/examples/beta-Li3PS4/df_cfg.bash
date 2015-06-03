#!/bin/bash

cur_path=${PWD}

for Li2p in {0..4}
do
  Li3p=$((4-${Li2p}))
  pref=Li3PS4_Li2_${Li2p}
  mkdir ${pref}/
  cd ${pref}/
  supercell -v 2 -i Li3PS4_PNMA_2011_icsd180319.cif -m -q -o ${pref} -p "Li2:p=${Li2p}" -p "Li3:p=${Li3p}" -s 1x1x1 > log.out &
  cd ${cur_path}
done

wait
