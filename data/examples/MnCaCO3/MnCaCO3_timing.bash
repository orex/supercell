#!/bin/bash

#go through all the combination

function get_table_row {
  supercell -d -i MnCaCO3.cif -s $1 -m -p "Mn:p=$2" -p "Ca:p=$3" | ts %.s > pp.tmp

  N_sym=`sed -rn 's/^[0-9\.]+ ([0-9]+) symmetry operation found for supercell.$/\1/gp' pp.tmp`
  N_tot=`sed -rn 's/^[0-9\.]+ The total number of combinations is ([0-9]+).*$/\1/gp' pp.tmp`
  N_unq=`sed -rn 's/^[0-9\.]+ Combinations after merge: (.*)$/\1/gp' pp.tmp`
  
  x=`echo $2 $3 | awk '{printf "%.4f", $2 / ($1+$2)}'`
  CO=`echo $2 $3 | awk '{print $1+$2}'`
  T1=`sed -rn 's/^([0-9\.]+) [0-9]+ symmetry operation found for supercell.$/\1/gp' pp.tmp`
  T2=`sed -rn 's/^([0-9\.]+) Combinations after merge: .*$/\1/gp' pp.tmp`
  RT=`echo $T1 $T2 | awk '{printf "%.2f", $2-$1}'`
  
  echo -e "A${2}B${3}(CO)$CO\t$x\t${N_tot}\t${N_unq}\t${RT}"
}

echo -e "Composition\tx\tW\tM\truntime"

for i in {0..12}
do
  get_table_row 2x2x1 $((24-$i)) $i
done

for i in {1..6}
do
  get_table_row 3x3x1 $((54-$i)) $i
done



