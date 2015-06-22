#!/bin/bash

cells=( 1x1x1 1x1x2 1x1x3 1x1x4 2x2x1 2x2x2 2x2x3 2x3x3 3x3x3 )

#go through all the combination
mkdir -p out
for i in "${cells[@]}"
do
  #The next line can be run in parallel (with & in the end), but for precise timing it run sequentially.
  /usr/bin/time -o out/cell_$i.time -f "%e" supercell -d -i FeSbO4.cif -s $i -m > out/cell_$i.out 2>/dev/null 
done

#Calculating symmetry operations number. For case of huge number of permutations 
#(cells 2x3x3 3x3x3) the calculations can be done only with permutation group fixed. 
#Otherwise, the program exit before symmetry calculations.
for i in "${cells[@]}"
do
  #The program runs in parallel.
  supercell -d -i FeSbO4.cif -s $i -p "r((Fe1|Sb1)):fixed" -m > out/cell_${i}_sym.out 2>/dev/null &
done

wait

echo -e "cell\tN_sym\tN_tot\t\t\tN_unq\tRun time, s"

for i in "${cells[@]}"
do
  N_sym=`grep -Po '[0-9]* symmetry operation found for supercell.' out/cell_${i}_sym.out | awk '{print $1}'`
  N_tot=`grep -Po '(?<=The total number of combinations is ).*' out/cell_$i.out`
  N_unq=`grep -Po '(?<=Combinations after merge: ).*' out/cell_$i.out`
  [[ -z "$N_unq" ]] && N_unq="N/A"
  RT=`head -n 1 out/cell_$i.time | grep '^[0-9\.][0-9\.]*$'`
  [[ -z "$RT" ]] && RT="N/A"
  echo -e "$i\t${N_sym}\t${N_tot}\t\t\t${N_unq}\t${RT}"
done




