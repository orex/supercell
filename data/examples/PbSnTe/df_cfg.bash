#!/bin/bash

cells=( 1x1x1 1x1x2 1x1x3 1x2x2 2x2x2 )
xm1=( 16 8 4 2 )

mkdir -p out

for i in "${cells[@]}"
do
  M=`echo $i | sed 's/x/ /g' | awk '{print 4*$1*$2*$3}'`
  for j in "${xm1[@]}"
  do
    P=`echo "$M $j" | awk '{print $1/$2}'`
    if [[ -z `echo $P | grep "\." ` ]]; then
      /usr/bin/time -f "%e" -o out/cell_${i}_Pb1d${j}.time supercell -d -i PbSnTe2.cif -s $i -m -p "Pb1:p=$P" -p "Sn1:p=$(($M-$P))" -v 2  > out/cell_${i}_Pb1d${j}.out 2>&1 &
    fi
  done
done

echo "The process is long, please wait"

sleep 10

while [[ -z `grep "Combinations after merge" out/cell_2x2x2_Pb1d2.out` ]]; do
  sleep 5
  prg=`grep -a "Finished.*Stored .* configurations. Left .*" out/cell_2x2x2_Pb1d2.out | tail -n 1`
  echo -ne "$prg"
done
echo ""
echo ""

wait

echo -ne "N"
for i in "${cells[@]}"
do
  echo -ne "\t"`echo $i | sed 's/x/ /g' | awk '{print 8*$1*$2*$3}'`
done
echo ""

echo -ne "cell"
for i in "${cells[@]}"
do
  echo -ne "\t$i"
done
echo ""

echo -ne "N_smop"

for i in "${cells[@]}"
do
  N_sym_op=`grep -a -Po '[0-9]* symmetry operation found for supercell.' out/cell_${i}_Pb1d2.out | awk '{print $1}'`
  echo -ne "\t${N_sym_op}"
done
echo ""


for j in "${xm1[@]}"
do
  echo -ne "1/$j"
  for i in "${cells[@]}"
  do
    if [[ -f out/cell_${i}_Pb1d${j}.out ]]; then
      N_tot=`grep -a -Po '(?<=The total number of combinations is )[0-9]*' out/cell_${i}_Pb1d${j}.out`
      N_unq=`grep -a -Po '(?<=Combinations after merge: ).*' out/cell_${i}_Pb1d${j}.out`
      echo -ne "\t${N_tot}(${N_unq})"
    else
      echo -ne "\tN/A"
    fi  
  done 
  echo ""
done



#echo -e "cell\tN_sym\tN_tot\t\t\tN_unq"

#for i in "${cells[@]}"
#do
#  N_sym=`grep -Po '[0-9]* symmetry operation found for supercell.' out/cell_$i.out | awk '{print $1}'`
#  N_tot=`grep -Po '(?<=The total number of combinations is ).*' out/cell_$i.out`
#  N_unq=`grep -Po '(?<=Combinations after merge: ).*' out/cell_$i.out`
#  echo -e "$i\t${N_sym}\t${N_tot}\t\t\t${N_unq}"
#done




