#!/bin/bash

cur_path=$PWD

if [[ -z "$PRG_GULP" ]]; then
  PRG_GULP="gulp"
fi

cells=( 1x1x2 1x2x2 )
#cells=( 1x1x2 )
rcuts=( 4.51  6.39 7.83 9.037)

for i in "${cells[@]}"
do
  rm -rf cell_$i
  mkdir cell_$i
  supercell -i PbSnTe2.cif -s $i -m -v 2 -o cell_$i/PbSnTe_c${i} > cell_$i/cell_${i}.out 2>&1  &
done

wait

for i in "${cells[@]}"
do
  cd ${cur_path}
  sqsf="${cur_path}/SQS-$i"
  echo -ne "cfg" > $sqsf
  for k in "${rcuts[@]}"
  do 
    echo -ne "\t$k" >> $sqsf
  done
  echo "" >> $sqsf
  
  for j in cell_${i}/*.cif
  do
    pth=`dirname $j`/`basename $j .cif`
    name=`basename $j .cif`

    echo -n $j >> $sqsf

    cd ${cur_path}
    rm -rf $pth/gulp
    mkdir -p $pth/gulp
    cd $pth/gulp

    frac=`sed -nr 's/^\s*\w+\s+(\w+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+).*$/\1 core \2 \3 \4 \\\n/p' ${cur_path}/$j | tr -d "\n" `
    cells=`sed -nr 's/Size a=([0-9\.]+), b=([0-9\.]+), c=([0-9\.]+)/\1 \2 \3/pg' ${cur_path}/cell_$i/cell_${i}.out`

    sed "s/@CELLABC@/$cells/g" ${cur_path}/PbSnTe-SQS.gin_template | sed "s/@NAME@/$name/g" | \
    sed "s/@FRAC@/$frac/g" > input_r.gulp 
    
    for k in "${rcuts[@]}"
    do 
      rmin=`echo $k | awk '{print $1-0.05}'`
      rmax=`echo $k | awk '{print $1+0.05}'`
      cat input_r.gulp | sed "s/@RMIN@/$rmin/g" |  sed "s/@RMAX@/$rmax/g" > input_r${k}.gulp
      ${PRG_GULP} < input_r${k}.gulp > g_r${k}.out
      cr=`grep "Interatomic potentials     =" g_r${k}.out | awk '{print int($4+0.1) };'`
      echo -ne "\t$cr" >> $sqsf
    done
    echo "" >> $sqsf
  done 
done