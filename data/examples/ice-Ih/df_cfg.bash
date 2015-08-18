#!/bin/bash

cur_path=$PWD

rm -rf ice-Ih-cfgs
mkdir ice-Ih-cfgs

if [[ -z "$PRG_GULP" ]]; then
  PRG_GULP="gulp"
fi  

supercell -i ice-Ih.cif -m -q -p "O:c=-2" -p "H*:c=+1" -o ice-Ih-cfgs/ice-Ih

echo -e "cfg\tCNO\tH-H" > ice-Ih-cfgs.anlz

for i in ice-Ih-cfgs/ice-Ih*.cif
do
  echo $i
  pth=`dirname $i`/`basename $i .cif`
  name=`basename $i .cif`

  cd ${cur_path}
  rm -rf $pth/gulp
  mkdir -p $pth/gulp
  cd $pth/gulp

  frac=`sed -nr 's/^\s*\w+\s+(\w+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+).*$/\1 core \2 \3 \4 \\\n/p' ${cur_path}/$i | tr -d "\n" `
  #echo $frac
  sed "s/@FRAC@/$frac/g" ${cur_path}/ice-Ih.gin_template | sed "s/@NAME@/$name/g" > input.gulp
  $PRG_GULP < input.gulp > g.out
  
  x=`cat g.out | grep -A 40 "Distance calculation :" | grep "   H     core        H     core         1"`
  if [[ -z "$x" ]]; then
    HH="OK"
  else
    HH="Fail"
  fi
  x=`cat g.out | grep -A 40 "Distance calculation :" | grep "   O     core        H     core" | awk '{print $6+$8 == 2}' | grep 0`
  if [[ -z "$x" ]]; then
    CNO="OK"
  else
    CNO="Fail"
  fi

  echo -e "$i\t$CNO\t$HH" >> ${cur_path}/ice-Ih-cfgs.anlz

done
