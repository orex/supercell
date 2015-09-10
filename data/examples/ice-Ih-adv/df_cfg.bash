#!/bin/bash

cur_path=$PWD

#Prepearing folders for data
rm -rf ice-Ih-cfgs-l1
mkdir ice-Ih-cfgs-l1

rm -rf ice-Ih-cfgs-l2
mkdir ice-Ih-cfgs-l2

rm -rf ice-Ih-cfgs-final
mkdir ice-Ih-cfgs-final


if [[ -z "$PRG_GULP" ]]; then
  PRG_GULP="gulp"
fi  

echo "Level 1: Permutations for H1-H4 groups" 

#Generating partially disordered configurations. Groups H5-H8 are excluded from permuation.
supercell -i ice-Ih-121-adv.cif -m -q -p "O:c=-2" -p "H*:c=+1" -p "r(H[5-8]):fixed" -o ice-Ih-cfgs-l1/ice-Ih-l1 > ice-Ih-cfgs-l1/log.out

cd ${cur_path}

# #Loop over all structures obtained before with supercell program
for i in ice-Ih-cfgs-l1/ice-Ih-l1*.cif
do
  pth=`dirname $i`/`basename $i .cif`
  name=`basename $i .cif`

  cd ${cur_path}
  rm -rf $pth/gulp
  mkdir -p $pth/gulp
  cd $pth/gulp

  #Prepearing input file for GULP. Only fully occupied positions will be included to the file.
  frac=`sed -nr 's/^\s*\w+\s+(\w+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+).*1.000$/\1 core \2 \3 \4 \\\n/p' ${cur_path}/$i | tr -d "\n" `
  sed "s/@FRAC@/$frac/g" ${cur_path}/ice-Ih.gin_template | sed "s/@NAME@/$name/g" > input.gulp
  $PRG_GULP < input.gulp > g.out
  
  #Checking H-H distance. If less than 1.0 Angstrom (see ice-Ih.gin_template file) the configuration is exculded.
  x=`cat g.out | grep -A 40 "Distance calculation :" | grep "   H     core        H     core         1"`
  if [[ -z "$x" ]]; then
    cp -a ${cur_path}/$i ${cur_path}/ice-Ih-cfgs-l2/.
  fi
done

echo "Level 2: Permutations for H5-H8 groups" 

cd ${cur_path}

#Loop over all "good" structures obtained in previous step.
for i in ice-Ih-cfgs-l2/ice-Ih*.cif
do
  pth=`dirname $i`/`basename $i .cif`
  name=`basename $i .cif`
  
  cd ${cur_path}
  rm -rf $pth
  mkdir -p $pth  
  
  #Generate new structures, based on "good".
  supercell -i $i -m -q -p "O:c=-2" -p "H*:c=+1"  -o $pth/${name}-l2 > $pth/log.out
  
  cd ${cur_path}/$pth  
  
  for j in ice*.cif
  do
    cd ${cur_path}/$pth  
    dm=`basename $j .cif`
    rm -rf $dm
    mkdir -p $dm/gulp
    cd $dm/gulp
    frac=`sed -nr 's/^\s*\w+\s+(\w+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+).*1.000$/\1 core \2 \3 \4 \\\n/p' ${cur_path}/$pth/$j | tr -d "\n" `
    sed "s/@FRAC@/$frac/g" ${cur_path}/ice-Ih.gin_template | sed "s/@NAME@/$dm/g" > input.gulp
    $PRG_GULP < input.gulp > g.out
    
    
    #Checking H-H distance. If less than 1.0 Angstrom (see ice-Ih.gin_template file) the configuration is exculded.
    x=`cat g.out | grep -A 40 "Distance calculation :" | grep "   H     core        H     core         1"`
    if [[ -z "$x" ]]; then
      cp -a ${cur_path}/$pth/$j ${cur_path}/ice-Ih-cfgs-final/.
    fi
    
  done
done




