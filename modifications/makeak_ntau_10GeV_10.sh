#!/bin/bash

datadir=../gpfs/data/skimmed/ntau_10GeV_10

mkdir -p ${datadir}/log/
mkdir -p ${datadir}/awkd/
mkdir -p job/skimmed/ntau_10GeV_10


for file in `cat filelists/ntau_10GeV_10.txt`; do
  # echo $file
  filename_d=${file%.*}
  filename=${filename_d##*/}
  echo ${filename}

  # bsub -q s -o job/skimmed/ntau_10GeV_10/output.%J -e job/skimmed/ntau_10GeV_10/errors.%J "python LCIO2ak2_edit_skimmed.py ${file} ${datadir}/awkd/${filename}.h5 1000 0 > ${datadir}/log/${filename}.log"
   echo "python LCIO2ak2_edit_skimmed.py ${file} ${datadir}/awkd/${filename}.h5 10 5 > ${datadir}/log/${filename}.log"
done

