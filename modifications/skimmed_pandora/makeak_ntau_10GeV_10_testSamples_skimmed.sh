#!/bin/bash

cd ..

datadir=../gpfs/data/skimmed/pandora/ntau_10GeV_10

mkdir -p ${datadir}/log_testSamples/
mkdir -p ${datadir}/awkd_testSamples/
mkdir -p job/skimmed/pandora/ntau_10GeV_10

lcio_path=/home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/


for file in `cat filelists/ntau_10GeV_10_testSamples.txt`; do
  # echo $file
  filename_d=${file%.*}
  filename=${filename_d##*/}
  echo ${filename}

  for i in `seq 0 19`; do
    a=$((i*5))
    S=$(printf "%03d\n" "${a}")
    # filename=${name}${S}"_reco_REC"
    # echo $filename
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 1000 0"
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 10 ${a} > log/ntau_10GeV_10/${filename}_${a}.log"
    # bsub -q s "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file ../gpfs/data/ntau_10GeV_10/awkd/${filename}.h5 10 ${a} > ../gpfs/data/ntau_10GeV_10/log/${filename}.log"
    
    # if [ -e ${datadir}/awkd/${filename}/${filename}_${S}.h5 ]; then 
    #   echo "${datadir}/awkd/${filename}/${filename}_${S}.h5 already exist "
    #   continue;
    # fi

    echo ${filename}_${S}.h5
    # echo python LCIO2ak2_edit_skimmed_pandora.py ${lcio_path}${file} ${datadir}/awkd_testSamples/${filename}_${S}.h5 5 ${a}
    echo bsub -q s -o job/skimmed/pandora/ntau_10GeV_10/output.%J -e job/skimmed/pandora/ntau_10GeV_10/errors.%J "python LCIO2ak2_edit_skimmed_pandora.py ${lcio_path}${file} ${datadir}/awkd_testSamples/${filename}_${S}.h5 5 ${a} > ${datadir}/log_testSamples/${filename}_${S}.log"

    # bsub -q s -o job/${type}/${filename}/output.%J -e job/${type}/${filename}/errors.%J "python LCIO2ak2_edit.py $file ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log"
    # echo python LCIO2ak2_edit.py $file ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log
    # python LCIO2ak2_edit.py $file ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log
    # echo "finish"
  done

  # echo bsub -q s -o job/skimmed/ntau_10GeV_10/output.%J -e job/skimmed/ntau_10GeV_10/errors.%J "python LCIO2ak2_edit_skimmed.py ${lcio_path}${file} ${datadir}/awkd_testSamples/${filename}.h5 1000 0 > ${datadir}/log_testSamples/${filename}.log"
  # python LCIO2ak2_edit_skimmed.py ${lcio_path}${file} ${datadir}/awkd/${filename}.h5 1000 0 > ${datadir}/log/${filename}.log
done

