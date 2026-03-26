#!/bin/bash

cd ..

datadir=../gpfs/data/skimmed/pandora/neutron

last="000_reco_REC"

a=0
for file in `cat filelists/neutron.list`; do
  echo $file
  filename_d=${file%.*}
  filename=${filename_d##*/}
  # name=${filename%${last}}
  # echo $filename_d
  # echo $filename
  # ${KUMA} | cut -c 1-${LN}
  # name=${filename} | cut -c 1-${LN}

  # python LCIO2ak2_edit.py $file ../gpfs/data/uds/awkd/${filename}.h5 1 0 > ../gpfs/data/uds/log/${filename}.log
  # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py $file ../gpfs/data/uds/awkd/${filename}.h5 10 0 > ../gpfs/data/uds/log/${filename}.log"

  # a1=`echo $filename | cut -d '00' -f 1`
  # a2=`echo $filename | cut -d '00' -f 2`
  # echo $a1 $a2
  # name=${filename%%${last}}
  # echo $name


  # mkdir -p ${datadir}/log/${filename}
  # mkdir -p ${datadir}/awkd/${filename}
  # jobdir=job/skimmed/pandora/neutron/${filename}
  # mkdir -p ${jobdir}
  mkdir -p ${datadir}/log/concat
  mkdir -p ${datadir}/awkd/concat
  jobdir=job/skimmed/pandora/neutron/concat
  mkdir -p ${jobdir}

  name=${filename%${last}}
  echo ${filename}

  if [ -e ${datadir}/awkd/concat/${name}.h5 ]; then 
    echo "   ${datadir}/awkd/concat/${name}.h5 already exist "
    continue;
  fi

  bsub -q s -o ${jobdir}/output.%J -e ${jobdir}/errors.%J "python LCIO2ak2_brems.py ${file} ${datadir}/awkd/concat/${name}.h5 1000 ${a} > ${datadir}/log/concat/${name}.log"
  # bsub -q s -o ${jobdir}/output.%J -e ${jobdir}/errors.%J "python LCIO2ak2_single_particle.py ${file} ${datadir}/awkd/concat/${name}.h5 1000 ${a} > ${datadir}/log/concat/${name}.log"

<< COMMENTOUT
  for i in `seq 0 99`; do
    a=$((i*5))
    S=$(printf "%03d\n" "${a}")
    # filename=${name}${S}"_reco_REC"
    # echo $filename
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 1000 0"
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 10 ${a} > log/ntau_10GeV_10/${filename}_${a}.log"
    # bsub -q s "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file ../gpfs/data/ntau_10GeV_10/awkd/${filename}.h5 10 ${a} > ../gpfs/data/ntau_10GeV_10/log/${filename}.log"
    echo ${filename}_${S}.h5
    if [ -e ${datadir}/awkd/${filename}/${filename}_${S}.h5 ]; then 
        echo "   ${datadir}/awkd/${filename}/${filename}_${S}.h5 already exist "
        continue;
    fi
    echo bsub -q s -o ${jobdir}/output.%J -e ${jobdir}/errors.%J "python LCIO2ak2_neutron.py ${file} ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log"
    # echo "python LCIO2ak2_edit_skimmed_pandora.py $file ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log"
    # echo "finish"
  done
  # sleep 100
COMMENTOUT
done
