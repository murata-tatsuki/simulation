#!/bin/bash

cd ..

partType=bb           # uu, dd, ss, uu_brems, dd_brems, ss_brems
datadir=../gpfs/data/skimmed/pandora/fixed_uds/${partType}

####  40,  91, 200, 350, 500 GeV
#### 500  500  200  200  100 events/file

a=0
for file in `cat filelists/pfa_${partType}.list`; do
  # echo $file
  filename_d=${file%.*}
  filename=${filename_d##*/}
  

  mkdir -p ${datadir}/log/${filename}
  mkdir -p ${datadir}/awkd/${filename}
  jobdir=job/skimmed/pandora/fixed_uds/${partType}/${filename}
  mkdir -p ${jobdir}
  
  echo ${filename}
  ene_d=${filename##*${partType}}
  ene=${ene_d%%.*}
  echo ${ene}

  nevent=0
  if [ 40 -eq $ene ]; then
    nevent=$((100-1))
  elif [ 91 -eq $ene ]; then
    nevent=$((100-1))
  elif [ 200 -eq $ene ]; then
    nevent=$((40-1))
  elif [ 350 -eq $ene ]; then
    nevent=$((40-1))
  elif [ 500 -eq $ene ]; then
    nevent=$((20-1))
  fi
  echo ${nevent}

  event_number=$(lcio_event_counter $file)
  if [ ${event_number} -ne $(((${nevent}+1)*5)) ]; then
    echo "wrong event numbers"
    continue
  fi


  for i in `seq 0 ${nevent}`; do
    a=$((i*5))
    S=$(printf "%03d\n" "${a}")
    # filename=${name}${S}"_reco_REC"
    # echo $filename
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 1000 0"
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 10 ${a} > log/ntau_10GeV_10/${filename}_${a}.log"
    # bsub -q s "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file ../gpfs/data/ntau_10GeV_10/awkd/${filename}.h5 10 ${a} > ../gpfs/data/ntau_10GeV_10/log/${filename}.log"
    # echo ${filename}_${S}.h5
    if [ -e ${datadir}/awkd/${filename}/${filename}_${S}.h5 ]; then 
        echo "   ${datadir}/awkd/${filename}/${filename}_${S}.h5 already exist "
        continue;
    fi
    
    bsub -q s -o ${jobdir}/output.%J -e ${jobdir}/errors.%J "python LCIO2ak2_brems.py ${file} ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log"
    # echo "python LCIO2ak2_edit_skimmed_pandora.py $file ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log"
    # echo "finish"
  done
  # sleep 100
done
