#!/bin/bash

#echo makeak.sh basename nevent_per_file
#
#eventperloop=50
#basename=$1
#neventfile=$2
#nloopmax=$(($neventfile/$eventperloop-1))
#basedir=/home/ilc/suehara/gpfs/doublePG/example
#akdir=$basedir/data/$basename/awkd
datadir=../gpfs/data/uds
#inputfile=$basename.list
#
#mkdir $akdir
#mkdir $akdir/log
#
#a=0
#for file in `cat $inputfile`; do
#    a=$(($a+1))
#    echo $file, $a
#    for i in `seq 0 $nloopmax`; do
#	bsub -q s -o $akdir/log/${basename}_${a}_${i}.bsublog "./makeak_indiv.sh $basename $file $a $i $eventperloop >& $akdir/log/${basename}_${a}_${i}_bsub.log"
#	sleep 1
#    done
#done

last="00_reco_REC"

a=0
for file in `cat filelists/uds91.list`; do
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

  mkdir ${datadir}/log/${filename}
  mkdir ${datadir}/awkd/${filename}
  mkdir job/${filename}

  for i in `seq 0 99`; do
    a=$((i*5))
    S=$(printf "%03d\n" "${a}")
    # filename=${name}${S}"_reco_REC"
    # echo $filename
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 1000 0"
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 10 ${a} > log/ntau_10GeV_10/${filename}_${a}.log"
    # bsub -q s "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file ../gpfs/data/ntau_10GeV_10/awkd/${filename}.h5 10 ${a} > ../gpfs/data/ntau_10GeV_10/log/${filename}.log"
    echo ${filename}_${S}.h5
    bsub -q s -o job/${filename}/output.%J -e job/${filename}/errors.%J "python LCIO2ak2_edit.py $file ${datadir}/awkd/${filename}/${filename}_${S}.h5 5 ${a} > ${datadir}/log/${filename}/${filename}_${S}.log"
    # echo "finish"
  done
  # sleep 100
done
