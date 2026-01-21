#!/bin/bash

#echo makeak.sh basename nevent_per_file
#
#eventperloop=50
#basename=$1
#neventfile=$2
#nloopmax=$(($neventfile/$eventperloop-1))
#basedir=/home/ilc/suehara/gpfs/doublePG/example
#akdir=$basedir/data/$basename/awkd
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
for file in `cat filelists/ntau_10GeV_10.txt`; do
  echo $file
  filename=${file%.*}
  name=${filename%${last}}
  # ${KUMA} | cut -c 1-${LN}
  # name=${filename} | cut -c 1-${LN}

  # a1=`echo $filename | cut -d '00' -f 1`
  # a2=`echo $filename | cut -d '00' -f 2`
  # echo $a1 $a2
  # name=${filename%%${last}}
  # echo $name
  for i in `seq 0 9`; do
    a=$((i*10))
    S=$(printf "%02d\n" "${a}")
    filename=${name}${S}"_reco_REC"
    echo $filename
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 1000 0"
    # bsub -q s -o job/output.%J -e job/errors.%J "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file h5/ntau_10GeV_10/${filename}.h5 10 ${a} > log/ntau_10GeV_10/${filename}_${a}.log"
    echo bsub -q s "python LCIO2ak2_edit.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/$file ../gpfs/data/skimmed/ntau_10GeV_10/awkd/${filename}.h5 10 ${a} > ../gpfs/data/skimmed/ntau_10GeV_10/log/${filename}.log"
    # echo "finish"
  done
done
