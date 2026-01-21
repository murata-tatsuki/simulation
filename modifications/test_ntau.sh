#!/bin/bash

basename=ntau_10GeV_10
basedir=/home/ilc/suehara/gpfs/doublePG/example
akdir=data/ntau_10GeV_10/awkd
inputfile=$basename.list

#data/ntau_10GeV_10/reco
#
#a=0
#for file in `cat $inputfile`; do
#    a=$(($a+1))
#    echo $file, $a
#    for i in `seq 0 9`; do
#	bsub -q s -o /gpfs/group/ilc/users/murata/$akdir/log/${basename}_${a}_${i}.bsublog $python LCIO2ak2_edit.py ${basedir}/$file /gpfs/group/ilc/users/murata/$akdir/${basename}_${a}_${n}.h5 10 $i >& /gpfs/group/ilc/users/murata/$akdir/log/${basename}_${a}_${i}.log
#    done
#done


#python /home/ilc/suehara/gpfs/doublePG/example/LCIO2ak2.py /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/reco/ntau_10GeV_10_75000_reco_REC.slcio /home/ilc/suehara/gpfs/doublePG/example/data/ntau_10GeV_10/awkd/ntau_10GeV_10_${a}_${n}.h5 10 $i

#python LCIO2ak2.py ntau_10GeV_10_75000_reco_REC.slcio ntau_10GeV.h5 1 5

#python LCIO2ak2_edit.py ntau_10GeV_10_75000_reco_REC.slcio ntau_10GeV.h5 1 0
#python LCIO2ak2_edit.py ntau_10GeV_10_75000_reco_REC.slcio ntau_10GeV.h5 100 0
python LCIO2ak2_edit.py ntau_10GeV_10_75000_reco_REC.slcio test.h5 1 0

