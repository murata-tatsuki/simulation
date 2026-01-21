basename=$1
file=$2
a=$3
i=$4

basedir=/home/ilc/suehara/gpfs/doublePG/example
akdir=data/$basename/awkd

for j in `seq 0 9`; do
    n=$(($i*10+$j))
    python $basedir/LCIO2ak2.py $file $basedir/$akdir/${basename}_${a}_${n}.h5 1 $n >& $basedir/$akdir/log/${basename}_${a}_${n}.log
done
