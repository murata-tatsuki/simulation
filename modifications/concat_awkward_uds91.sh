# basedir=$1
# basename=$2
# nstart=$3
# nend=$4
# mkdir $basedir/concat

# for i in `seq $nstart $nend`; do
#     echo processing ${basename}_${i}.h5 ...
#     ls $basedir/${basename}_${i}_*.h5 > temp/${i}.list
#     python concat_awkward.py temp/${i}.list $basedir/concat/${basename}_${i}.h5
# done

basedir=../gpfs/data/skimmed/uds/awkd

mkdir -p ${basedir}/concat
mkdir -p job/skimmed/uds91/concat

num=1
a=0
for file in `cat filelists/uds91.list`; do
    # echo $file
    filename_d=${file%.*}
    filename=${filename_d##*/}
    # name=${filename%${last}}
    # echo $filename_d
    # echo $filename
    # a=$((i*10))

    fileNum=`ls ${basedir}/${filename}/ | wc -w`

    if [ $fileNum -ne 100 ]; then 
        echo "did not correctly generate ${filename}"
        break
    fi

    if [ -e ${basedir}/concat/uds91_${num}.h5 ]; then 
        echo "${basedir}/concat/uds91_${num}.h5 already exist "
        let num++
        continue;
    fi

    echo "processing ${filename}.h5 => uds91_${num}.h5"
    ls ${basedir}/${filename}/${filename}_*.h5 > temp/temp.list

    # bsub -q s -o job/skimmed/uds91/concat/output.%J -e job/skimmed/uds91/concat/errors.%J "python concat_awkward.py temp/temp.list ${basedir}/concat/uds91_${num}.h5"
    python concat_awkward.py temp/temp.list ${basedir}/concat/uds91_${num}.h5
    # python concat_awkward.py temp/temp.list ${basedir}/concat/uds91_${num}.h5
    # bsub -q s "python LCIO2ak2_edit.py $file ${basedir}/${filename}_${S}.h5 10 ${a} > ../gpfs/data/uds/log/${filename}_${S}.log"
    let num++
done

