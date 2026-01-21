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
cd .. 

type=ntau_10GeV_10

basedir=../gpfs/data/skimmed/pandora/${type}/awkd_testSamples

mkdir -p $basedir/concat

num=900

for file in `cat filelists/ntau_10GeV_10_testSamples.txt`; do
    # echo $file
    filename_d=${file%.*}
    filename=${filename_d##*/}
    # name=${filename%${last}}
    # echo $filename_d
    # echo $filename
    # a=$((i*10))

    # fileNum_5=`find ${basedir}/${filename} -name ${filename}_???.h5 | wc -l`
    # fileNum_1=`find ${basedir}/${filename} -name ${filename}_???_?.h5 | wc -l`
    fileNum_5=`ls ${basedir}/${filename}_???.h5 | wc -w`
    # fileNum_1=`ls ${basedir}/${filename}/${filename}_???_?.h5 | wc -w`
    fileNum=$((fileNum_5*5))

    if [ $fileNum -ne 100 ]; then ## event数に応じて変更が必要
        echo "did not correctly generate ${filename} event number : ${fileNum}"
        # break
        let num++
        continue
    fi

    # if [ -e ${basedir}/concat/${type}_${num}.h5 ]; then 
    #     echo "${basedir}/concat/${type}_${num}.h5 already exist "
    #     let num++
    #     continue
    # fi

    echo "processing ${filename}.h5 => ${type}_${num}.h5"
    ls ${basedir}/${filename}_???.h5 > temp/temp.list

    bsub -q s -o job/skimmed/${type}/concat/output.%J -e job/skimmed/${type}/concat/errors.%J "python concat_awkward.py temp/temp.list ${basedir}/concat/${type}_${num}.h5"
    # echo bsub -q s -o job/${type}/concat/output.%J -e job/${type}/concat/errors.%J "python concat_awkward.py temp/temp.list ${basedir}/concat/${type}_${num}.h5"
    # python concat_awkward.py temp/temp.list ${basedir}/concat/uds91_${num}.h5
    # bsub -q s "python LCIO2ak2_edit.py $file ${basedir}/${filename}_${S}.h5 10 ${a} > ../gpfs/data/uds/log/${filename}_${S}.log"
    let num++
done

