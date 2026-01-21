#!/bin/bash

cd ..

partType=nnqq_brems
basedir=../gpfs/data/skimmed/pandora/${partType}/awkd
jobdir=job/skimmed/pandora/${partType}

mkdir -p ${basedir}/concat
mkdir -p ${jobdir}/concat

# =Pn23n23h_

num=101
a=0
for file in `cat filelists/nnqq_moreStats.txt`; do
    # echo $file
    filename_d=${file%.*}
    filename=${filename_d##*/}
    # name=${filename%${last}}
    # echo $filename_d
    # echo $filename
    # a=$((i*10))
    name=${filename#*"n23n23h_"}
    name="${name:0:5}_${name:5:2}"
    # name=`echo ${name//./_}`


    fileNum=`ls ${basedir}/${filename}/ | wc -w`

    if [ $fileNum -ne 20 ]; then 
        echo "did not correctly generate ${filename}"
        let num++
        continue
    fi

    if [ -e ${basedir}/concat/nnqq_${num}_${name}.h5 ]; then 
        echo "${basedir}/concat/nnqq_${num}_${name}.h5 already exist "
        let num++
        continue;
    fi

    echo "processing ${filename}.h5 => nnqq_${num}_${name}.h5"
    ls ${basedir}/${filename}/${filename}_*.h5 > temp/temp.list

    bsub -q s -o ${jobdir}/concat/output.%J -e ${jobdir}/concat/errors.%J "python concat_awkward.py temp/temp.list ${basedir}/concat/nnqq_${num}_${name}.h5"
    # python concat_awkward.py temp/temp.list ${basedir}/concat/uds91_${num}.h5
    # bsub -q s "python LCIO2ak2_edit.py $file ${basedir}/${filename}_${S}.h5 10 ${a} > ../gpfs/data/uds/log/${filename}_${S}.log"
    let num++
done

