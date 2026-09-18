#!/bin/bash

cd ..

partType=nnqq_2M
basedir=../gpfs/data/skimmed/pandora/${partType}/awkd
jobdir=job/skimmed/pandora/${partType}

mkdir -p ${basedir}/concat
mkdir -p ${jobdir}/concat
mkdir -p ${jobdir}/concat/filelists

# =Pn23n23h_

num=0
a=0
for file in `cat filelists/nnqq_2M_0_dd.txt`; do
    # echo $file
    filename_d=${file%.*}
    filename=${filename_d##*/}
    # name=${filename%${last}}
    # echo $filename_d
    # echo $filename
    # a=$((i*10))
    name=${filename#*"n23n23h_"}
    ft=${name:5:1}
    name="${name:0:6}"
    # name="${name:3:3}_${name:0:2}"
    # name=`echo ${name//./_}`

    # echo ${name}
    # echo ${ft}

    fileNum=`ls ${basedir}/${filename}/ | wc -w`

    if [ $ft -eq 0 -a $fileNum -ne 20 ]; then 
        echo "did not correctly generate ${filename}"
        let num++
        continue
    elif [ $ft -ne 0 -a $fileNum -ne 200 ]; then 
        echo "did not correctly generate ${filename}"
        let num++
        continue
    fi

    if [ -e ${basedir}/concat/nnqq_${num}_${name}.h5 ]; then 
        echo "${basedir}/concat/nnqq_${num}_${name}.h5 already exist "
        let num++
        continue;
    fi

    rm ${jobdir}/concat/filelists/nnqq_${num}_${name}.txt
    touch ${jobdir}/concat/filelists/nnqq_${num}_${name}.txt

    echo "processing ${filename}.h5 => nnqq_${num}_${name}.h5"
    ls ${basedir}/${filename}/${filename}_*.h5 > ${jobdir}/concat/filelists/nnqq_${num}_${name}.txt

    bsub -q s -o ${jobdir}/concat/output.%J -e ${jobdir}/concat/errors.%J "python concat_awkward.py ${jobdir}/concat/filelists/nnqq_${num}_${name}.txt ${basedir}/concat/nnqq_${num}_${name}.h5"
    # python concat_awkward.py ${jobdir}/concat/filelists/nnqq_${num}_${name}.txt ${basedir}/concat/uds91_${num}.h5
    # bsub -q s "python LCIO2ak2_edit.py $file ${basedir}/${filename}_${S}.h5 10 ${a} > ../gpfs/data/uds/log/${filename}_${S}.log"
    let num++
done

