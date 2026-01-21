#!/bin/bash

cd ..

basedir=../gpfs/data/skimmed/pandora/neutron/awkd
jobdir=job/skimmed/pandora/neutron

mkdir -p ${basedir}/concat
mkdir -p ${jobdir}/concat

num=1
a=0
for file in `cat filelists/neutron.list`; do
    echo $file
    filename_d=${file%.*}
    filename=${filename_d##*/}
    name=${filename%${last}}
    echo $filename_d
    echo $filename
    # a=$((i*10))

    fileNum=`ls ${basedir}/${filename}/ | wc -w`

    if [ $fileNum -ne 100 ]; then 
        echo "did not correctly generate ${filename}"
        let num++
        continue
    fi

    if [ -e ${basedir}/concat/neutron_${num}.h5 ]; then 
        echo "${basedir}/concat/neutron_${num}.h5 already exist "
        let num++
        continue;
    fi

    echo "processing ${filename}.h5 => neutron_${num}.h5"
    ls ${basedir}/${filename}/${filename}_*.h5 > temp/temp.list

    bsub -q s -o ${jobdir}/concat/output.%J -e ${jobdir}/concat/errors.%J "python concat_awkward.py temp/temp.list ${basedir}/concat/neutron_${num}.h5"
    # python concat_awkward.py temp/temp.list ${basedir}/concat/uds91_${num}.h5
    # bsub -q s "python LCIO2ak2_edit.py $file ${basedir}/${filename}_${S}.h5 10 ${a} > ../gpfs/data/uds/log/${filename}_${S}.log"
    let num++
done

