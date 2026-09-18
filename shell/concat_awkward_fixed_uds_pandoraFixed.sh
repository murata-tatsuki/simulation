#!/bin/bash

cd ..

partType=ss     # uu dd ss
# basedir=/gpfs/group/ilc/users/murata/pfa/data/skimmed/pandora/fixed_uds_pandoraFixed/${partType}/awkd
# jobdir=job/skimmed/pandora/fixed_uds_pandoraFixed/${partType}
basedir=/gpfs/group/ilc/users/murata/pfa/data/skimmed/pandora/fixed_uds_pandoraFixed_nobrems/${partType}/awkd
jobdir=job/skimmed/pandora/fixed_uds_pandoraFixed_nobrems/${partType}

mkdir -p ${basedir}/concat
mkdir -p ${jobdir}/concat

# =Pn23n23h_

rm -rf ${jobdir}/concat/filelists
mkdir -p ${jobdir}/concat/filelists

a=0
for file in `cat filelists/pfa_${partType}.list`; do
    # echo $file
    filename_d=${file%.*}
    filename=${filename_d##*/}
    # name=${filename%${last}}
    # echo $filename_d
    # echo $filename
    # a=$((i*10))
    name=${filename#*"n001_"}
    num=${name:0:3}
    # name="${name:0:5}_${name:5:2}"
    # name=`echo ${name//./_}`

    # echo ${filename}
    ene_d=${filename##*${partType}}
    ene=${ene_d%%.*}
    # echo ${ene}

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
    # echo ${nevent}
    let nevent++
    # nevent=$((nevent*5))
    # echo ${nevent}
    
    mkdir -p ${basedir}/concat/${ene}GeV

    fileNum=$(/bin/ls -1 "${basedir}/${filename}/" | wc -l)

    if [ $fileNum -ne $nevent ]; then 
        echo "did not correctly generate ${filename}"
        continue
    fi

    if [ -e ${basedir}/concat/${ene}GeV/${partType}_${num}.h5 ]; then 
        echo "${basedir}/concat/${ene}GeV/${partType}_${num}.h5 already exist " $fileNum
        continue;
    fi

    mkdir -p ${jobdir}/concat/filelists/${ene}GeV
    listfile=${jobdir}/concat/filelists/${ene}GeV/${partType}_${num}.txt

    echo "processing ${filename}.h5 => ${partType}_${num}.h5"
    printf '%s\n' "${basedir}/${filename}/${filename}_"*.h5 > "${listfile}"

    nlist=$(wc -l < "${listfile}")
    if [ "$nlist" -ne "$nevent" ] || grep -q '^$' "${listfile}"; then
        echo "filelist incomplete for ${filename}: ${nlist} != ${nevent}, skip bsub"
        continue
    fi

    bsub -q s -o ${jobdir}/concat/output.%J -e ${jobdir}/concat/errors.%J "python concat_awkward.py ${listfile} ${basedir}/concat/${ene}GeV/${partType}_${num}.h5"
    # python concat_awkward.py ${jobdir}/concat/filelists/${ene}GeV/${partType}_${num}.txt ${basedir}/concat/uds91_${num}.h5
    # bsub -q s "python LCIO2ak2_edit.py $file ${basedir}/${filename}_${S}.h5 10 ${a} > ../gpfs/data/uds/log/${filename}_${S}.log"
done

