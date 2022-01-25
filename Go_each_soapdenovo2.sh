#!/bin/bash
# 201013
# 각각 시료당 contig. 이렇게 하는 방법인 줄 알고 만들었지만, co-aeemble을 해야 하고 그것을 사용하여 mapping 해야 한다.

# Go_soapdenovo2.sh -t=project -f=merged_uniq -a=pattern -c=configFile
# requrement:configFile
# vi configFile

# i
# #maximal read length
# #max_rd_len=150
# [LIB]
# #average insert size of the library
# avg_ins=300
# #if sequences are forward-reverse of reverse-forward
# reverse_seq=0
# #in which part(s) the reads are used (only contigs, only scaffolds, both contigs and scaffolds, only gap closure)
# asm_flags=3
# #cut the reads to the given length
# rd_len_cutoff=100
# #in which order the reads are used while scaffolding
# rank=1
# # cutoff of pair number for a reliable connection (at least 3 for short insert size)
# pair_num_cutoff=3
# #minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
# map_len=32
# #fastq file for single reads
# q=single_read

source ~/miniconda3/etc/profile.d/conda.sh

#    conda activate soapdenovo2
#    conda activate shotgun
title=""
dir=""
pattern=""
configFile=""

for arg in "$@"
do
    case $arg in
        -t=*|--title=*)
        title="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -f=*|--file=*)
        dir="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -a=*|--pattern=*)
        pattern="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -c=*|--configFile=*)
        configFile="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo ""
echo "Usage: `basename $0` -t=project -f=merged_uniq -a=pattern -c=configFile"
echo ""
exit 0
fi


# dir path
out=$title



#--- path set ---#
rm -r $out


filenames=$(ls $dir/*$pattern | xargs -n1 -I{} basename "{}" )
samples=$(for filename in $filenames; do i=$(echo $filename | sed 's/\(.*\)\'"$pattern"'/\1/');  echo ${i}; done|sort|uniq)
echo $samples

echo ""
echo "#===   Assemble for single end   ===#"
echo ""
for i in ${samples[*]}
do
echo ${i}
    echo ""
    echo "### creating configFile for $i ###"
    echo ""
    mkdir -p $out/${i}.assemble $out/contigs

    # option1 잘됨
    # sed 's/single_read/'"${i}${pattern}"'/' $configFile > $PWD/$out/${i}.assemble/configFile_${i}
    
    # option2 잘됨

    file=$PWD/$dir/${i}${pattern}
    echo $file
    sed "s|single_read|$file|g" $configFile > $PWD/$out/${i}.assemble/configFile_${i}
    
    conda activate soapdenovo2
    SOAPdenovo-63mer all -s  $PWD/$out/${i}.assemble/configFile_${i} -K 31 -o $out/${i}.assemble/${i}
    
    conda activate shotgun
    seqkit seq -m 200 $out/${i}.assemble/${i}.contig > $out/contigs/${i}.200.contig

done

conda deactivate
