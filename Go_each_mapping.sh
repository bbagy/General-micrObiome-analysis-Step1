#!/bin/bash
# 201013
# Go_mapping.sh -t=project -f=merged_uniq -a=pattern -g=contigs -b=contigs_pattern
# requirement for samtools
# sudo apt-get install libncurses5
# 각각 다른 contig를 이용한 mapping. 이렇게 하는 방법인 줄 알고 만들었지만, co-aeemble을 해야 하고 그것을 사용하여 mapping 해야 한다.

source ~/miniconda3/etc/profile.d/conda.sh
conda activate shotgun



title=""
dir=""
pattern=""
contigsDir=""
contigsPattern=""

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
        -g=*|--contigsDir=*)
        contigsDir="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -b=*|--contigsPattern=*)
        contigsPattern="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo ""
echo "Usage: `basename $0` -t=project -f=merged_uniq -a=pattern -g=contigs -b=contigs_pattern"
echo ""
exit 0
fi


# dir path
out=$title

#--- path set ---#
rm -r $out

mkdir -p $out/sam $out/index $out/bam $out/count

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
    echo "### Creating bowtie2-index for $i ###"
    echo ""
    bowtie2-build $contigsDir/${i}$contigsPattern $out/index/${i}
    
    echo ""
    echo "### Mapping using bowtie2 for $i ###"
    echo ""
    
    bowtie2 -x $out/index/${i} -U $dir/${i}${pattern} -S $out/sam/${i}.sam
    
    echo ""
    echo "### Sam to bam for $i ###"
    echo ""
    samtools import $contigsDir/${i}$contigsPattern $out/sam/${i}.sam  $out/bam/${i}.bam
    samtools sort $out/bam/${i}.bam $out/bam/${i}.bam.sorted
    samtools index $out/bam/${i}.bam.sorted.bam
    
    echo ""
    echo "### Getting count for $i ###"
    echo ""
    
    samtools idxstats $out/bam/${i}.bam.sorted.bam > $out/bam/${i}.idxstats
    cut -f1,3 $out/bam/${i}.idxstats > $out/count/${i}.counts
    
    
done

conda deactivate
