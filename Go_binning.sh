#!/bin/bash
# 201013
# Go_binning.sh -t=project -c=contigsfile -f=abundantTab -b=bamDIR -a=pattern
# requirement for samtools
# sudo apt-get install libncurses5
# 각각 다른 contig를 이용한 mapping. 이렇게 하는 방법인 줄 알고 만들었지만, co-aeemble을 해야 하고 그것을 사용하여 mapping 해야 한다.

source ~/miniconda3/etc/profile.d/conda.sh
conda activate MAGs



title=""
contigsfile=""
dir=""
pattern=""



for arg in "$@"
do
    case $arg in
        -t=*|--title=*)
        title="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -c=*|--contigsfile=*)
        contigsfile="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -f=*|--file=*)
        dir="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -b=*|--bamDIR=*)
        bamDIR="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -a=*|--pattern=*)
        pattern="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo ""
echo "Usage: `basename $0` -t=[project] -c=[contigsfile] -f=[abundantTab_location] -a=[pattern]"
echo ""
exit 0
fi


# dir path
out=${title}_binning
maxbin=${title}_binning/maxbin
metabat=${title}_binning/metabat
concoct=${title}_binning/concoct
das=${title}_binning/das


#--- path set ---#
rm -r $out

mkdir -p $out $maxbin $metabat $concoct $das
    

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
    echo "### Binning using bowtie2 for $i ###"
    echo ""
    echo "### metabat2 for $i ###"
    echo ""
    jgi_summarize_bam_contig_depths --outputDepth $metabat/metabat2_depth.txt $bamDIR/${i}.bam.sorted.bam
    mkdir -p $metabat/${i}/
    metabat2 -i $contigsfile -a $metabat/metabat2_depth.txt -o $metabat/${i}/${i}
    Fasta_to_Scaffolds2Bin.sh -i $metabat/${i} -e fa > $metabat/${i}_metabat.scaffolds2bin.tsv
    echo ""
    echo "### MaxBin for $i ###"
    echo ""
    mkdir -p $maxbin/${i}/
    run_MaxBin.pl -thread 20 -contig $contigsfile  -abund $dir/${i}${pattern} -out $maxbin/${i}/${i}
    Fasta_to_Scaffolds2Bin.sh -i $maxbin/${i} -e fasta > $maxbin/${i}_maxbin.scaffolds2bin.tsv
    echo ""
    echo "### concoct for $i ###"
    echo ""
    mkdir -p $concoct/${i}/
    concoct -t 20 --composition_file $contigsfile --coverage_file $dir/${i}${pattern} -b $concoct/${i}/${i}
    perl -pe "s/,/\tconcoct./g;" $concoct/${i}/${i}_clustering_gt1000.csv > $concoct/${i}_concoct.scaffolds2bin.tsv
    
    echo ""
    echo "### DAS for $i ###"
    echo ""
    mkdir -p $das/${i}/
    DAS_Tool  -i $metabat/${i}_metabat.scaffolds2bin.tsv,$maxbin/${i}_maxbin.scaffolds2bin.tsv,$concoct/${i}_concoct.scaffolds2bin.tsv -l metabat,maxbin,concoct -c $contigsfile -o $das/${i}/${i} -t 20 --score_threshold 0
    
    echo ""
    echo ""
    
done

conda deactivate
