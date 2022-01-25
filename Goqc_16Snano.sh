#!/bin/bash
# 20211220
# Goqc_16Snano.sh [input dir name]

source ~/miniconda3/etc/profile.d/conda.sh
conda activate shotgun


# set usage
if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi


#input
maindir=$(basename $1)

# run scripts

rm 1_pooled_fastq 2_unique_fastq 3_1500bp_fastq
mkdir 1_pooled_fastq 2_unique_fastq 3_1500bp_fastq
dirs=$(ls $maindir | xargs -n1 -I{} basename "{}" )


for subdir in $dirs
do
    cat $maindir/$subdir/*.gz > 1_pooled_fastq/${subdir}_fastq.gz
    
    clumpify.sh in=1_pooled_fastq/${subdir}_fastq.gz out=2_unique_fastq/${subdir}_unique.fastq.gz dedupe
    
    seqkit subseq -r 21:1600 2_unique_fastq/${subdir}_unique.fastq.gz > 3_1500bp_fastq/${subdir}_unique.1500bp.fastq.gz
done

conda deactivate



