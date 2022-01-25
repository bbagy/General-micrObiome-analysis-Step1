#!/bin/bash
# 201023
# no bowtie filter
# Goqc_illumina.sh [db path] [paired 1/2]
# requrement:
# conda install -c bioconda bowtie2
# conda install -c bioconda cutadapt
# sudo conda install -c bioconda trim-galore
# update 200331
source ~/miniconda3/etc/profile.d/conda.sh
conda activate shotgun


# set usage
if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi


#input
dir=$(basename $1)
filenames=$(ls $dir/*.gz | xargs -n1 -I{} basename "{}" )
samples=$(for filename in $filenames; do i=$(echo $filename | sed 's/\(.*\)_R[12]_001\.fastq.gz/\1/');  echo ${i}; done|sort|uniq)
echo $samples

# determind pairing
paired=$2

# run scripts
rm -r QC
mkdir -p QC/filtered_homopolymer


# run for single end
if [ "$paired" -eq 1 ]; then
for i in ${samples[*]}
do
    mkdir -p QC/qc_temp/${i}.filter_temp
    echo "starting trim_galore"
    trim_galore --fastqc_args "--outdir=QC/qc_temp/${i}.filter_temp" --stringency 1 --output_dir QC/qc_temp/${i}.filter_temp $dir/${i}_R1_001.fastq.gz >> log.txt
    
    # homopolymer
    bbduk.sh in=QC/merged_uniq/${i}_merged_unique.fastq.gz  out=QC/filtered_homopolymer/${i}_homo_filtered_merged_unique.fastq literal=AAAAAA,CCCCCC,GGGGGG,TTTTTT k=6 mm=f

    gzip QC/filtered_homopolymer/${i}_homo_filtered_merged_unique.fastq

    done
   
# run for paired end
elif [ "$paired" -eq 2 ]; then
for i in ${samples[*]}
do
    mkdir -p QC/qc_temp/${i}.filter_temp
    echo "starting trim_galore"
    trim_galore --paired --fastqc_args "--outdir=QC/qc_temp/${i}.filter_temp" --stringency 1 --output_dir QC/qc_temp/${i}.filter_temp $dir/${i}_R1_001.fastq.gz $dir/${i}_R2_001.fastq.gz >> log.txt
    
    # homopolymer

    bbduk.sh in1=QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz in2=QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz out1=QC/filtered_homopolymer/${i}_R1_homo_filtered.fastq out2=QC/filtered_homopolymer/${i}_R2_homo_filtered.fastq literal=AAAAAA,CCCCCC,GGGGGG,TTTTTT k=6 mm=f

    gzip QC/filtered_homopolymer/${i}_R1_homo_filtered.fastq
    gzip QC/filtered_homopolymer/${i}_R2_homo_filtered.fastq
    
    done

else
    echo "Add paired information 1 or 2"
fi
conda deactivate
