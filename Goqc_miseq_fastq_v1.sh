#!/bin/bash
# 190822
# Goqc_nextra_fastq [input dir name] [db path] [paired 1/2]
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

# filter index
DB=$2
echo $DB
#DB=$($DBdir/*.bt2)

# determind pairing
paired=$3

# run scripts
rm -r QC
mkdir -p QC/filtered_fastq QC/filtered_homopolymer

# run for single end
if [ "$paired" -eq 1 ]; then
for i in ${samples[*]}
do
    echo "starting bowtie2"
    mkdir -p QC/qc_temp/${i}.filter_temp
    bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $DB QC/qc_temp/${i}.filter_temp/${i}_R1_001_trimmed.fq.gz -S QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam --un-gz QC/filtered_fastq/${i}_R1_filtered.fastq.gz --al-gz QC/qc_temp/${i}.filter_temp/${i}_R1_host.fastq.gz >> log.txt
    
    rm QC/qc_temp/${i}.filter_temp/${i}_R1_001_trimmed.fq.gz

    grep -v "^@" QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam | cut -f1 > QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list
    
    seqkit grep -f QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list -v QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz -o QC/filtered_fastq/${i}_R1_filtered.fastq.gz
    
    rm QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz
    
    
    done
   
# run for paired end
elif [ "$paired" -eq 2 ]; then
for i in ${samples[*]}
do
    mkdir -p QC/qc_temp/${i}.filter_temp
    bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $DB -1 QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz  -2 QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz -S QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam

    
    grep -v "^@" QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam | cut -f1 > QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list
    
    seqkit grep -f QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list -v QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz -o QC/filtered_fastq/${i}_R1_filtered.fastq.gz
    
    seqkit grep -f QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list -v QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz -o QC/filtered_fastq/${i}_R2_filtered.fastq.gz

     rm QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz
    
    done
else
    echo "Add paired information 1 or 2"
fi
conda deactivate






