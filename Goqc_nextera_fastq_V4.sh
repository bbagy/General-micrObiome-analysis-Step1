#!/bin/bash
# 190822
# GoGrowth_rate_V3 [input dir name] [db path] [paired 1/2]
# requrement:
# conda install -c bioconda bowtie2
# conda install -c bioconda cutadapt
# sudo conda install -c bioconda trim-galore
# update 200331


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
mkdir -p QC/filtered_fastq

# run for single end
if [ "$paired" -eq 1 ]; then
for i in ${samples[*]}
do
    mkdir -p QC/qc_temp/${i}.filter_temp
    echo "starting trim_galore"
    trim_galore --fastqc_args "--outdir=QC/qc_temp/${i}.filter_temp" --stringency 1 --output_dir QC/qc_temp/${i}.filter_temp $dir/${i}_R1_001.fastq.gz >> log.txt

    echo "starting bowtie2"
    bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $DB QC/qc_temp/${i}.filter_temp/${i}_R1_001_trimmed.fq.gz -S QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam --un-gz QC/filtered_fastq/${i}_R1_filtered.fastq.gz --al-gz QC/qc_temp/${i}.filter_temp/${i}_R1_host.fastq.gz >> log.txt
    
    rm QC/qc_temp/${i}.filter_temp/${i}_R1_001_trimmed.fq.gz

    grep -v "^@" QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam | cut -f1 > QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list
    done


    # merge file
    mkdir QC/merged_fastq
    mkdir -p QC/merged_uniq


    for i in $(for i in $(ls QC/filtered_fastq/* | xargs -n1 -I{} basename "{}")
    do
    echo $(echo $i | sed 's/\(.*\)_L00[1234]_R1\_filtered.fastq.gz/\1/');
    done | sort | uniq)
    do
    cat QC/filtered_fastq/${i}_L001_R1_filtered.fastq.gz QC/filtered_fastq/${i}_L002_R1_filtered.fastq.gz QC/filtered_fastq/${i}_L003_R1_filtered.fastq.gz QC/filtered_fastq/${i}_L004_R1_filtered.fastq.gz > QC/merged_fastq/${i}_merged.fastq.gz


    clumpify.sh in=QC/merged_fastq/${i}_merged.fastq.gz out=QC/merged_uniq/${i}_merged_unique.fastq.gz dedupe

    done
   
# run for paired end
elif [ "$paired" -eq 2 ]; then
    for filename in $filenames
    do

    i=$(echo $filename | sed 's/\(.*\)_R[12]_00[1234]\.fastq.gz/\1/')

    echo $filename
    echo ${i}

    mkdir -p QC/qc_temp/${i}.filter_temp
    trim_galore --paired --fastqc_args "--outdir=QC/qc_temp/${i}.filter_temp" --stringency 1 --output_dir QC/qc_temp/${i}.filter_temp $dir/${i}_R1_001.fastq.gz $dir/${i}_R2_001.fastq.gz >> log.txt
    
    
    #bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $DB -1 QC/qc_temp/${i}.filter_temp/${i}_R1_001_trimmed.fq.gz -2 QC/qc_temp/${i}.filter_temp/${i}_R2_001_trimmed.fq.gz -S QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam
    
    bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $DB -1 QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz  -2 QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz -S QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam

    
    grep -v "^@" QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam | cut -f1 > QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list
    
    seqkit grep -f QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list -v QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz -o QC/filtered_fastq/${i}_R1_filtered.fastq.gz
    
    seqkit grep -f QC/qc_temp/${i}.filter_temp/Homo_sapiens.contam_list -v QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz -o QC/filtered_fastq/${i}_R2_filtered.fastq.gz

     rm QC/qc_temp/${i}.filter_temp/${i}_R1_001_val_1.fq.gz QC/qc_temp/${i}.filter_temp/${i}_R2_001_val_2.fq.gz
    
    done


    # merge file
    mkdir QC/merged_fastq
    mkdir -p QC/merged_uniq


    for i in $(for i in $(ls QC/filtered_fastq/* | xargs -n1 -I{} basename "{}")
    do
    echo $(echo $i | sed 's/\(.*\)_L00[1234]_R[12]\_filtered.fastq.gz/\1/');
    done | sort | uniq)
    do
    cat QC/filtered_fastq/${i}_L001_R1_filtered.fastq.gz QC/filtered_fastq/${i}_L002_R1_filtered.fastq.gz QC/filtered_fastq/${i}_L003_R1_filtered.fastq.gz QC/filtered_fastq/${i}_L004_R1_filtered.fastq.gz > QC/merged_fastq/${i}_R1_merged.fastq.gz
    
    cat QC/filtered_fastq/${i}_L001_R2_filtered.fastq.gz QC/filtered_fastq/${i}_L002_R2_filtered.fastq.gz QC/filtered_fastq/${i}_L003_R2_filtered.fastq.gz QC/filtered_fastq/${i}_L004_R2_filtered.fastq.gz > QC/merged_fastq/${i}_R2_merged.fastq.gz
    

    clumpify.sh in1=QC/merged_fastq/${i}_R1_merged.fastq.gz in2=QC/merged_fastq/${i}_R2_merged.fastq.gz out1=QC/merged_uniq/${i}_R1_merged_unique.fastq.gz out2=QC/merged_uniq/${i}_R2_merged_unique.fastq.gz dedupe
    


    done
else
    echo "Add paired information 1 or 2"
fi







