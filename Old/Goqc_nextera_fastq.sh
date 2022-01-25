#!/bin/bash
# 190822
# Goqc_for_nextera_fastq.sh [input dir name]
# requrement:
# conda install -c bioconda bowtie2
# conda install -c bioconda cutadapt
# sudo conda install -c bioconda trim-galore



# set usage
if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi



dir=$(basename $1)
filenames=$(ls $dir/*.gz | xargs -n1 -I{} basename "{}" )
echo $filenames

DB=$2
touch -c "$DB"/*
ls -la "$DB"
echo $DBdir
#DB=$($DBdir/*.bt2)






#DB=(~/DB/human_index/GRCh38/GRCh38)
echo $DB
#echo $filenames


#DB=(~/DB/human_index/GRCh38/GRCh38)
#DB=(/Users/heekukpark/DB/mouse_index/mm10)

rm -r QC
mkdir -p QC/filtered_fastq

for filename in $filenames
do

i=$(echo $filename | sed 's/\(.*\)_R1_00[1234]\.fastq.gz/\1/')

#echo $filename
#echo ${i}

mkdir -p QC/qc_temp/${i}.filter_temp
trim_galore --fastqc_args "--outdir=QC/qc_temp/${i}.filter_temp" --stringency 1 --output_dir QC/qc_temp/${i}.filter_temp $dir/${i}_R1_001.fastq.gz >> log.txt


bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $DB QC/qc_temp/${i}.filter_temp/${i}_R1_001_trimmed.fq.gz -S QC/qc_temp/${i}.filter_temp/Homo_sapiens.sam --un-gz QC/filtered_fastq/${i}_R1_filtered.fastq.gz --al-gz QC/qc_temp/${i}.filter_temp/${i}_R1_host.fastq.gz >> log.txt

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



