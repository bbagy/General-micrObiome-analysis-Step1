#!/bin/bash
### Starting from joined_fasta
echo ""
echo "### Step1 collecting viral and phage seq###"
echo ""
rm -r virus_sorting/

# dir path
NTsam=(virus_sorting/sorted_viral_sequence/sam)
alNTDIR=(virus_sorting/sorted_viral_sequence/aligned_viral_reads_by_NT)
unNTDIR=(virus_sorting/sorted_viral_sequence/unaligned_viral_reads_by_NT)

#mkdir -p virus_sorting/sorted_viral_sequence/sorted_phageFromal
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_virFromal
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_TotalvirFromun
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_Nonvir
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_phageFromun
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_virFromun
#mkdir -p virus_sorting/candidatePhage_seq
#mkdir -p virus_sorting/candidateVir_seq
#mkdir -p virus_sorting/sorted_finalvir



test input
while getopts ":f:v:" opt; do
  case $opt in
    f) dir="$OPTARG"
    ;;
    v) viralNT="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


### viral index
#viralIndex=(~/DB/VirusSeeker/index/VirusDBNT)
### phage db
#phageDB=(~/DB/phast_db/diamond/phageDBNR_20180727)
### viral db
#viralDB=(~/DB/VirusSeeker/diamond/VirusDBNR_20160802)
### bacterial db
#bacteriaDB=(~/DB/phast_db/diamond/bacteria)



# input
#dir=$(basename $1)
filenames=$(ls $dir/*.gz | xargs -n1 -I{} basename "{}" )
mkdir -p ${alNTDIR} ${unNTDIR} ${NTsam}

# run bowtie2
for filename in $filenames
do
i=$(echo $filename | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/')
#echo $filename
#echo ${i}
echo ""
echo "### $i Gathering total viral sequences ###"

    bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $viralNT -1 $dir/${i}_R1_merged_unique.fastq.gz  -2 $dir/${i}_R2_merged_unique.fastq.gz -S ${NTsam}/${i}.sam
    
    grep -v "^@" ${NTsam}/${i}.sam | cut -f1 > ${NTsam}/${i}_viral_reads_by_NT_list
    
    # remove by list (unalign)
    echo ""
    echo "### $i Sorting by viral sequences list ###"
    seqkit grep -f ${NTsam}/${i}_viral_reads_by_NT_list -v $dir/${i}_R1_merged_unique.fastq.gz -o ${unNTDIR}/${i}_R1_filtered.fastq.gz
    seqkit grep -f ${NTsam}/${i}_viral_reads_by_NT_list -v $dir/${i}_R2_merged_unique.fastq.gz -o ${unNTDIR}/${i}_R2_filtered.fastq.gz
    # extract by list (align)
    seqkit grep -f ${NTsam}/${i}_viral_reads_by_NT_list $dir/${i}_R1_merged_unique.fastq.gz -o ${alNTDIR}/${i}_R1_filtered.fastq.gz
    seqkit grep -f ${NTsam}/${i}_viral_reads_by_NT_list $dir/${i}_R2_merged_unique.fastq.gz -o ${alNTDIR}/${i}_R2_filtered.fastq.gz
    

done


