#!/bin/bash
### Starting from joined_fasta
echo ""
echo "### Step1 collecting viral and phage seq###"
echo ""

#mkdir -p virus_sorting/sorted_viral_sequence/sorted_phageFromal
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_virFromal
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_TotalvirFromun
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_Nonvir
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_phageFromun
#mkdir -p virus_sorting/sorted_viral_sequence/sorted_virFromun
#mkdir -p virus_sorting/candidatePhage_seq
#mkdir -p virus_sorting/candidateVir_seq
#mkdir -p virus_sorting/sorted_finalvir

### viral index
#viralIndex=(~/DB/VirusSeeker/index/VirusDBNT)
### phage db
#phageDB=(~/DB/phast_db/diamond/phageDBNR_20180727)
### viral db
#viralDB=(~/DB/VirusSeeker/diamond/VirusDBNR_20160802)
### bacterial db
#bacteriaDB=(~/DB/phast_db/diamond/bacteria)



##-------- Input arg ---------##
# Default values of arguments
title=""
dir=""
ViralNT=""
ViralNR=""
phageNR=""
bacterialNR=""

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
        -v=*|--ViralNT=*)
        viralNT="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
        -r=*|--ViralNR=*)
        viralNR="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
        -p=*|--phageNR=*)
        phageNR="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
        -b=*|--bacterialNR=*)
        bacterialNR="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
    esac
done




# dir path
NTsam=(sorted_viral_sequence/sam)
alNTDIR=(sorted_viral_sequence/aligned_viral_reads_by_NT)
unNTDIR=(sorted_viral_sequence/unaligned_viral_reads_by_NT)
phageoutNRDIR=(unsorted_viral_sequence/phageNRdiamond/phageNRdiamond_out)
alphageNRDIR=(unsorted_viral_sequence/phageNRdiamond/aligned_phage_reads_by_NT)
unphageNRDIR=(unsorted_viral_sequence/phageNRdiamond/unaligned_phage_reads_by_NT)

#--- path set ---#
rm -r $title
mkidir $title

mkdir -p ${title}/${alNTDIR} $title/${unNTDIR} $title/${NTsam}

# run bowtie2
#dir=$(basename $1)
filenames=$(ls $dir/*.gz | xargs -n1 -I{} basename "{}" )
for filename in $filenames
do
i=$(echo $filename | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/')
#echo $filename
#echo ${i}
echo ""
echo "### $i Gathering total viral sequences ###"
echo ""
    bowtie2 --threads 16 --very-sensitive -k 1 --no-unal --phred33 -x $viralNT -1 $dir/${i}_R1_merged_unique.fastq.gz  -2 $dir/${i}_R2_merged_unique.fastq.gz -S $title/${NTsam}/${i}.sam
    
    grep -v "^@" $title/${NTsam}/${i}.sam | cut -f1 > $title/${NTsam}/${i}_viral_reads_by_NT_list
    
    # remove by list (unalign)
    echo ""
    echo "### $i Sorting by viral sequences list ###"
    seqkit grep -f $title/${NTsam}/${i}_viral_reads_by_NT_list -v $dir/${i}_R1_merged_unique.fastq.gz -o $title/${unNTDIR}/${i}_R1_filtered.fastq.gz
    seqkit grep -f $title/${NTsam}/${i}_viral_reads_by_NT_list -v $dir/${i}_R2_merged_unique.fastq.gz -o $title/${unNTDIR}/${i}_R2_filtered.fastq.gz
    # extract by list (align)
    seqkit grep -f $title/${NTsam}/${i}_viral_reads_by_NT_list $dir/${i}_R1_merged_unique.fastq.gz -o $title/${alNTDIR}/${i}_R1_filtered.fastq.gz
    seqkit grep -f $title/${NTsam}/${i}_viral_reads_by_NT_list $dir/${i}_R2_merged_unique.fastq.gz -o $title/${alNTDIR}/${i}_R2_filtered.fastq.gz
    
done



echo "### $i Gathering total pahge sequences from unaligned against viral NT sequences ###"
echo "### $i loading pahgeNR DB ..... ###"
mkdir -p virus_sorting/sorted_viral_sequence/phageNRdiamond/${i}.blastx_vir

diamond blastx  --evalue 1e-5 --query-cover 80 --id 50 -q virus_sorting/sorted_viral_sequence/al_viral_index/${i}_al_viral.fasta -d ${phageDB} -o virus_sorting/sorted_viral_sequence/phageNRdiamond/${i}.blastx_vir/${i}_phageDB-outfmtjuan.blastx -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

awk '!x[$1]++' virus_sorting/sorted_viral_sequence/phageNRdiamond/${i}.blastx_vir/${i}_phageDB-outfmtjuan.blastx > virus_sorting/sorted_viral_sequence/phageNRdiamond/${i}.blastx_vir/${i}_phageDB-outfmtjuan_tophits.txt

grep -v "^@" virus_sorting/sorted_viral_sequence/phageNRdiamond/${i}.blastx_vir/${i}_phageDB-outfmtjuan_tophits.txt | cut -f1 > virus_sorting/sorted_viral_sequence/phageNRdiamond/${i}.blastx_vir/${i}_phageDB-outfmtjuan_tophits_list.txt



