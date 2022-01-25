#!/bin/bash
# 190905
#---------------#
#-- run kaiju --#
#---------------#
# GoKaiju_vir_filter_V1_WS.sh -t=Virome_single -f=merged_uniq -a=pattern -c=4 -p=1
# 200510
# Work station으로 DB loacation 이 설정.
# looping 기존에는 migered_uniq 는 kraken2 output을 이용하기 위해 _classified.fastq.gz  로 설정 되어 있음.
# #echo $(echo $i | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/');
# 200528
# filter prophage reads from viral read



title=""
dir=""
pattern=""
ViralNT=""
evalue=0.00001


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
        -c=*|--cpu=*)
        cpu="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
        -p=*|--paired=*)
        paired="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi


# dir path
out=$title/kaiju_phage_out
phage=$title/1_phage/1_phage_reads
phage_out=$title/1_phage/2_phage_out
phage_mpa=$title/1_phage/3_phage_mpa
virus=$title/2_virus/1_virus_reads
virus_filtered=$title/2_virus/2_virus_filtered
bacteriatem=$title/2_virus/3_bacteria_temp
virus_out=$title/2_virus/4_virus_out
virus_mpa=$title/2_virus/5_virus_mpa

#--- path set ---#
rm -r $title
mkdir -p $title $out $mpa $phage $phage_out $phage_mpa $virus $virus_filtered $bacteriatem $virus_out $virus_mpa
# set input dir for fastq or fasta


#dir=$(basename $1)
#threads=$(basename $2)

filenames=$(ls $dir/*_R1_*.gz | xargs -n1 -I{} basename "{}" )
echo $filenames

### location of DB
# mac
#DB=(~/DB/kaijuDB/viruses/kaiju_db_viruses.fmi)
#NODE=(~/DB/kaijuDB/nodes.dmp)
#NAME=(~/DB/kaijuDB/names.dmp)
#minikraken1=(~/DB/kraken_db/minikraken1)

# Work station
DB=(/media/uhlemannlab/6bcc25c9-efd7-4cc5-ad7e-98ca99596cf41/heekuk/DB/kaijuDB/viruses/kaiju_db_viruses.fmi)
BacteriaDB=(/media/uhlemannlab/6bcc25c9-efd7-4cc5-ad7e-98ca99596cf41/heekuk/DB/kaiju_custom/bacteria/bacteria.fmi)
PhageDB=(/media/uhlemannlab/6bcc25c9-efd7-4cc5-ad7e-98ca99596cf41/heekuk/DB/kaiju_custom/prophage/prophage.fmi)
NODE=(/media/uhlemannlab/6bcc25c9-efd7-4cc5-ad7e-98ca99596cf41/heekuk/DB/kaijuDB/nodes.dmp)
NAME=(/media/uhlemannlab/6bcc25c9-efd7-4cc5-ad7e-98ca99596cf41/heekuk/DB/kaijuDB/names.dmp)
minikraken1=(/media/uhlemannlab/6bcc25c9-efd7-4cc5-ad7e-98ca99596cf41/heekuk/DB/minikraken1)


if [ "$paired" -eq 1 ]; then
    echo ""
    echo "#===   Run for single end   ===#"
    echo ""
    for filename in $filenames
    do

        i=$(echo $filename | sed 's/\(.*\)_R[1]\'"$pattern"'/\1/')
    
    echo ${i}
        echo ""
        echo "### $i Gathering total phage sequences ###"
        echo ""
        echo "kaiju -t node -f DB -a greedy -e 3 -E $evalue -i ${i}_R1$pattern  -o ${i}_kaiju.txt -v -z $cpu"
        kaiju -t ${NODE} -f ${PhageDB} -a greedy -e 3 -E $evalue -i $dir/${i}_R1$pattern   -o ${out}/${i}_kaiju.txt -z $cpu
        echo ""
        kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i $out/${i}_kaiju.txt -o  $out/${i}_kaiju_phage_classified.txt -u
#keep classified list # -wE를 하면 정확히 C 만 골라낼수 있다.
        grep -v "^@" $out/${i}_kaiju_phage_classified.txt | cut -f2 > $out/${i}_kaiju_phage_classified_list

    # remove by list (unalign= eukaryotic virus)
        echo ""
        echo "### $i Sorting by phage sequences list ###"
        echo ""
        seqkit grep -f $out/${i}_kaiju_phage_classified_list -v $dir/${i}_R1$pattern -o $virus/${i}_R1_virus.fastq.gz
        echo "### $i removing bacterial reads ###"
        kaiju -t ${NODE} -f ${BacteriaDB} -a greedy -e 3 -E $evalue -i $virus/${i}_R1_virus.fastq.gz -o ${bacteriatem}/${i}_kaiju_bacteria.txt -z $cpu
        kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i ${bacteriatem}/${i}_kaiju_bacteria.txt -o  ${bacteriatem}/${i}_kaiju_bacteria_classified.txt -u
        grep -v "^@" ${bacteriatem}/${i}_kaiju_bacteria_classified.txt | cut -f2 > ${bacteriatem}/${i}_kaiju_bacteria_classified_list
        seqkit grep -f ${bacteriatem}/${i}_kaiju_bacteria_classified_list -v $virus/${i}_R1_virus.fastq.gz -o $virus_filtered/${i}_R1_filtered_virus.fastq.gz
        
        echo "### $i get virus taxa ###"
        kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E $evalue -i ${virus_filtered}/${i}_R1_filtered_virus.fastq.gz -o ${virus_out}/${i}_kaiju_virus.txt -z $cpu
        kraken-mpa-report --db $minikraken1 ${virus_out}/${i}_kaiju.txt > ${virus_mpa}/${i}.mpa.txt
        
    # extract by list (align= phage)
        echo "### $i Sorting phage reads ###"
        seqkit grep -f $out/${i}_kaiju_phage_classified_list $dir/${i}_R1$pattern  -o ${phage}/${i}_R1_phage.fastq.gz
        kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E $evalue -i ${phage}/${i}_R1_phage.fastq.gz -o ${phage_out}/${i}_kaiju_phage.txt -z $cpu
        kraken-mpa-report --db $minikraken1 ${phage_out}/${i}_kaiju_phage.txt > ${phage_mpa}/${i}.mpa.txt
    done
    echo "# Merge mpa files"
    merge_metaphlan_tables.py ${virus_mpa}/*.mpa.txt > $title/kaiju_virus_mpa.txt
    merge_metaphlan_tables.py ${phage_mpa}/*.mpa.txt > $title/kaiju_phage_mpa.txt

#=== paired end ====
elif [ "$paired" -eq 2 ]; then
    echo ""
    echo "#===   Run for paired end   ===#"
    echo ""
    for filename in $filenames
    do
        i=$(echo $filename | sed 's/\(.*\)_R[12]\'"$pattern"'/\1/')
        echo ""
        echo ${i}
        echo ""
        echo "### $i Gathering total viral sequences ###"
        echo ""
        echo "kaiju -t NODE -f DB -a greedy -e 3 -E $evalue -i ${i}_R1$pattern  -j ${i}_R2$pattern  -o ${i}_kaiju.txt -v -z $cpu"
        echo ""
        kaiju -t ${NODE} -f ${PhageDB} -a greedy -e 3 -E $evalue -i $dir/${i}_R1$pattern  -j $dir/${i}_R2$pattern  -o ${out}/${i}_kaiju.txt -z $cpu
        kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i $out/${i}_kaiju.txt -o  $out/${i}_kaiju_phage_classified.txt -u
        #keep classified list
        grep -v "^@" $out/${i}_kaiju_phage_classified.txt | cut -f2 > $out/${i}_kaiju_phage_classified_list
        
        # remove by list (unalign= eukaryotic virus)
        echo ""
        echo "### $i Sorting by phage sequences list ###"
        echo ""
        seqkit grep -f $out/${i}_kaiju_phage_classified_list -v $dir/${i}_R1$pattern  -o ${virus}/${i}_R1_virus.fastq.gz
        seqkit grep -f $out/${i}_kaiju_phage_classified_list -v $dir/${i}_R2$pattern  -o ${virus}/${i}_R2_virus.fastq.gz
        
        echo "### $i removing bacterial reads ###"
        kaiju -t ${NODE} -f ${BacteriaDB} -a greedy -e 3 -E $evalue -i ${virus}/${i}_R1_virus.fastq.gz -j ${virus}/${i}_R2_virus.fastq.gz -o ${bacteriatem}/${i}_kaiju_bacteria.txt -z $cpu
        kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i ${bacteriatem}/${i}_kaiju_bacteria.txt -o  ${bacteriatem}/${i}_kaiju_bacteria_classified.txt -u
        grep -v "^@" ${bacteriatem}/${i}_kaiju_bacteria_classified.txt | cut -f2 > ${bacteriatem}/${i}_kaiju_bacteria_classified_list
        
        seqkit grep -f ${bacteriatem}/${i}_kaiju_bacteria_classified_list -v ${virus}/${i}_R1_virus.fastq.gz  -o ${virus_filtered}/${i}_R1_virus_filtered.fastq.gz
        seqkit grep -f ${bacteriatem}/${i}_kaiju_bacteria_classified_list -v ${virus}/${i}_R2_virus.fastq.gz  -o ${virus_filtered}/${i}_R2_virus_filtered.fastq.gz
        
        echo "### $i get virus taxa ###"
        kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E $evalue -i ${virus_filtered}/${i}_R1_virus_filtered.fastq.gz  -j ${virus_filtered}/${i}_R2_virus_filtered.fastq.gz -o ${virus_out}/${i}_kaiju_virus.txt -z $cpu
        kraken-mpa-report --db $minikraken1 ${virus_out}/${i}_kaiju_virus.txt > ${virus_mpa}/${i}.mpa.txt
        
        # extract by list (align= phage)
        seqkit grep -f $out/${i}_kaiju_phage_classified_list $dir/${i}_R1$pattern  -o ${phage}/${i}_R1_phage.fastq.gz
        seqkit grep -f $out/${i}_kaiju_phage_classified_list $dir/${i}_R2$pattern  -o ${phage}/${i}_R2_phage.fastq.gz
        
        kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E $evalue -i ${phage}/${i}_R1_phage.fastq.gz  -j ${phage}/${i}_R2_phage.fastq.gz -o ${phage_out}/${i}_kaiju_phage.txt -z $cpu
        kraken-mpa-report --db $minikraken1 ${phage_out}/${i}_kaiju_phage.txt > ${phage_mpa}/${i}.mpa.txt
        
    done
    
    echo ""
    echo "# Merge mpa files"
    echo ""
    merge_metaphlan_tables.py ${virus_mpa}/*.mpa.txt > $title/kaiju_virus_mpa.txt
    merge_metaphlan_tables.py ${phage_mpa}/*.mpa.txt > $title/kaiju_phage_mpa.txt
    
else
    echo "Add paired information 1 or 2"
fi
