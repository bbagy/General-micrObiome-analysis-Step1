#!/bin/bash
# 190905
#---------------#
#-- run kaiju --#
#---------------#
# kaiju 부터 kraken mpa까지
# 190905
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# GoKaiju.sh [input dir name]

title=""
dir=""
ViralNT=""
paired=""

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
out=(kaiju_out)
mpa=(mpa)
krona=(krona)
classified=(classified)
unclassified=(unclassified)



#--- path set ---#
rm -r $title

# set input dir for fastq or fasta


#dir=$(basename $1)
#threads=$(basename $2)

filenames=$(ls $dir/* | xargs -n1 -I{} basename "{}" )


# location of DB
DB=(~/DB/kaijuDB/viruses/kaiju_db_viruses.fmi)
NODE=(~/DB/kaijuDB/nodes.dmp)
NAME=(~/DB/kaijuDB/names.dmp)
minikraken1=(~/DB/kraken_db/minikraken1)



mkdir -p $title/$out $title/$mpa $title/$krona



for filename in $filenames
do
filenames=$(ls $dir/*.gz | xargs -n1 -I{} basename "{}" )
for filename in $filenames
do
i=$(echo $filename | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/')
echo $filename
echo ${i}
echo ""
echo "### $i Gathering total viral sequences ###"
echo ""
kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E 0.05 -i $dir/${i}_R1_merged_unique.fastq.gz -j $dir/${i}_R1_merged_unique.fastq.gz -o $title/$out/${i}_kaiju.txt -v -z $cpu
kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i $title/$out/${i}_kaiju.txt -o  $title/$out/${i}_kaiju_classified.txt -u

#keep classified list
grep -v "^@" $title/$out/${i}_kaiju_classified.txt | cut -f2 > $title/$out/${i}_kaiju_classified_list

# remove by list (unalign)
echo ""
echo "### $i Sorting by viral sequences list ###"
seqkit grep -f $title/$out/${i}_kaiju_classified_list -v $dir/${i}_R1_merged_unique.fastq.gz -o $title/${unclassified}/${i}_R1_unclassified.fastq.gz
seqkit grep -f $title/$out/${i}_kaiju_classified_list -v $dir/${i}_R2_merged_unique.fastq.gz -o $title/${unclassified}/${i}_R2_unclassified.fastq.gz
# extract by list (align)
seqkit grep -f $title/$out/${i}_kaiju_classified_list $dir/${i}_R1_merged_unique.fastq.gz -o $title/${classified}/${i}_R1_classified.fastq.gz
seqkit grep -f $title/$out/${i}_kaiju_classified_list $dir/${i}_R2_merged_unique.fastq.gz -o $title/${classified}/${i}_R2_classified.fastq.gz
echo "# Create krona file of $i"
kaiju2krona -t ${NODE} -n ${NAME} -i $title/$out/${i}_kaiju.txt -o $title/$krona/${i}_kaiju.krona
echo "# Create krona plot of $i"
ktImportText -o $title/$krona/${i}_kaiju.html $title/$krona/${i}_kaiju.krona
echo "# Create mpa file of ${i%_merged.fastq}_kaiju.txt"
kraken-mpa-report --db ${minikraken1} $title/$out/${i}_kaiju.txt > $title/$mpa/${i}.mpa.txt
done

echo "# Merge mpa files"
merge_metaphlan_tables.py kaiju/mpa/*.mpa.txt > kaiju/kaiju_mpa.txt
