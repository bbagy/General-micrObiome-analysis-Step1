#!/bin/bash
# 190905
#---------------#
#-- run kaiju --#
#---------------#
# kaiju 부터 kraken mpa까지
# 190905
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# GoKaiju_vir_V3.sh -t=Virome_single -f=merged_uniq -c=4 -p=1
# 200510
# Work station으로 DB loacation 이 설정.
# looping 기존에는 migered_uniq 는 kraken2 output을 이용하기 위해 _classified.fastq.gz  로 설정 되어 있음.
# #echo $(echo $i | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/');

#GoKaiju_vir_V3.sh

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
out=$title/kaiju_out
mpa=$title/mpa
krona=$title/krona
classified=$title/classified
unclassified=$title/unclassified

#--- path set ---#
rm -r $title
mkdir -p $title $out $mpa $krona
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
DB=(/media/uhlemann/core1/DB/kaijuDB/viruses/kaiju_db_viruses.fmi)
#DB=(/media/uhlemann/core1/DB/kaijuDB/progenomes/kaiju_db_progenomes.fmi)

NODE=(/media/uhlemann/core1/DB/kaijuDB/nodes.dmp)
NAME=(/media/uhlemann/core1/DB/kaijuDB/names.dmp)
minikraken1=(/media/uhlemann/core1/DB/kaijuDB/minikraken1)


if [ "$paired" -eq 1 ]; then
    echo ""
    echo "#===   Run for single end   ===#"
    echo ""
    for filename in $filenames
    do

        i=$(echo $filename | sed 's/\(.*\)_R[1]\'"$pattern"'/\1/')
    
    echo ${i}
        echo ""
        echo "### $i Gathering total viral sequences ###"
        echo ""
        echo "kaiju -t node -f DB -a greedy -e 3 -E $evalue -i ${i}_R1$pattern  -o ${i}_kaiju.txt -v -z $cpu"
        kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E $evalue -i $dir/${i}_R1$pattern   -o ${out}/${i}_kaiju.txt -v -z $cpu
        echo ""
    
        kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i $out/${i}_kaiju.txt -o  $out/${i}_kaiju_classified.txt -u

#keep classified list # -wE를 하면 정확히 C 만 골라낼수 있다.
        grep -v "^@" $out/${i}_kaiju_classified.txt | cut -f2 > $out/${i}_kaiju_classified_list


    # remove by list (unalign)
        echo ""
        echo "### $i Sorting by viral sequences list ###"
        echo ""
        seqkit grep -f $out/${i}_kaiju_classified_list -v $dir/${i}_R1$pattern -o $unclassified/${i}_R1_unclassified.fastq.gz

    # extract by list (align)
        seqkit grep -f $out/${i}_kaiju_classified_list $dir/${i}_R1$pattern  -o $classified/${i}_R1_classified.fastq.gz
        echo ""
        echo "# Create krona file of $i"
        echo ""
        kaiju2krona -t ${NODE} -n ${NAME} -i $out/${i}_kaiju.txt -o $krona/${i}_kaiju.krona
        echo ""
        echo "# Create krona plot of $i"
        echo ""
        ktImportText -o $krona/${i}_kaiju.html $krona/${i}_kaiju.krona
        echo ""
        echo "# Create mpa file of ${i}_kaiju.txt"
        echo ""
        kraken-mpa-report --db $minikraken1 $out/${i}_kaiju.txt > $mpa/${i}.mpa.txt
    done

    echo "# Merge mpa files"
    merge_metaphlan_tables.py $mpa/*.mpa.txt > $title/kaiju_mpa.txt


# paired end
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
        echo "kaiju -t NODE -f DB -a greedy -e 3 -E $evalue-i ${i}_R1$pattern  -j ${i}_R2$pattern  -o ${i}_kaiju.txt -v -z $cpu"
        echo ""
        kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E $evalue -i $dir/${i}_R1$pattern  -j $dir/${i}_R2$pattern  -o ${out}/${i}_kaiju.txt -v -z $cpu
        kaiju-addTaxonNames -t ${NODE} -n ${NAME} -i $out/${i}_kaiju.txt -o  $out/${i}_kaiju_classified.txt -u
        #keep classified list
        grep -v "^@" $out/${i}_kaiju_classified.txt | cut -f2 > $out/${i}_kaiju_classified_list
        
        echo ""
        echo "### $i Sorting by viral sequences list ###"
        echo ""
        # remove by list (unalign)
        seqkit grep -f $out/${i}_kaiju_classified_list -v $dir/${i}_R1$pattern  -o ${unclassified}/${i}_R1_unclassified.fastq.gz
        seqkit grep -f $out/${i}_kaiju_classified_list -v $dir/${i}_R2$pattern  -o ${unclassified}/${i}_R2_unclassified.fastq.gz
        # extract by list (align)
        seqkit grep -f $out/${i}_kaiju_classified_list $dir/${i}_R1$pattern  -o ${classified}/${i}_R1_classified.fastq.gz
        seqkit grep -f $out/${i}_kaiju_classified_list $dir/${i}_R2$pattern  -o ${classified}/${i}_R2_classified.fastq.gz
        echo ""
        echo "# Create krona file of $i"
        echo ""
        kaiju2krona -t ${NODE} -n ${NAME} -i $out/${i}_kaiju.txt -o ${krona}/${i}_kaiju.krona
        echo ""
        echo "# Create krona plot of $i"
        echo ""
        ktImportText -o ${krona}/${i}_kaiju.html ${krona}/${i}_kaiju.krona
        echo ""
        echo "# Create mpa file of ${i%_merged.fastq}_kaiju.txt"
        echo ""
        kraken-mpa-report --db $minikraken1 $out/${i}_kaiju.txt > $mpa/${i}.mpa.txt
    done
    
    echo ""
    echo "# Merge mpa files"
    echo ""
    merge_metaphlan_tables.py $mpa/*.mpa.txt > $title/kaiju_mpa.txt
    
else
    echo "Add paired information 1 or 2"
fi
