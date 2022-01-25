#!/bin/bash
#-------------------------------#
#-- get certain taxa sequence --#
#-------------------------------#
# kaiju 부터 kraken mpa까지
# 200624
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# GoGet_taxa_V1.sh  --fasta=classified --pattern== --krakenFiles=kraken_out --taxaID== --out=Salmonella --paired=1(2)
# 200510
# Work station으로 DB loacation 이 설정.
# looping 기존에는 migered_uniq 는 kraken2 output을 이용하기 위해 _classified.fastq.gz  로 설정 되어 있음.
# #echo $(echo $i | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/');





fasta=""
krakenFiles=""
out=""

for arg in "$@"
do
    case $arg in
        -o=*|--out=*)
        out="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -f=*|--fasta=*)
        fasta="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -a=*|--pattern=*)
        pattern="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -k=*|--krakenFiles=*)
        krakenFiles="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -t=*|--taxaID=*)
        taxaID="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -p=*|--paired=*)
        paired="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` --fasta=[location of fastq dir] --krakenFiles=[location of kraken_out dir] --taxaID=28901 --out=[out dir] --paired=1(2)"
exit 0
fi


# dir path
out=$out/$taxaID
tem=$out/tem

#--- path set ---#
rm -r $out 
mkdir -p $tem $out
# set input dir for fastq or fasta


#dir=$(basename $1)
#threads=$(basename $2)

filenames=$(ls $fasta/*_R1_*.gz | xargs -n1 -I{} basename "{}" )
echo $filenames


if [ "$paired" -eq 1 ]; then
    echo ""
    echo "#===   Run for single end   ===#"
    echo ""
    for filename in $filenames
    do
        i=$(echo $filename | sed 's/\(.*\)\'"$pattern"'/\1/')
    
    echo ${i}
        echo ""
        echo "### $i Gathering target taxaID ###"
        echo ""
        awk '$3 == '"$taxaID"'' $krakenFiles/${i}_out.txt > $out/tem/${i}_taxaID.txt
        cut -f2 $tem/${i}_taxaID.txt > $tem/${i}_taxaID_list.txt

        seqkit grep -f  $out/tem/${i}_taxaID_list.txt $fasta/${i}$pattern -o $out/${i}_$taxaID.fastq.gz
    done
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
        echo "### $i Gathering target taxaID ###"
        echo ""
        awk '$3 == '"$taxaID"'' $krakenFiles/${i}_out.txt > $out/tem/${i}_taxaID.txt
        cut -f2 $tem/${i}_taxaID.txt > $tem/${i}_taxaID_list.txt

        # extract by list (align)
        seqkit grep -f $tem/${i}_taxaID_list.txt $fasta/${i}_R1$pattern  -o ${classified}/${i}_R1_$taxaID.fastq.gz
        seqkit grep -f $tem/${i}_taxaID_list.txt $fasta/${i}_R2$pattern  -o ${classified}/${i}_R2_$taxaID.fastq.gz
    done
else
    echo "Add paired information 1 or 2"
fi
