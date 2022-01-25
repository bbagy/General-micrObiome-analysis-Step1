#!/bin/bash
#-------------------------------#
#-- get certain taxa sequence --#
#-------------------------------#
# fastq.gz로부터 특정 taxa 서열을 추출 한다.taxaID를 기본으로 추출이 가능 하다.
# 200624
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# GoGet_taxa.sh  --fasta=classified --pattern== --krakenFiles=kraken_out --taxaID== --out=Salmonella --paired=1(2)
# 200510
# Work station으로 DB loacation 이 설정.
# looping 기존에는 migered_uniq 는 kraken2 output을 이용하기 위해 _classified.fastq.gz  로 설정 되어 있음.
# #echo $(echo $i | sed 's/\(.*\)_R[12]_\merged_unique.fastq.gz/\1/');


source ~/miniconda3/etc/profile.d/conda.sh
conda activate shotgun


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
out=extract_$out
taxaIDdir=$out/$taxaID
tem=$out/tem

#--- path set ---#
#rm -r $out
mkdir -p $out $tem $taxaIDdir
# set input dir for fastq or fasta


#dir=$(basename $1)
#threads=$(basename $2)


#echo $filenames


if [ "$paired" -eq 1 ]; then
    filenames=$(ls $fasta/*.gz | xargs -n1 -I{} basename "{}" )
    echo ""
    echo "#===   Run for single end   ===#"
    echo ""
    for filename in $filenames
    do
        i=$(echo $filename | sed 's/\(.*\)\'"$pattern"'/\1/')
    
    echo ${i}
        echo ""
        echo "### $i Gathering target taxaID ###"
        echo "awk '$3 == '"$taxaID"'' $krakenFiles/${i}_out.txt > $tem/${i}_$taxaID.txt"
        awk '$3 == '"$taxaID"'' $krakenFiles/${i}_out.txt > $tem/${i}_${taxaID}_out.txt
        cut -f2 $tem/${i}_${taxaID}_out.txt > $tem/${i}_${taxaID}_list.txt

        seqkit grep -f  $tem/${i}_${taxaID}_list.txt $fasta/${i}$pattern -o $taxaIDdir/${i}_${taxaID}_fastq.gz
    done
# paired end
elif [ "$paired" -eq 2 ]; then
    filenames=$(ls $fasta/*_R1_*.gz | xargs -n1 -I{} basename "{}" )
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
        echo "awk '$3 == '"$taxaID"'' $krakenFiles/${i}_out.txt > $tem/${i}_$taxaID.txt"
        awk '$3 == '"$taxaID"'' $krakenFiles/${i}_out.txt > $tem/${i}_${taxaID}_out.txt
        cut -f2 $tem/${i}_${taxaID}_out.txt > $tem/${i}_${taxaID}_list.txt


        echo "### $i sorting ###"
        # extract by list (align)
        seqkit grep -f $tem/${i}_${taxaID}_list.txt $fasta/${i}_R1$pattern  -o ${taxaIDdir}/${i}_R1_${taxaID}_fastq.gz
        seqkit grep -f $tem/${i}_${taxaID}_list.txt $fasta/${i}_R2$pattern  -o ${taxaIDdir}/${i}_R2_${taxaID}_fastq.gz
    done
else
    echo "Add paired information 1 or 2"
fi
conda deactivate
