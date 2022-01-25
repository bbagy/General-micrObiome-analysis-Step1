#!/bin/bash
# 190905
#---------------#
#-- run kaiju --#
#---------------#
# kaiju 부터 kraken mpa까지
# 190905
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# GoKaiju.sh [input dir name]


if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi


source activate shotgun

# set input dir for fastq or fasta

dir=$(basename $1)
threads=$(basename $2)

filenames=$(ls $dir/* | xargs -n1 -I{} basename "{}" )

/Users/heekukpark/DB/kaijuDB/progenomes/kaiju_db_progenomes.fmi
# location of DB
DB=(~/DB/kaijuDB/progenomes/kaiju_db_progenomes.fmi)
NODE=(~/DB/kaijuDB/nodes.dmp)
NAME=(~/DB/kaijuDB/names.dmp)
minikraken1=(~/DB/kraken_db/minikraken1)


rm -r kaiju
mkdir -p kaiju/kaiju_out kaiju/mpa kaiju/krona


for filename in $filenames
do

i=$(echo $filename)

echo "# Run kaiju of $i"
kaiju -t ${NODE} -f ${DB} -a greedy -e 3 -E 0.05 -i $dir/${i} -o kaiju/kaiju_out/${i%_merged.fastq}_kaiju.txt -v -z $threads

echo "# Create krona file of $i"
kaiju2krona -t ${NODE} -n ${NAME} -i kaiju/kaiju_out/${i%_merged.fastq}_kaiju.txt -o kaiju/krona/${i%_merged.fastq}_kaiju.krona

echo "# Create krona plot of $i"
ktImportText -o kaiju/krona/${i%_merged.fastq}_kaiju.html kaiju/krona/${i%_merged.fastq}_kaiju.krona


echo "# Create mpa file of ${i%_merged.fastq}_kaiju.txt"
kraken-mpa-report --db ${minikraken1} kaiju/kaiju_out/${i%_merged.fastq}_kaiju.txt > kaiju/mpa/${i%_merged.fastq}.mpa.txt
done

echo "# Merge mpa files"
merge_metaphlan_tables.py kaiju/mpa/*.mpa.txt > kaiju/kaiju_mpa.txt



#grep -c "U  N" kaiju/kaiju_out/* > kaiju/reads_unclassified.txt
#grep -c "C    N" kaiju/kaiju_out/* > kaiju/reads_classified.txt




echo "#-------- done -------#"

conda deactivate





