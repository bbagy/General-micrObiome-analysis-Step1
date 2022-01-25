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


if [ "$paired" -eq 1 ]; then
    echo "single"


# paired end
elif [ "$paired" -eq 2 ]; then

    echo "paired"

#grep -c "U  N" kaiju/kaiju_out/* > kaiju/reads_unclassified.txt
#grep -c "C    N" kaiju/kaiju_out/* > kaiju/reads_classified.txt
    
else
    echo "Add paired information 1 or 2"
fi
