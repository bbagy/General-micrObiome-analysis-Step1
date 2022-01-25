#!/bin/bash
#-----------------#
#-- run metahit --#
#-----------------#
# megahit runniung command
# 20201014
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# Go_megahit.sh -t=[title] -f=[files path] -p=1
# co-assemble 이 가능 하도록 만들었다.

source ~/miniconda3/etc/profile.d/conda.sh
conda activate MAGs


title=""
dir=""
mincontiglen=200
numthreads=20
memory=0.85

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
        -p=*|--paired=*)
        paired="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` -t=[title] -f=[files path] -p=1(single end)/2(paired end)"
exit 0
fi


# dir path
out=${title}_assemple


#--- path set ---#
rm -r $out
#mkdir -p $out
# set input dir for fastq or fasta

if [ "$paired" -eq 1 ]; then
    echo ""
    echo "#===   Run for single end   ===#"
    echo ""
    R1s=`ls -d $dir/*.gz | ls $dir/ -p | grep -v / | tr '\n' ','`
    echo $dir
    echo $R1s
    echo ""
    echo "megahit -r $R1s --min-contig-len $mincontiglen -m $memory -o $out/ -t $numthreads --k-min 63 --k-max 91 --k-step 28"
    megahit -r $R1s --min-contig-len $mincontiglen -m $memory -o $out -t $numthreads --k-min 63 --k-max 91 --k-step 28

# paired end
elif [ "$paired" -eq 2 ]; then
    echo ""
    echo "#===   Run for paired end   ===#"
    echo ""
    R1s=`ls -d $dir/*R1*.gz | ls $dir/ -p | grep -v / | tr '\n' ','`
    R2s=`ls -d $dir/*R2*.gz | ls $dir/ -p | grep -v / | tr '\n' ','`
    echo ""
    echo $R1s
    echo ""
    echo $R2s
    echo ""
    echo "megahit -r $R1s --min-contig-len $mincontiglen -m $memory -o $out/ -t $numthreads --k-min 63 --k-max 91 --k-step 28"
    megahit -1 $R1s -2 $R2s --min-contig-len $mincontiglen -m $memory -o $out/ -t $numthreads --k-min 63 --k-max 91 --k-step 28
    
else
    echo "Add paired information 1 or 2"
fi
