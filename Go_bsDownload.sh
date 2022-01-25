#!/bin/bash
# 20210408 / Heekuk Park
# Go_bsDownload.sh -m=masterDIR -i=IDs -n=names
## requirement
##  bs in path

mastarDIR=""
IDs=""
names=""

for arg in "$@"
do
    case $arg in
        -m=*|--mastarDIR=*)
        mastarDIR="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -i=*|--IDs=*)
        IDs="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -n=*|--names=*)
        names="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
    esac
done

if [ "$1" == "-h" ]; then
echo ""
echo "Usage: `basename $0` -m=[mastarDIR] -i=[IDs] -n=[names]"
echo ""
exit 0
fi

# dir path
cd $mastarDIR

#--- path set ---#
set -- $names
for ID in $IDs
do
bs download project -i $ID  -o ${mastarDIR}/$1 --extension=fastq.gz

rm ${mastarDIR}/$1/*.json

mkdir -p ${mastarDIR}/$1_fastq
cd ${mastarDIR}/$1
for subdir in *
do
cp $subdir/*.gz ${mastarDIR}/$1_fastq
done
cd $mastarDIR

shift
done
