#!/bin/bash
#----------------#
#-- get genome --#
#----------------#
# 190905
# https://www.ncbi.nlm.nih.gov/books/NBK179288/
# source ./install-edirect.sh

if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input file name] [input file name]"
exit 0
fi


# set input file and out file

input=$(basename $1)
output=$(basename $2)
#filenames=$(ls $dir/* | xargs -n1 -I{} basename "{}" )


samples=(`awk '{ print $1}' $input`)

for org in ${samples[*]}
do
echo $org
esearch -db assembly -query "$org [ORGN]" |
efilter -query "representative [PROP]" |
elink -target nuccore -name assembly_nuccore_refseq |
efetch -format docsum |
xtract -pattern DocumentSummary -element  AccessionVersion TaxId Completeness Slen Title |
sed 's/,.*//' |
grep -v -i -e scaffold -e plasmid -e sequence -e patch |
sort -t $'\t' -k 2,2nr >> $output
done


echo "#-------- done -------#"
