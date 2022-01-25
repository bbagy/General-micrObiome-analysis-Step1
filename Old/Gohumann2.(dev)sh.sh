#!/bin/bash
# 190822
# GoHumann2.sh [input dir name]
# requrement:
# conda install -c bioconda humann2
# sudo conda install -c conda-forge glpk


dir=$(basename $1)
rm -r humann2_out
filenames=$(ls $dir/* | xargs -n1 -I{} basename "{}" )

#echo $filenames

chocophlan=(/Users/heekukpark/DB/humann_db/db/chocophlan)
uniref=(/Users/heekukpark/DB/humann_db/db/db/uniref)


# run humann2
for i in $filenames
do
humann2 --input ${i} --threads 4 --output humann2_out  --metaphlan-options='--mpa_pkl /Users/heekukpark/DB/mpa_v20_m200/mpa_v20_m200.pkl --bowtie2db /Users/heekukpark/DB/mpa_v20_m200' --nucleotide-database ${chocophlan} --protein-database ${uniref}
done


# merge file
mkdir mkdir humann2_finalout
filenames=(genefamilies pathcoverage pathabundance)
for filename in ${filenames[*]}
do
humann2_join_tables -i humann2_out/ -o humann2_finalout/merged_${filename}.txt --file_name ${filename}
done


# for fish taco (merged_genefamilies_kegg.txt)
for i in $(ls humann2_finalout/*.txt)
do
humann2_regroup_table --input ${i} --output ${i%.txt}_kegg.txt --groups uniref50_rxn --ungrouped N --protected N
done


# merge metaphlan2_out
mkdir -p metaphlan2_out/raw
cp humann2_out/*/*metaphlan_bugs_list.tsv metaphlan2_out/raw
merge_metaphlan_tables.py metaphlan2_out/raw/*metaphlan_bugs_list.tsv > metaphlan2_out/metaphlan2_merged.txt
