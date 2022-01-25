
#!/bin/bash
# humann2 table normalization
# on the iMAC
# 20200805
# Loading mapping file from: /Users/anne-catrinuhlemann/miniconda2/envs/shotgun/lib/python2.7/site-packages/humann2/data/misc/map_ko_uniref90.txt.gz
# groups="uniref90_go uniref90_ko uniref90_eggnog uniref90_pfam uniref90_level4ec uniref90_infogo1000 uniref90_rxn"
# names="go kegg-orthology eggnog metacyc-rxn pfam ec infogo1000 kegg-module kegg-pathway metacyc-pwy uniref90"

source ~/opt/miniconda3/etc/profile.d/conda.sh
if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi


conda activate shotgun


# set input dir for seq.fasta

i=$(basename $1)

##### Gene Ontology
# dir path
out1="gene_ontology"
stratified_out1="gene_ontology/stratified_out"
rm -r $out1
mkdir -p $out1 $stratified_out1


groups="uniref90_go"
names="go"
set -- $names
echo $groups

for group in ${groups[*]}
do
# normalization cpm (counts per million)
humann2_renorm_table --input $i --output $out1/merged_genefamilies_cpm.txt --units cpm --update-snames;

# call kegg-orthology number
humann2_regroup_table --input $out1/merged_genefamilies_cpm.txt --output $out1/merged_genefamilies_${group}_cpm.txt --groups ${group};
# add names based on kegg-orthology number
humann2_rename_table --input $out1/merged_genefamilies_${group}_cpm.txt  --output $out1/merged_genefamilies_${group}_$1_cpm.txt --names $1;
# stratified_table
humann2_split_stratified_table --input $out1/merged_genefamilies_${group}_$1_cpm.txt --output $stratified_out1;
head -n 1 $stratified_out1/merged_genefamilies_${group}_$1_cpm_unstratified.txt >> $stratified_out1/merged_genefamilies_${group}_$1_cpm_unstratified_filtered.txt;
grep ":" $stratified_out1/merged_genefamilies_${group}_$1_cpm_unstratified.txt >> $stratified_out1/merged_genefamilies_${group}_$1_cpm_unstratified_filtered.txt;
shift
done









### for kegg-orthology
groups="uniref90_ko"
names="kegg-orthology"
echo $groups
# dir path
out2="kegg-orthology"
stratified_out2="kegg-orthology/stratified_out"
rm -r $out2
mkdir -p $out2 $stratified_out2

# normalization
humann2_renorm_table --input $i --output $out2/merged_genefamilies_cpm.txt --units cpm --update-snames;
# call kegg-orthology number
humann2_regroup_table --input $out2/merged_genefamilies_cpm.txt --output $out2/merged_genefamilies_${groups}_cpm.txt --groups ${groups};
# only for kegg orthology
conda activate fishtaco;
run_musicc.py $out2/merged_genefamilies_${groups}_cpm.txt  -n -c use_generic -v -o $out2/merged_genefamilies_KEGG_MUSiCC_corrected.txt; # learn_model 0 이 있으면 안된다.

# add names based on kegg-orthology number
conda activate shotgun;
humann2_rename_table --input $out2/merged_genefamilies_KEGG_MUSiCC_corrected.txt  --output $out2/merged_genefamilies_${groups}_${names}_cpm.txt --names ${names};
# stratified_table
humann2_split_stratified_table --input $out2/merged_genefamilies_${groups}_${names}_cpm.txt --output $stratified_out2;
head -n 1 $stratified_out2/merged_genefamilies_${groups}_${names}_cpm_unstratified.txt >> $stratified_out2/merged_genefamilies_${groups}_${names}_cpm_unstratified_filtered.txt;
grep ":" $stratified_out2/merged_genefamilies_${groups}_${names}_cpm_unstratified.txt >> $stratified_out2/merged_genefamilies_${groups}_${names}_cpm_unstratified_filtered.txt;


conda deactivate

### for fishtaco install problem
# sudo conda install -n fishtaco fishtaco
# sudo conda install -n fishtaco scikit-learn==0.18.1 -c conda-forge
