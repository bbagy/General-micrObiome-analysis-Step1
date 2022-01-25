
#!/bin/bash
# DADA2 out file을 이용한 tree file을 만들기 위한 qiime2 명령어 종합
#source ~/miniconda2/etc/profile.d/conda.sh

#source ~/miniconda3/etc/profile.d/conda.sh
source ~/opt/miniconda3/etc/profile.d/conda.sh

if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` [input dir name]"
exit 0
fi


conda activate qiime2


# set input dir for seq.fasta

i=$(basename $1)


mkdir -p ${i}_tree/
qiime tools import --input-path $i --output-path ${i}_tree/${i%.fna}.qza --type 'FeatureData[Sequence]'

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences ${i}_tree/${i%.fna}.qza   \
--o-alignment ${i}_tree/${i%.fna}_aligned-rep-seqs.qza \
--o-masked-alignment ${i}_tree/${i%.fna}_masked-aligned-rep-seqs.qza \
--o-tree ${i}_tree/${i%.fna}_unrooted-tree.qza   \
--o-rooted-tree ${i}_tree/${i%.fna}_rooted-tree.qza

qiime tools export --input-path ${i}_tree/${i%.fna}_rooted-tree.qza --output-path ${i}_tree/exported-tree



conda deactivate





#for i in *.fasta
#do
#mkdir -p ${i}_tree/
#qiime tools import --input-path $i --output-path ${i}_tree/${i%.fasta}.qza --type 'FeatureData[Sequence]'

#qiime phylogeny align-to-tree-mafft-fasttree \
#--i-sequences ${i}_tree/${i%.fasta}.qza   \
#--o-alignment ${i}_tree/${i%.fasta}_aligned-rep-seqs.qza \
#--o-masked-alignment ${i}_tree/${i%.fasta}_masked-aligned-rep-seqs.qza \
#--o-tree ${i}_tree/${i%.fasta}_unrooted-tree.qza   \
#--o-rooted-tree ${i}_tree/${i%.fasta}_rooted-tree.qza

#qiime tools export --input-path ${i}_tree/${i%.fasta}_rooted-tree.qza --output-path exported-tree
#done
