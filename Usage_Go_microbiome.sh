# ============================= #
#   Useage for Go_microbiome    #
# ============================= #

### qc for nextseq output
# Goqc_nextera_fastq_V4.sh [fastq] [/media/uhlemann/core1/DB/human_index/GRCh38/GRCh38] 1
# 1/2 로 paired and single 을 지정할수 있다.
mmindex=(/media/uhlemann/core1/DB/mouse_index/mm10/mm10)
hmindex=(/media/uhlemann/core1/DB/human_index/GRCh38/GRCh38)
Goqc_nextera_fastq_V5.sh fastq_all $mmindex 1




### Soapdeveno2-63mer illumina assembly for single end (20201013)
# configFile is required, it is in "/home/uhlemann/heekuk_path/"
Go_soapdenovo2.sh -t=project -f=merged_uniq -a=pattern -c=/home/uhlemann/heekuk_path/configFile

Go_megahit.sh -t=Arpaia -f=./ -p=1 # present DIR에서만 작동

### bowtie2 mapping for single end (20201013)
# Go_mapping.sh -t=project -f=merged_uniq -a=pattern -g=contigs -b=contigs_pattern
Go_mapping.sh -t=Arpaia -f=filtered_homopolymer/ -a=_homo_filtered_merged_unique.fastq.gz -c=filtered_homopolymer/Arpaia_assemple/final.contigs.fa


Go_binning.sh -t=Arpaia -c=filtered_homopolymer/Arpaia_assemple/final.contigs.fa -b=Arpaia_mapping/bam -f=Arpaia_mapping/count -a=.counts


### run kraken2 (core2 linux)
Humankraken2=(/media/uhlemann/core1/DB/kraken2_human)
kraken2DB=(/media/uhlemann/core1/DB/kraken2_STDB)
# single
Gokraken2.pl -thread 16 -o kraken_out -b bracken_out -c classified -u unclassified -r kraken_report -log kraken2/kraken2_log.txt --db ${Humankraken2} merged_uniq/*.fastq.gz
# paired
Gokraken2.pl -thread 16 -o kraken_out -b bracken_out -c classified -u unclassified -r kraken_report -log kraken2/kraken2_log.txt --db ${Humankraken2} --paired filtered_homopolymer/*.fastq.gz



### sequences extract
# GoGet_taxa_V2.sh --fasta=[dir] --pattern=[files pattern] --krakenFiles=[kraken_out] --taxaID="28901" --out=[tatget name] --paired=1/2

GoGet_taxa_V2.sh --fasta=classified_out --pattern=_classified.fastq.gz --krakenFiles=kraken_out --taxaID="28901" --out=Salmonella --paired=2







### run fishtaco
# use this script in the same dir
Gofishtaco.sh


### dada2 for 16S microbiome
# Go_dada2.R -p [project] -f [fastq] -t 5 -d [db dir]
DB=(/media/uhlemann/core1/DB/DADA2/Bacteria/silva_nr99_v138_wSpecies_train_set.fa.gz)
Go_dada2.R -p ${project} -f fastq/ -t 5 -d $DB

### generation tree file for dada2
# qiime_tree.sh seqs.fna # require run on same dir with seqs.fna
qiime_tree.sh



conda activate MAGs

Go_megahit.sh -t=Arpaia -f=./ -p=1 # present DIR에서만 작동
Go_mapping.sh -t=Arpaia -f=filtered_homopolymer/ -a=_homo_filtered_merged_unique.fastq.gz -c=filtered_homopolymer/Arpaia_assemple/final.contigs.fa

Go_binning.sh -t=Arpaia -c=filtered_homopolymer/Arpaia_assemple/final.contigs.fa -b=Arpaia_mapping/bam -f=Arpaia_mapping/count -a=.counts
