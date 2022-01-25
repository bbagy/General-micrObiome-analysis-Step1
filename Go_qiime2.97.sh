#!/bin/bash
#-------------------------------#
#-- qiime2 non stop pipeline --#
#-------------------------------#
# qiime2- dada2-vsearch
# 20200908
# dir을 변수로 설정 했다면 실제 명령어에 $dir/${i} 이런 변수를 넣어야 한다.
# Go_qiime2.sh  --project=[] --fastq=[dir] --pattern=[_001.fastq.gz] --samNum=[50] --db=[db dir]  --paired=[1/2]
# require conda install -c anaconda scikit-learn -n qiime2
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2

fastq=""
DB=""


for arg in "$@"
do
    case $arg in
        -p=*|--project=*)
        project="${arg#*=}"
        shift # project name
        ;;
        -f=*|--fastq=*)
        fastq="${arg#*=}"
        shift # fastq dir
        ;;
        -a=*|--pattern=*)
        pattern="${arg#*=}"
        shift # files pattern
        ;;
        -n=*|--samNum=*)
        samNum="${arg#*=}"
        shift # number of samples
        ;;
        -d=*|--db=*)
        db="${arg#*=}"
        shift # Remove --initialize from processing
        ;;
        -p=*|--paired=*)
        paired="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
    esac
done


if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` --project==[] --fastq=[dir] --pattern=[_001.fastq.gz] --fileNum=[50] --db=[db dir]  --paired=[1/2]"
exit 0
fi


# dir path
out=qiime2_$project
qiime2=$out/qiime2_out
vsearch=$out/vsearch

#--- path set ---#
mkdir -p $out $qiime2 $vsearch

    echo "creating manifast file"
    # sample-id
    filenames=$(ls $fastq/*.gz | xargs -n1 -I{} basename "{}" )
    for filename in $filenames
    do
    echo $(echo $filename | sed 's/\(.*\)_L001_R[12]\'"$pattern"'/\1/') >> $out/sample-id_tem
    done

    echo -e "sample-id" | cat - $out/sample-id_tem > $out/sample-id

    #absolute-filepath
    ls -d $PWD/$fastq/* > $out/absolute-filepath_tem
    echo -e "absolute-filepath" | cat - $out/absolute-filepath_tem > $out/absolute-filepath
    
    # direction
if [ "$paired" -eq 1 ]; then

    for i in `seq $samNum`
    do
    echo "forward" >> $out/direction_tem
    done
    echo -e "direction" | cat - $out/direction_tem > $out/direction
    rm $out/sample-id_tem $out/absolute-filepath_tem $out/direction_tem

elif [ "$paired" -eq 2 ]; then
    for i in `seq $samNum`
    do
    echo "forward" >> $out/direction_tem
    echo "reverse" >> $out/direction_tem
    done
    echo -e "direction" | cat - $out/direction_tem > $out/direction
    rm $out/sample-id_tem $out/absolute-filepath_tem $out/direction_tem
else
    echo "Add paired information 1 or 2"
fi

    # merge column
    pr -mts, $out/sample-id $out/absolute-filepath $out/direction > $out/manifast
    rm $out/sample-id $out/absolute-filepath $out/direction
    
    ### run qiime2
    # step 1 import
    
if [ "$paired" -eq 1 ]; then
   qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path $out/manifast --output-path $qiime2/1_${project}.demux.qza \
    --input-format SingleEndFastqManifestPhred33

elif [ "$paired" -eq 2 ]; then
    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $out/manifast --output-path $qiime2/1_${project}.demux.qza \
    --input-format PairedEndFastqManifestPhred33
else
    echo "Add paired information 1 or 2"
fi
    
    

    
    # step 2 denoise
    
if [ "$paired" -eq 1 ]; then

    qiime dada2 denoise-single \
    --i-demultiplexed-seqs $qiime2/1_${project}.demux.qza \
    --p-trim-left 5 \
    --p-n-threads 0 \
    --p-trunc-len 200 \
    --o-representative-sequences $qiime2/2_${project}_rep-seqs-dada2.qza \
    --o-table $qiime2/3_${project}_table-dada2.qza \
    --o-denoising-stats $qiime2/4_${project}_denoising-dada2.qza \
    --verbose

elif [ "$paired" -eq 2 ]; then
    qiime dada2 denoise-paired \
    --i-demultiplexed-seqs $qiime2/1_${project}.demux.qza \
    --p-trim-left-f 5 \
    --p-trim-left-r 5 \
    --p-n-threads 0 \
    --p-trunc-len-f 250 \
    --p-trunc-len-r 250 \
    --o-representative-sequences $qiime2/2_${project}_rep-seqs-dada2.qza \
    --o-table $qiime2/3_${project}_table-dada2.qza \
    --o-denoising-stats $qiime2/4_${project}_denoising-dada2.qza \
    --verbose
else
    echo "Add paired information 1 or 2"
fi

    qiime feature-table tabulate-seqs \
    --i-data $qiime2/2_${project}_rep-seqs-dada2.qza \
    --o-visualization $qiime2/2_${project}_rep-seqs-dada2.qzv
    
     # step 3 classification
    qiime feature-classifier classify-sklearn \
    --i-classifier $db \
    --i-reads $qiime2/2_${project}_rep-seqs-dada2.qza \
    --o-classification $qiime2/5_${project}_rep-seqs-dada2_classified.qza
    
    # step 4 Make Phylogenetic Trees
    qiime alignment mafft \
      --i-sequences $qiime2/2_${project}_rep-seqs-dada2.qza \
      --o-alignment $qiime2/6_${project}_aligned_rep-seqs.qza
    qiime alignment mask \
      --i-alignment $qiime2/6_${project}_aligned_rep-seqs.qza \
      --o-masked-alignment $qiime2/7_${project}_aligned_masked_rep-seqs.qza
    qiime phylogeny fasttree \
      --i-alignment $qiime2/7_${project}_aligned_masked_rep-seqs.qza \
      --o-tree $qiime2/8_${project}_unrooted-tree.qza
    qiime phylogeny midpoint-root \
      --i-tree $qiime2/8_${project}_unrooted-tree.qza \
      --o-rooted-tree $qiime2/9_${project}_unrooted-tree.qza
      
    ## step 5 exporting
    # tree
    qiime tools export \
    --input-path $qiime2/9_${project}_unrooted-tree.qza \
    --output-path $qiime2/10_${project}_exported
    # taxa table
    qiime tools export \
    --input-path $qiime2/5_${project}_rep-seqs-dada2_classified.qza \
    --output-path $qiime2/10_${project}_exported
    
    # count table
    qiime tools export \
    --input-path $qiime2/3_${project}_table-dada2.qza \
    --output-path $qiime2/10_${project}_exported

    biom convert -i $qiime2/10_${project}_exported/feature-table.biom -o $qiime2/10_${project}_exported/otutable.txt --to-tsv

    sed '1d' $qiime2/10_${project}_exported/otutable.txt > $qiime2/10_${project}_exported/otutable_cleaned.txt

    ## step 6 qiime2 to phyloseq

    Go_qiime2ps.R -p ${project} -b $qiime2/10_${project}_exported/otutable_cleaned.txt -x $qiime2/10_${project}_exported/taxonomy.tsv -t $qiime2/10_${project}_exported/tree.nwk -o $qiime2/10_${project}_exported
    conda deactivate

    
#====== vsearch after dsds2 =========#
    # step 3 vsearch
    conda activate qiime2
    qiime vsearch cluster-features-de-novo \
    --i-sequences $qiime2/2_${project}_rep-seqs-dada2.qza \
    --i-table $qiime2/3_${project}_table-dada2.qza \
    --p-perc-identity .97 --p-threads 6 --output-dir vsearch \
    --o-clustered-table $vsearch/2_${project}_vsearch_clustered_tab.qza \
    --o-clustered-sequences $vsearch/1_${project}_vsearch_clustered_seq.qza \
    --verbose
    rm -r vsearch

    # step 4 classification
    qiime feature-classifier classify-sklearn \
    --i-classifier $db \
    --i-reads $vsearch/1_${project}_vsearch_clustered_seq.qza \
    --o-classification $vsearch/3_${project}_vsearch_clustered_classified.qza
    
    # step 5 Make Phylogenetic Trees
    qiime alignment mafft \
      --i-sequences $vsearch/1_${project}_vsearch_clustered_seq.qza \
      --o-alignment $vsearch/4_${project}_vsearch_clustered_aligned_rep-seqs.qza

    qiime alignment mask \
      --i-alignment $vsearch/4_${project}_vsearch_clustered_aligned_rep-seqs.qza \
      --o-masked-alignment $vsearch/5_${project}_vsearch_clustered_aligned_masked_rep-seqs.qza

    qiime phylogeny fasttree \
      --i-alignment $vsearch/5_${project}_vsearch_clustered_aligned_masked_rep-seqs.qza \
      --o-tree $vsearch/6_${project}_vsearch_clustered_aligned_masked_unrooted-tree.qza

    qiime phylogeny midpoint-root \
      --i-tree $vsearch/6_${project}_vsearch_clustered_aligned_masked_unrooted-tree.qza \
      --o-rooted-tree $vsearch/7_${project}_vsearch_clustered_aligned_masked_rooted-tree.qza
    
    ## step 6 exporting
    # tree
    qiime tools export \
    --input-path $vsearch/7_${project}_vsearch_clustered_aligned_masked_rooted-tree.qza \
    --output-path $vsearch/8_${project}_vsearch_exported
    # taxa table
    qiime tools export \
    --input-path $vsearch/3_${project}_vsearch_clustered_classified.qza \
    --output-path $vsearch/8_${project}_vsearch_exported
    
    # count table
    qiime tools export \
    --input-path $vsearch/2_${project}_vsearch_clustered_tab.qza \
    --output-path $vsearch/8_${project}_vsearch_exported

    biom convert -i $vsearch/8_${project}_vsearch_exported/feature-table.biom -o $vsearch/8_${project}_vsearch_exported/otutable.txt --to-tsv

    sed '1d' $vsearch/8_${project}_vsearch_exported/otutable.txt > $vsearch/8_${project}_vsearch_exported/otutable_cleaned.txt

    ## step 7 qiime2 to phyloseq
    
    Go_qiime2ps.R -p ${project} -b $vsearch/8_${project}_vsearch_exported/otutable_cleaned.txt -x $vsearch/8_${project}_vsearch_exported/taxonomy.tsv -t $vsearch/8_${project}_vsearch_exported/tree.nwk -o $vsearch/8_${project}_vsearch_exported
    conda deactivate
