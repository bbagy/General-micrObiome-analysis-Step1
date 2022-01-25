#!/bin/bash
# Aug 22 2019
# 어려운 건 아니지만 쉽게 쓰기 위해 만들어 보았다.
# sudo conda install -n fishtaco scikit-learn==0.18.1 -c conda-forge
# source ~/miniconda2/etc/profile.d/conda.sh
source ~/opt/miniconda3/etc/profile.d/conda.sh

conda activate fishtaco

echo "#---- Creat DA_functions.tab ----#"
python ~/Dropbox/04_Scripts/R_source/fishtaco/fishtaco/compute_differential_abundance.py --class label.txt -o DA_functions.tab function.txt

echo "#--------  Run fishtaco --------#"
rm -r fishtaco
mkdir fishtaco
run_fishtaco.py -ta taxa.txt -fu function.txt -l label.txt -en DA_functions.tab  -inf --output_prefix fishtaco/fishtaco -log --multiple_hypothesis_correction none -map_function_level none -functional_profile_already_corrected_with_musicc


echo "#--------  Done --------#"

conda deactivate

# from the tutorial
# run_fishtaco.py -ta bv_qpcr.txt -fu bv_metagenomes.txt -l bv_metadata_fishtaco.txt -inf -log -assessment single_taxa
