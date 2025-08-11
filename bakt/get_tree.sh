#!/usr/bin/bash
eval "$(conda shell.bash hook)"
conda activate rna_fold

touch naRNAfastas/all_naRNA.fasta
for dir in bakta_annotation_*/; do
    if [ -f "${dir}"*.ffn ]; then
        parent_dir=$(basename "$dir")
        grep -A 1 'Nucleoid-associated noncoding RNA 4' "${dir}"*.ffn | grep -v '^--$' > "naRNAfastas/${parent_dir}_naRNA.fasta"
        cat naRNAfastas/${parent_dir}_naRNA.fasta >> naRNAfastas/all_naRNA.fasta
        RNAfold -i naRNAfastas/${parent_dir}_naRNA.fasta -T 37.5 --noPS > "naRNAfastas/${parent_dir}_second_structure.txt"
        RNAdistance -Xm < naRNAfastas/${parent_dir}_second_structure.txt > "naRNAfastas/${parent_dir}_dist.mat"
    fi
done

conda activate R_env

mkdir hc_and_heatmap
touch hc_and_heatmap/create_tree.log
for dismat in naRNAfastas/*.mat; do
    matrix_name_dir=$(basename "$dismat")
    Rscript dist_trees_script.R ${dismat} ${matrix_name_dir} >> hc_and_heatmap/create_tree.log 2>&1
done

zip hc_and_heatmap/hc_heat_plots.zip hc_and_heatmap/*.pdf