#!/usr/bin/bash
eval "$(conda shell.bash hook)"

conda activate R_env
mkdir pdf_plots
touch pdf_plots/create_pdf.log
cat 'NEW RUN!' >> pdf_plots/create_pdf.log
date >> pdf_plots/create_pdf.log
for dir in bakta_annotation_*/; do
    if [ -f "${dir}"*.ffn ]; then
        parent_dir=$(basename "$dir")
        clean_name=${parent_dir#bakta_annotation_}
        Rscript try_find_genes.R ${parent_dir}/${clean_name}.tsv pdf_plots/${clean_name}.pdf >> pdf_plots/create_pdf.log 2>&1
    fi
done