conda activate rnafold
RNAfold -i ../cleaned_nucleo_narna.fasta -T 37.5 > second_struc.txt
RNAdistance -Xm < second_struc.txt > dist_matrix.txt
