for i in {1..22}; do
    CHIT_IN_FILE=./RNAseq/cht_input_file_atopy_ext.ranphe.txt.$i
    nohup python ./CHIT/CHIT.py \
        --bnb_disp ./RNAseq/cht_bnb_coef_atopy.txt \
        --as_disp ./RNAseq/cht_as_coef_atopy.txt \
        --cov_file ./pca/cov.matrix.txt \
        $CHIT_IN_FILE ./RNAseq/output/chit_atopy_results.txt.$i &
    sleep 1
done
