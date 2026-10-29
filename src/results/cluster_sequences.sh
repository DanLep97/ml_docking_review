# Create separate databases for train and test
mmseqs createdb ../../databases/PDBbindv2020/proteins.fasta trainDB
mmseqs createdb ../../databases/QBioLiP/enzymes.fasta testDB

# Search test sequences against train sequences
mmseqs search testDB trainDB test_vs_trainDB tmp_folder \
    --search-type 3 \
    -s 7.5 \
    -e 1 \
    --min-seq-id 0 \
    -c 0.0 \
    --cov-mode 0 \
    --threads 8 \
    --max-seqs 20000000 \

# Convert alignment results to tabular format
mmseqs convertalis testDB trainDB test_vs_trainDB test_vs_train.tsv \
  --format-output "query,target,pident,alnlen,evalue,bits"

rm test_vs_trainDB*
rm testDB*
rm trainDB*