#!/usr/bin/env bash

python3 ./scripts/generate.py ./example/reference_chr1_chr2.fasta \
        --min_distance 200000 \
        --output ./example/query.fasta \
        --csv ./example/query.csv \
        --trl_count 20 \
        --trl_min_length 100000 \
        --trl_max_length 1000000 \
        --dup_count 20 \
        --dup_min_length 100000 \
        --dup_max_length 1000000 \
        --inv_count 20 \
        --inv_min_length 100000 \
        --inv_max_length 1000000 \
        --invtrl_count 20 \
        --invtrl_min_length 100000 \
        --invtrl_max_length 1000000 \
        --trl_interchr 20 \
        --trl_interchr_inv 20 \
        --ins_count 20 \
        --ins_min_length 10000 \
        --ins_max_length 100000 \
        --del_count 20 \
        --del_min_length 10000 \
        --del_max_length 100000

python3 ./scripts/check.py ./example/reference_chr1_chr2.fasta \
                                   ./example/query.fasta \
                                   ./example/query.csv

LC_ALL=C java -jar ./OMTools/OMTools.jar FastaToOM \
    --fastain ./example/reference_chr1_chr2.fasta \
    --enzyme BspQI \
    --refmapout ./example/reference.cmap

LC_ALL=C java -jar ./OMTools/OMTools.jar FastaToOM \
    --fastain ./example/query.fasta \
    --enzyme BspQI \
    --refmapout ./example/query.cmap

LC_ALL=C java -jar ./OMTools/OMTools.jar OptMapDataGenerator \
    --refmapin ./example/query.cmap \
    --optmapout ./example/query.sdata \
    --moleno 3000 \
    --scalesd 0.001 \
    --subound 1.001 \
    --slbound 0.999 \
    --fsize 1000000 \
    --fubound 1100000 \
    --flbound 900000 \
    --rsln 600 \
    --meas 500 \
    --fpr 0.000002

python3 ./scripts/overlap_sdata.py ./example/query.sdata ./example/query.csv \
    --out_overlap ./example/query_overlap.sdata \
    --out_nonoverlap ./example/query_non_overlap.sdata \
    --connection_out ./example/query_connection.csv

LC_ALL=C java -jar ./OMTools/OMTools.jar DataTools \
    --optmapin ./example/query_overlap.sdata \
    --optmapout ./example/query_overlap.cmap

cd coma
source coma/bin/activate
coma -r ../example/reference.cmap -q ../example/query_overlap.cmap \
-o ../results/alignment.xmap -oM split-alignments -sm ../results/out.smap -ms 9000 -er 5000

cd ..
source ./bionano/tools/pipeline/Solve3.8.2_main-04-30-24-rel/bionano_packages/VCFConverter/src/SMAP/bin/activate

./bionano/tools/pipeline/Solve3.8.2_main-04-30-24-rel/SMSV/20240506/variant_clustering \
--smapFile ./results/out.smap \
--xmapFile ./results/alignment.xmap \
--rCmapFile ./example/reference.cmap \
--out ./results \
--sample 1 \
--minNumMembers 2 \
--maxTransClusterDistance 50000 \
--no-confidence

python3 ./bionano/tools/pipeline/Solve3.8.2_main-04-30-24-rel/bionano_packages/VCFConverter/src/bionano_vcf_converter.py \
-s ./results/out.smap \
-r ./example/reference.cmap \
-x ./results/alignment.xmap --species_reference human_hg38 \
-C ./results/1_cluster_molecule_variant.txt