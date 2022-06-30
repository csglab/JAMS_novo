#!/bin/bash
which R
EXPERIMENT_ID=CTCF_HEK293_GSM2026781_smallest
DATA=./data/CTCF_demo
DATA_DIR=${DATA}/02_formatted_data/smallest_demo
OUT_DIR=${DATA}/05_motif_discovery/runs

set -o xtrace
./JAMS --task NOVO \
     --experiment ${EXPERIMENT_ID} \
     --flanking 20 \
     --novo_max_iterations 100 \
     --data_dir ${DATA_DIR} \
     --novo_start_motif_length 8 \
     --novo_skip_data \
     --novo_exclude_meth \
     --output_dir ${OUT_DIR}

# --novo_exclude_meth
# --novo_start_motif_length

exit
