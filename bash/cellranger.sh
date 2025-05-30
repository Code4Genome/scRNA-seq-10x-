#!/bin/bash
REFERENCE="/home/c/cellranger/refdata_genome"
FASTQ_DIR="/home/c/fastqs"
OUT_DIR="/home/c/Analysis/test"
SAMPLES=("SA" "SB" "SC" "SD", "SE")  # Add your sample names

# RUN CELLRANGER COUNT FOR EACH SAMPLE
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing ${SAMPLE}..."

    cellranger count \
        --id="${SAMPLE}" \
        --transcriptome="${REFERENCE}" \
        --fastqs="${FASTQ_DIR}/${SAMPLE}" \
        --sample="${SAMPLE}" \
        --localcores=8 \
        --localmem=64

    # Move results to output directory
    mv "${SAMPLE}" "${OUT_DIR}/${SAMPLE}"
done

    chmod +x run_cellranger.sh
    ./run_cellranger.sh
