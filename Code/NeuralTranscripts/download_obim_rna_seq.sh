#!/usr/bin/env bash

OUTPUT_DIR=O.bim/Data/PRJNA285380/Raw

mkdir -p $OUTPUT_DIR

accessions=$(tail -n +2 Code/ObimNeuralGenes/SraRunTable.txt | cut -f 1 -d "," | less)

for acc in ${accessions[@]}; do prefetch.2.10.8 -C yes -p $acc -O $OUTPUT_DIR && fastq-dump.2.10.8 --split-e --skip-technical -O $OUTPUT_DIR $OUTPUT_DIR/$acc.sra && rm $OUTPUT_DIR/$acc.sra; done

gzip $OUTPUT_DIR/*
