# Alignment

```bash
mkdir -p /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment

# python /private7/projects/Combinatorics/Code/align.py \
# --genome /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna \
# --in_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Data \
# --recursive \
# --out_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment \
# --processes 6 \
# --threads 6 \
# whole_transcriptome_illumina \
# --postfix ".fastq.bz2" \
# --decompress_cmd "bzip2 -dkc" \
# --known_sites_bed_file /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/HSPA1L_sites_in_transcript_zero_based.bed


# nohup \
# python /private7/projects/Combinatorics/Code/align.py \
# --genome /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna \
# --in_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Data \
# --recursive \
# --out_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment \
# --processes 3 \
# --threads 15 \
# whole_transcriptome_illumina \
# --postfix ".fastq.bz2" \
# --decompress_cmd "bzip2 -dkc" \
# --known_sites_bed_file /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/HSPA1L_sites_in_transcript_zero_based.bed \
# --separate_by_chrom_only_chroms_with_known_editing \
# > /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/align.test.9.10.24.out &


# samtools view -@ 10 -h \
# -o /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment/GTEX-1122O_WholeBlood_GTEX-1122O-0005-SM-5O99J.ByChrom/GTEX-1122O_WholeBlood_GTEX-1122O-0005-SM-5O99J.sorted.aligned.filtered.lcl|NC_000001.11_cds_XP_011540840.1_4.bam \
# /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment/GTEX-1122O_WholeBlood_GTEX-1122O-0005-SM-5O99J.sorted.aligned.filtered.bam \
# lcl|NC_000001.11_cds_XP_011540840.1_4


# samtools view -@ 10 -h \
# -o "/private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment/GTEX-1122O_WholeBlood_GTEX-1122O-0005-SM-5O99J.ByChrom/GTEX-1122O_WholeBlood_GTEX-1122O-0005-SM-5O99J.sorted.aligned.filtered.lcl|NC_000001.11_cds_XP_011540840.1_4.bam" \
# /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment/GTEX-1122O_WholeBlood_GTEX-1122O-0005-SM-5O99J.sorted.aligned.filtered.bam \
# "lcl|NC_000001.11_cds_XP_011540840.1_4"

nohup \
python /private7/projects/Combinatorics/Code/align.py \
--genome /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/GCF_000001405.40_GRCh38.p14_rna.fna \
--in_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Data \
--recursive \
--out_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment \
--processes 6 \
--threads 6 \
whole_transcriptome_illumina \
--postfix ".fastq.bz2" \
--decompress_cmd "bzip2 -dkc" \
--known_sites_bed_file /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/3genes_sites.bed \
--separate_by_chrom_only_chroms_with_known_editing \
> /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/align.15.10.24.out &
```
* alu 18
* 11:52



# Editing detection & distinct proteins finding

## Illumina

### Pileup

```bash
nohup python Code/pileup.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--data_table Code/Alignment/DataTable.Squid.ShortReads.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/Illumina \
--top_x_noisy_positions 3 \
--assurance_factor 1.5 \
--min_percent_of_max_coverage 0.1 \
--snp_noise_level 0.1 \
--include_flags 3 \
--exclude_flags 2304 \
--min_bq 30 \
--parity PE \
--processes 4 \
--threads 10 \
--keep_pileup_files \
>> D.pealeii/MpileupAndTranscripts/Illumina/pileup.out & 
```
 
### Distinct proteins

:warning: **Pay attention!** Running the distinct proteins & expression levels simulations for some of these Illumina transcripts will require a computer with **~2TB of memory**.

#### Finding isoforms

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina/*.unique_proteins.csv)

JULIA_PROJECT=.

nohup julia \
--threads 50 --proc 3 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/Illumina \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/Illumina/maximal_independent_set_5.Proteins.out &
```

#### Calculating expression levels

```bash
python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir D.pealeii/MpileupAndTranscripts/Illumina \
--proteins_prefix "reads.sorted.aligned.filtered."

DISTINCTFILES=$(cat D.pealeii/MpileupAndTranscripts/Illumina/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat D.pealeii/MpileupAndTranscripts/Illumina/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat D.pealeii/MpileupAndTranscripts/Illumina/ChromsNamesForExpressionLevels.txt)

nohup \
julia \
--project=. \
--threads 40 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir D.pealeii/MpileupAndTranscripts/Illumina \
> D.pealeii/MpileupAndTranscripts/Illumina/expressionlevels.out &
``` 


## O.vul polished IsoSeq data

total_mapped_reads 50

### Pileup

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/pileup.out & 
```


### Distinct proteins

#### Finding isoforms

```bash
mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/DistinctProteins.regular.log
```


## Human

### Pileup

```bash
mkdir -p /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/MpileupAndTranscripts

nohup python Code/pileup_with_subparsers.py \
--transcriptome  /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/GCF_000001405.40_GRCh38.p14_rna.fna \
--known_editing_sites /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/3genes_sites.bed \
--include_flags 3 \
--exclude_flags 2304 \
--min_bq 30 \
--parity PE \
--out_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/MpileupAndTranscripts \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Alignment/ByChrom \
--interfix_start ".sorted.aligned.filtered" \
--cds_regions /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/Annotations/HSPA1L_and3genes_coding_regions_from_transcript.bed \
--min_mapped_reads_per_position 0 \
> /private10/Projects/Maor/combinatorial_analysis/collab_with_kobi/pileup.29.10.24.out &
```
* alu 18
* 17:44