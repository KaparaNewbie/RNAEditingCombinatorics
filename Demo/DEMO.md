

# Running the demo

It should take about 5 minutes to run this demo.

```bash
# alignment

python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir Demo \
--out_dir Demo \
pacbio

# pileup (editing per position, trasncript and predicted protein)

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir Demo \
--processes 1 \
--threads 10 \
--gz_compression \
directed_sequencing_data \
--data_table Demo/DataTable.Squid.LongReads.Demo.csv \
> Demo/pileup.out &

# finding distinct proteins

INFILES=$(echo Demo/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 10 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir Demo \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded

# observing the distinct proteins' algorithm results

DISTINCT_FILE="Demo/PCLO-CNS-RESUB.DistinctUniqueProteins.<time-stamp>.csv" # replace the <time-stamp> with that of the actual file you just created
head -n1 <(cut -f 1-5 $DISTINCT_FILE) && tail -n6 <(cut -f 1-5 $DISTINCT_FILE)
```
The output of the above command should look like this:


| Fraction | FractionRepetition | Algorithm | AlgorithmRepetition |	NumUniqueSamples |
| -------- | ------- | -------- | ------- | -------- | 
| 1.0 | 3 | Descending | 1 |	3102 |
| 1.0 | 3 | Descending | 2 |	3102 |
| 1.0 | 4 | Ascending | 1 |	3093 |
| 1.0 | 4 | Ascending | 2 |	3094 |
| 1.0 | 4 | Descending | 1 |	3103 |
| 1.0 | 4 | Descending | 2 |	3102 |



<br>
<br>

# BTS: setting the demo

```bash
mkdir -p Demo/Data/CCS

samtools view \
--with-header \
--subsample 0.1 \
--subsample-seed 1892 \
--threads 20 \
--bam \
--output Demo/PCLO-CNS-RESUB.C0x1291.ccs.bam \
D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam
```