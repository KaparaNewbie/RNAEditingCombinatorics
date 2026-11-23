> **A-to-I editing generates unparalleled complexity in the neural transcriptome of cephalopods**
>
> _Kobi Shapira, Ruti Balter, David Gorelik, Joshua J C Rosenthal, Erez Y. Levanon & Eli Eisenberg_

**Revised-manuscript workflow.** This document describes the manuscript analyses, their processing commands and the producer–consumer chain for figures and tables. Commands target the current tracked source where equivalent behavior is established; qualified historical steps and unresolved reproduction requirements are identified explicitly.

For additional historical processing commands, alternative analyses, and extended development records, see the original [README.md](README.md). That document preserves provenance; it is not an independently validated reproduction path.

The workflow is:

**Original reads and annotations → preprocessing and alignment → editing detection, translation and isoform selection → scientific analysis notebooks → generated scientific tables → `final_manuscript_figures` → publication figures.**

The primary squid analyses use the original targeted PacBio data (PB1), the original pooled Illumina data, and the additional PacBio biological replicates (PB3). Octopus Iso-Seq and single-cell data support the transcriptome-wide comparisons. The earlier PacBio Data 2 dataset is outside the revised manuscript workflow.

# Installation and software requirements

Use separate environments for Python analysis, PacBio processing and Julia. Run
processing commands from the repository root unless they explicitly change directory.
Install Julia and Docker separately; Conda alone does not provide this complete stack.

```bash
git clone https://github.com/KaparaNewbie/RNAEditingCombinatorics
cd RNAEditingCombinatorics
conda env create -f combinatorics2_simplified.yml
conda env create -f pacbiocomb.yml
conda activate combinatorics2
MANUSCRIPT_REPO="$PWD"
alias COMB='cd "$MANUSCRIPT_REPO" && conda activate combinatorics2'
alias PBCOMB='cd "$MANUSCRIPT_REPO" && conda activate pacbiocomb'
# With Julia 1.7.2 installed and on PATH:
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

[combinatorics2_simplified.yml](combinatorics2_simplified.yml) is the current broad
specification; [combinatorics2.yml](combinatorics2.yml) records the earlier environment.
The current specification pins scikit-learn 1.9.0 because the accepted S7 code checks that version; other listed versions describe the validated runtime. Neither file is an immutable lockfile. `pandas>=2.2` is needed by the selected notebook groupby cells. [pacbiocomb.yml](pacbiocomb.yml) supplies PacBio
executables and need not contain Python. Historical `*.explicit.txt` filename lists
are not complete Conda URL lockfiles. Do not update an existing shared environment
for a reproduction run. Use isolated prefixes when testing installation.

| Software / environment | Role and version requirements |
|---|---|
| combinatorics2 | Current Python processing and scientific notebooks. Tested installed stack: Python 3.13.14, NumPy 2.4.6, pandas 3.0.3, SciPy 1.18.0, scikit-learn 1.9.0, statsmodels 0.14.6. These identify the validation runtime, not arbitrary new package pins. |
| pacbiocomb | CCS/HiFi, pbmm2, pbmarkdup, extracthifi, lima and Iso-Seq; preserve the stage-specific versions below. |
| Julia 1.7.2 | Selection, reassignment and S19 simulation; use the tracked [Project.toml](Project.toml) and [Manifest.toml](Manifest.toml), including the custom BioAlignments source. |
| Graphics | Matplotlib 3.11.0, Plotly 6.8.0, Kaleido 1.3.0, Pillow 11.3.0; Kaleido needs an existing Chrome executable. `rsvg-convert`, librsvg/Cairo and Poppler (`pdfinfo`, `pdftoppm`) are external system requirements. Matching Open Sans/Arial files are explicit external assets. |
| Docker / PyMuPDF | Final physical cropping and Figure 2 assembly use `combinatorics-figure2:2026-10-13`, Python 3.11.16 / PyMuPDF 1.26.7. See the external-image limitation below. PyMuPDF is not installed in the inspected combinatorics2 environment. |
| Jupytext 1.19.4 | Reconstruct `.ipynb` from Git-trackable `.py:percent` sources; conversion does not execute cells. |
| Sequencing/annotation tools | samtools 1.18 for the recorded PB1/complete-noise stages; 1.22.1 for PB3 March 2026; Salmon 1.7.0 and OrthoFinder 2.5.5 for neural annotations. BWA, FastQC, MultiQC, PRINSEQ-lite, Trimmomatic/adapters and PacBio commands are stage-specific prerequisites. |
| Numerical controls | S7 requires seven initialized OpenBLAS threads for historical distance reduction and one fitting thread. The validated environment records libopenblas 0.3.28 (`pthreads_h94d23a6_1`); the broad Conda specification does not lock this build. S19 requires one Julia thread and seed 1892. Preserve fixed solution identities, sampling seeds and numerical reduction order. |

Adapt lab paths to the same named native inputs, preserving their reference versions,
cohorts and selected solutions. The original CentOS 7 environment and some historical
binaries may need a compatible host. Raw Illumina selection/reassignment can require
about 2 TB RAM; the retained-processed-input figure route does not rerun that inference.
Verify executable versions for each stage, even when the command uses the `COMB` alias.

## Processing command compatibility

The alignment parsers and legacy pileup argument handoff have bounded compatibility checks. These checks do not establish an end-to-end scientific reproduction. Preserve the stage-specific inputs and external-tool versions below.

| Stage | Current-source compatibility and required action |
|---|---|
| Read preparation | PB1 now supplies `pacbio` with the established CCS defaults explicitly. Illumina and polished Iso-Seq subcommands match the parser. Current setup resolves FastQC, MultiQC, PRINSEQ-lite, Trimmomatic and its adapter file while constructing the parser, including for PacBio; provide a complete environment. Preserve the recorded sequencing-tool versions and adapter identity. |
| Alignment | The duplicate PacBio option is repaired: `--separate_by_chrom` and `--merge_by_chrom_bams` are independent flags. Parser construction and callable argument binding passed for `pacbio`, `illumina` and `whole_transcriptome_isoseq` in `combinatorics2`. PB1 supplies the subcommand and required known-sites argument below. No alignment was rerun. |
| Original Illumina pileup | The legacy `pileup.py` now passes `override_existing_pileup_files=True`, preserving historical overwrite behavior. `--skip_existing_pileup_files` reuses existing pileups; downstream stages still run. Mock tests verified argument order, flags and both modes without native processing. The unused notebook-execution import is commented out; direct import and parser construction pass in `combinatorics2` without `papermill`. Other notebook-execution utilities retain their own dependencies. |
| Directed and pooled pileup | Listed options match `pileup_with_subparsers.py`. Its default samtools path is `~/anaconda3/envs/combinatorics2/bin/samtools`, regardless of the active shell environment. PB1 explicitly selects the recorded 1.18 executable below; PB3 requires 1.22.1; July complete-noise commands explicitly select 1.18. For other stages, establish the matching executable/version before a fresh run. |
| Isoform selection and expression | Julia option names match the tracked parsers. All-cell 1,000 and neuronal/non-neuronal 50-read expression recipes are no longer prerequisites: their notebook cohorts now use the same ordered distinct-solution inputs. The revised all-cell 50-read entropy-considered producer remains required. The older all-cell 50-read expression branch is historical and is not required by the revised manuscript. Preserve `firstcolpos` (15 squid, 16 octopus), fractions, algorithms and selected files. MIS defaults to seed 1892; expression defaults to fraction 1.0 and Ascending/Descending. |
| Annotation and auxiliary programs | ORF/BED conversion, known-site parsing and O. bim neural annotation match current interfaces. File-list preparation requires `--reads_dir`; the revised single-cell recipe supplies it. Cloud transfers and simulation/chimera scripts contain site-specific settings. |

The manuscript adapters require Python 3.11+ (including `hashlib.file_digest`); a historical environment name alone does not establish compatibility. The modern specification is the starting point for current Python source, while external tools still require the stage-specific versions.


For a small pipeline illustration, see [Demo/DEMO.md](Demo/DEMO.md). It is not a manuscript reproduction.

# Reference data and annotations

Reference versions are part of the input definition. The original squid transcriptome and the extended Jan2025 transcriptome are distinct local files. Preserve the version named in each command; a new download from a mutable URL is not evidence of identity with the original input.

| Reference or annotation | Local location and role |
|---|---|
| Original squid transcriptome | `D.pealeii/Annotations/orfs_squ.fa`; original PB1 and Illumina processing |
| Extended squid transcriptome | `D.pealeii/Annotations/Jan2025/orfs_squ.fa`; PB3 and July 2026 complete-noise processing |
| Squid editing sites | Original `D.pea.EditingSites.bed` and extended `Jan2025/D.pea.EditingSites.bed` / `.csv`; use the version paired with the reference |
| Octopus transcriptome and editing sites | `O.vulgaris/Annotations/orfs_oct.fa`, `O.vul.EditingSites.bed` / `.csv` |
| O. bimaculoides transcriptome | `O.bim/Annotations/orfs_bim.fa`; expression and orthology-based neural annotation |
| ORF intervals | Species-specific `orfs_*.bed`, generated by [Code/orfs_fasta_to_bed.py](Code/orfs_fasta_to_bed.py) from FASTA header coordinates |
| Single-cell transcript annotation | `O.vulgaris/Annotations/oct.bed` → `oct.gff` → `oct.gtf`, using the commands below; the exact reference/index handoff is also required from the upstream repository |
| O. bimaculoides GTF | `O.bim/Annotations/orfs_bim.gtf`; producer [prepare_salmon_gtf_input_from_transcriptome.py](Code/NeuralTranscripts/prepare_salmon_gtf_input_from_transcriptome.py). The current `neural_transcripts.py` route indexes the transcript FASTA directly. |

The original squid FASTA has SHA-256 `49f094d8f83aabb5961760353aca9d2d349c185a19a14451e78dd0c45f1c6e09`; the extended Jan2025 file has SHA-256 `b19dd2d1eb6aa9e8601fbc93b00145fe787b208bc7f5c7d119eaad08e97f1221`. These identify the inspected local copies. The acquisition/extension procedure for the Jan2025 FASTA and editing annotations remains **TO VERIFY**; the original download commands below do not generate that extended set.

### Original squid reference and editing annotations

```bash
mkdir -p D.pealeii/Annotations
cd D.pealeii/Annotations

wget https://www.tau.ac.il/~elieis/squid/orfs_squ.fa
conda activate combinatorics
samtools faidx orfs_squ.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx
wget https://www.tau.ac.il/~elieis/squid/ST4.txt.gz
wget https://www.tau.ac.il/~elieis/squid/ST5.txt.gz
gunzip ST4.txt.gz
gunzip ST5.txt.gz

COMB  # return to the main project' dir
```

Export the original editing-site sheet as CSV (also used by the simulation producer) and convert it into a six-column BED file. The start coordinate is zero-based; the end is the original one-based editing location.

```python
import pandas as pd
df = pd.read_excel("D.pealeii/Annotations/editing_sites.xlsx", sheet_name="D.pea Editing sites")
df.to_csv("D.pealeii/Annotations/D.pea.EditingSites.csv", index=False)
bed_cols = ["Trinity id", "Editing location (base1)", "SwissProt name"]
df2 = df.filter(bed_cols).rename(
    columns={
        "Trinity id": "#TrinityID",
        "Editing location (base1)": "End",
        "SwissProt name": "SwissProtName",
    }
)
df2.insert(1, "Start", df2["End"] - 1)
df2["Score"] = "."
df2["Strand"] = "+"
df2.to_csv("D.pealeii/Annotations/D.pea.EditingSites.bed", index=False, sep="\t")
```

Extract ORF coordinates from the original transcriptome:

```bash
python Code/orfs_fasta_to_bed.py \
--in_fasta D.pealeii/Annotations/orfs_squ.fa \
--out_bed D.pealeii/Annotations/orfs_squ.bed
```

### Octopus reference and editing annotations

```bash
mkdir O.vulgaris/Annotations
cd O.vulgaris/Annotations
wget https://www.tau.ac.il/~elieis/squid/orfs_oct.fa
samtools faidx orfs_oct.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx

COMB  # return to the main dir of the project

python \
Code/parse_known_sites.py \
--editing_sites_excel_file O.vulgaris/Annotations/editing_sites.xlsx \
--sheet_name "O.vul Editing sites" \
--out_csv_file O.vulgaris/Annotations/O.vul.EditingSites.csv \
--out_bed_file O.vulgaris/Annotations/O.vul.EditingSites.bed \
--ignore_score

python Code/orfs_fasta_to_bed.py \
--in_fasta O.vulgaris/Annotations/orfs_oct.fa \
--out_bed O.vulgaris/Annotations/orfs_oct.bed
```

### Octopus BED, GFF and GTF for the single-cell handoff

The original transcript BED construction and conversions are retained below. The conversion commands name the AGAT container tag `0.8.0--pl5262hdfd78af_0`; they do not establish an immutable image digest.

```python
from pathlib import Path

import pandas as pd
from Bio import SeqIO

in_fasta = Path("O.vulgaris/Annotations/orfs_oct.fa")
out_bed = Path("O.vulgaris/Annotations/oct.bed") # a bed file of the complete transcriptome

records = SeqIO.parse(in_fasta, "fasta")

transcripts = []
transcripts_starts = []
transcripts_ends = []
names = []
strands = []

orfs_starts = []
orfs_ends = []

for record in records:
    description = record.description.split("\t")

    transcript_start = 0
    transcript_end = len(record.seq)
    transcript_name = description[-1].split()[0].split("|")[-1]
    strand_index = description.index("Strand") + 1
    strand = description[strand_index]

    orf_start_index = description.index("OrfStart") + 1
    orf_end_index = description.index("OrfEnd") + 1
    orf_start = int(description[orf_start_index]) - 1
    orf_end = int(description[orf_end_index])

    transcripts.append(record.id)
    transcripts_starts.append(transcript_start)
    transcripts_ends.append(transcript_end)
    names.append(transcript_name)
    strands.append(strand)

    orfs_starts.append(orf_start)
    orfs_ends.append(orf_end)

transcripts_strands_df = pd.DataFrame(
    {
        "Chrom": transcripts,
        "Start": transcripts_starts,
        "End": transcripts_ends,
        "Name": names,
        "Strand": strands,
        "ThickStart": orfs_starts,
        "ThickEnd": orfs_ends,
        "itemRgb": ".",
        "blockCount": 1,
        "blockStarts": 0,
    }
)
transcripts_strands_df = transcripts_strands_df.sort_values(["Chrom", "Start"])
transcripts_strands_df.insert(
    transcripts_strands_df.columns.get_loc("Strand"), "Score", "."
)
transcripts_strands_df.insert(
    transcripts_strands_df.columns.get_loc("blockCount") + 1, "BlockSizes", transcripts_strands_df["ThickEnd"] - transcripts_strands_df["ThickStart"]
)

transcripts_strands_df.to_csv(out_bed, sep="\t", index=False, header=False)
```

```bash
cd O.vulgaris/Annotations

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_bed2gff.pl \
--bed /dir/oct.bed \
--primary_tag mRNA \
-o /dir/oct.gff

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_sp_gff2gtf.pl \
--gff /dir/oct.gff \
-o /dir/oct.gtf
```

Return to the repository root after the conversion block before continuing with the dataset commands.

### O. bimaculoides reference

```bash
mkdir O.bim/Annotations
cd O.bim/Annotations
wget https://www.tau.ac.il/~elieis/squid/orfs_bim.fa
samtools faidx orfs_bim.fa

COMB # back to the main directory of the project
```

# Analyses

## Original squid PacBio (PB1)

| Public accession | Manuscript target | Original local input |
|---|---|---|
| SRR28713506 | GRIA2 (called `GRIA` in older tables/commands) | `D.pealeii/Data/Raw/GRIA-CNS-RESUB.C0x1291.subreads.bam` |
| SRR28713505 | PCLO | `D.pealeii/Data/Raw/PCLO-CNS-RESUB.C0x1291.subreads.bam` |

These are **original PacBio subreads BAMs**. The accession-to-target assignments are author supplied; the filenames are verified against local metadata and [DataTable.Squid.LongReads.csv](Code/Alignment/DataTable.Squid.LongReads.csv). The native target regions are `comp141693_c0_seq1` and `comp141882_c0_seq14`.

### Prepare consensus reads

Place the original demultiplexed subreads in the input directory and prepare CCS reads. Native records support **CCS 6.2.0** for the original PB1 consensus generation. These inputs require subread-to-consensus processing; PB3 below starts from delivered HiFi reads.

```bash
mkdir -p D.pealeii/Data/Raw
```

```bash
mkdir D.pealeii/Data/CCS

nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Raw \
--out_dir D.pealeii/Data/CCS \
--processes 6 \
--threads 10 \
pacbio \
--postfix .subreads.bam \
--min_passes 3 \
--min_rq 0.99 \
>> D.pealeii/Data/CCS/prepare_data.out &
```

This invocation uses the current `pacbio` interface. Source comparison with the earliest tracked implementation (`71c2b76`) found identical `basic_ccs`, `extract_hifi` and `pacbio_main` functions and unchanged defaults: six processes, ten threads, `.subreads.bam`, three passes and CCS predicted accuracy 0.99. The command makes those defaults explicit. PacBio tools run in `pacbiocomb`; preserve CCS 6.2.0. The downstream alignment consumes **BasicCCS**, even though preparation also creates HiFi and deduplicated branches. This source comparison does not establish byte-for-byte identity of newly generated reads.

### Align

Original PB1 alignments used **pbmm2 1.7.0**. The resulting GRIA/PCLO `.aligned.sorted.bam` files are selected by the native long-read data table.

```bash
mkdir -p D.pealeii/Alignment/BestN1

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir D.pealeii/Data/CCS/BasicCCS \
--out_dir D.pealeii/Alignment/BestN1 \
pacbio \
--known_sites_bed_file D.pealeii/Annotations/D.pea.EditingSites.bed \
> D.pealeii/Alignment/BestN1/align.out &
```

The current command preserves the original CCS preset, best-n = 1, two processes and 25 threads. Separation and merging remain disabled. The parser requires a known-sites BED, but the PB1 execution path returns before using it when separation is disabled. Only the CLI spelling/required arguments were repaired; the alignment implementation is unchanged.

### Detect editing and translate reads

The original February 2024 PB1 pileup used **samtools 1.18**. Retain the original squid reference, SE handling, minimum read quality 0.998 and base quality 30.

```bash
COMB
SAMTOOLS_118=/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools
"$SAMTOOLS_118" --version # must identify the recorded samtools 1.18

mkdir -p D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--samtools_path "$SAMTOOLS_118" \
--out_dir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--processes 2 \
--threads 10 \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.LongReads.csv \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/pileup.out &
```

This stage produces native position, read, unique-read, protein and unique-protein tables in `RQ998.TopNoisyPositions3.BQ30`. Their read memberships and original multiplicities remain essential for downstream analyses.

### Select distinct isoforms and assign expression

Julia generates repeated heuristic independent-set solutions. The largest saved solution is a selected heuristic result, not a proof of the global maximum. Assigned expression includes eligible support reassigned from compatible reads; it is distinct from the count of original reads directly supporting an isoform.

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.regular.log
```

The entropy-considered assigned-expression run uses the specified native solution files and reads:

```bash
WD="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30"


nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $WD/GRIA-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:29:20.csv $WD/PCLO-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:46:24.csv \
--allprotsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames GRIA PCLO \
--firstcolpos 15 \
--fractions 1.0 \
--outdir $WD \
--considerentropy \
--postfix_to_add .EntropyConsidered \
> $WD/expressionlevels.regular.7.8.2025.out &
```

### Ten-percent comparison for Figure S10

Keep the 0.1-fraction selections and assigned-expression outputs separate from the full-data selection.

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)

# julia \
# --project=. \
# --threads 40 --proc 6 \
# Code/Simulations/maximal_independent_set_5.jl \
# --infiles $INFILES \
# --postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
# --postfix_to_add .Fraction0_1 \
# --idcol Protein \
# --firstcolpos 15 \
# --datatype Proteins \
# --outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
# --maxfrac 0.1 \
# --fracstep 0.1 \
# --fracrepetitions 4 \
# --algrepetitions 2 \
# --algs Ascending Descending \
# --run_solve_threaded \
# 2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.regular.Fraction0_1.log

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .Fraction0_1 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--maxfrac 0.1 \
--fracstep 0.1 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.regular.Fraction0_1.log
```

```bash
# WD="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30"

# nohup \
# julia \
# --project=. \
# --threads 40 \
# Code/Simulations/expressionlevels.jl \
# --distinctfiles $WD/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.*.csv $WD/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.*.csv \
# --allprotsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
# --allreadsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
# --samplenames GRIA PCLO \
# --postfix_to_add .Fraction0_1 \
# --firstcolpos 15 \
# --fractions 0.1 \
# --outdir $WD \
# > $WD/expressionlevels.regular.fraction01.28.7.2025.out &


WD="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30"

nohup \
julia \
--project=. \
--threads 80 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $WD/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.06.02.2024-10:55:32.csv $WD/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.06.02.2024-11:05:02.csv \
--allprotsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames GRIA PCLO \
--innerthreadedassignment \
--considerentropy \
--postfix_to_add .Fraction0_1.EntropyConsidered \
--firstcolpos 15 \
--fractions 0.1 \
--outdir $WD \
> $WD/expressionlevels.regular.fraction01.entropy_considered.26.2.2026.out &
```

### Amino-acid dissimilarity comparisons for Figure 3

The retained comparisons use Miyata amino-acid groups and the Grantham 100 threshold. Preserve the specified grouping, score cutoff and validator.

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .AAgroupsMiyata1979 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--aagroups AA_groups_Miyata1979 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.AAgroupsMiyata1979.log
```

```bash

DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.06.02.2024-13:43:37.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.06.02.2024-14:14:55.csv
"""
ALLROTSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz
"""
ALLREADSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz
"""
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 80 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--allreadsfiles $ALLREADSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--aagroups AA_groups_Miyata1979 \
--innerthreadedassignment \
--considerentropy \
--postfix_to_add .AAgroupsMiyata1979.EntropyConsidered \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.AAgroupsMiyata1979.EntropyConsidered.26.2.2026.out &
```

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .GRANTHAM1974-100 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.GRANTHAM1974-100.log
```

```bash

DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.06.02.2024-14:27:55.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.06.02.2024-15:51:50.csv
"""
ALLROTSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz
"""
ALLREADSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz
"""
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 50 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--allreadsfiles $ALLREADSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--innerthreadedassignment \
--considerentropy \
--postfix_to_add .GRANTHAM1974-100.EntropyConsidered \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.GRANTHAM1974-100.EntropyConsidered.26.2.2026.out &
```

The Grantham expression command retains the effective output suffix `.GRANTHAM1974-100.EntropyConsidered`. The earlier duplicate option was removed: ArgParse 1.2.0 (the tracked manifest version) stores the last value for this scalar argument. Matrix, cutoff, predicate, input selections and entropy settings are unchanged.

The scientific producer is [squids_pacbio_5.BQ30.py](Code/Notebooks/squids_pacbio_5.BQ30.py). Its selected-expression, site-correlation, recoding, entropy and original-read-support sections feed the tables described below.


Domain results use the existing `revision-20261109-domains-*` cells and the accepted external ScanProsite response; see the notebook execution table.

## Squid Illumina

**SRR28713507 — PooledShortReads** is paired-end FASTQ data for the original 19-transcript analysis. The verified local raw files are:

- `D.pealeii/Data/RutisReads/Raw/reads_1.fastq.gz`
- `D.pealeii/Data/RutisReads/Raw/reads_2.fastq.gz`

The original commands use `D.pealeii/Data/Illumina/Raw` and `D.pealeii/Data/Illumina/TrimmedWoDup`. That raw directory is absent in the inspected checkout. The directory handoff from `RutisReads` to the command layout is **TO VERIFY**; the original commands are preserved below. The repeated mate-1 filename in the original README prose is resolved here using the actual mate-1 and mate-2 filenames.

[DataTable.Squid.ShortReads.csv](Code/Alignment/DataTable.Squid.ShortReads.csv) maps target regions to per-transcript BAMs. It contains 20 rows, including `VGA_BPS13`; the final 19-transcript population is the explicit `conditions`/`chroms` selection in [squids_illumina_5.py](Code/Notebooks/squids_illumina_5.py), which excludes that row.

### Prepare and align

The original preparation used mean-quality filtering and trimming to 265 bases. Trimmomatic **0.39** and PRINSEQ-lite **0.20.4** are supported historical environment candidates; their exact execution versions were not independently printed by the inspected preparation record. Original alignment/filtering is directly recorded as **bwa-mem2 2.2.1 and samtools 1.14**.

```bash
mkdir -p D.pealeii/Data/Illumina/Raw
```

```bash
nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Illumina/Raw/ \
--out_dir D.pealeii/Data/Illumina \
illumina \
--min_qual_mean 30 \
--trim_to_len 265 \
> D.pealeii/Data/Illumina/prepare_data.out &
```

The current alignment parser accepts this command; preserve the recorded external-tool versions and input population.

```bash
mkdir -p D.pealeii/Alignment/Illumina

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir D.pealeii/Data/Illumina/TrimmedWoDup \
--out_dir D.pealeii/Alignment/Illumina \
--threads 40 \
illumina \
--require_flags 3 \
--exclude_flags 2304 \
--separate_by_chrom \
> D.pealeii/Alignment/Illumina/align.out &
```

### Detect editing and select distinct isoforms

The original Illumina pileup command is retained with its original flags and thresholds. The alignment version above does not establish the executable version for every subsequent pileup.

The internal argument mismatch is repaired. The command below keeps historical overwrite behavior; add `--skip_existing_pileup_files` only when intentionally reusing pileups. The script loads directly in `combinatorics2` after commenting its unused notebook-execution import; no environment change is required for that repair.

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

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina/*.unique_proteins.csv)

nohup julia \
--project=. \
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

### Revised assigned-expression results

The retained scientific analysis uses the fixed entropy-considered expression outputs in `D.pealeii/MpileupAndTranscripts/Illumina/FinalFixedExpressionFromCloud/<gene>/`. The repository's native expression producer is [expressionlevels.jl](Code/Simulations/expressionlevels.jl); the final notebook's input definitions select the retained per-gene solutions and assignments.

The original upload and download entry points are:

```bash
python Code/upload_illumina_samples_to_cloud.py
```

```bash
python Code/download_illumina_fixed_results_from_cloud.py
```

These commands depend on the author's Google Cloud bucket and authenticated access. The uploader explicitly skips RUSC2 and RIMS2, which were handled separately. The downloader retrieves 19 result directories; downloading existing results does not regenerate them. The exact cloud execution configuration and complete per-gene handoff are not documented in the tracked route. The earlier non-entropy-considered expression command is not a substitute for these retained assignments.

### Matched-depth comparison for Figure S6

This branch samples **80,000 read pairs**, with `--consider_parity_when_sampling` and seed 1892. It supplies the fixed-sampling comparison in `various_plots.BQ30`; it does not replace the full-data Illumina population.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/Illumina80KPairs

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--include_flags 3 \
--exclude_flags 2304 \
--top_x_noisy_positions 3 \
--assurance_factor 1.5 \
--min_percent_of_max_coverage 0.1 \
--snp_noise_level 0.1 \
--min_bq 30 \
--parity PE \
--out_dir D.pealeii/MpileupAndTranscripts/Illumina80KPairs \
--processes 10 \
--threads 10 \
--keep_pileup_files \
--gz_compression \
--sample_reads \
--num_sampled_reads 80000 \
--consider_parity_when_sampling \
--seed 1892 \
directed_sequencing_data \
--data_table D.pealeii/Alignment/Illumina/reads.ByChrom/data_table.csv \
> D.pealeii/MpileupAndTranscripts/Illumina80KPairs/pileup.4.1.26.out &
```

```bash
tmux new -s Illumina80K # alu 18

INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina80KPairs/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 60 --proc 10 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/Illumina80KPairs \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/Illumina80KPairs/DistinctProteins.regular.4.1.26.log
```

## Additional squid biological replicates (PB3)

| Public accession | Biological replicate | Original local input |
|---|---|---|
| SRR41297225 | JR1 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR1.r841072252411_C01.hifireads.bam` |
| SRR41297224 | JR2 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR2.r841072252411_C01.hifireads.bam` |
| SRR41297223 | JR3 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR3.r841072252411_C01.hifireads.bam` |

These are delivered **HiFi/CCS BAMs with UMIs**, not raw subreads. The accession-to-replicate assignments are author supplied and the exact filenames are verified locally. Provider records identify CCS **8.2.0**, lima **2.11.0**, pbtrim **1.1.0** and jasmine **2.3.0**. The local pipeline starts from those delivered files.

PB3 targets GRIA2 (`comp141693_c0_seq1`), ADAR1 (`comp134400_c0_seq1_extended`) and IQEC1 (`comp141565_c6_seq3`). Keep biological replicate identity through UMI processing, per-site comparisons and assigned-expression summaries.

### Align the delivered HiFi reads

March 2026 processing used local **pbmm2 26.1.99** and **samtools 1.22.1** with the extended Jan2025 reference. The native March log records `merge_samples_bams_per_chrom=True`; the three target BAM headers record `samtools merge` of the JR1/JR2/JR3 per-transcript alignments. Thus merging was already used before the April Git commit. The exact uncommitted March source snapshot has not been recovered. The current option spelling is `--merge_by_chrom_bams`.

The current alignment parser accepts this command; preserve the recorded external-tool versions and input population.

```bash
mkdir -p D.pealeii/Alignment/AdditionalUMILongReads

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--in_dir D.pealeii/Data/AdditionalRawWithUMIs/CCS \
--out_dir D.pealeii/Alignment/AdditionalUMILongReads \
--processes 3 \
--threads 15 \
--ansi_color \
pacbio \
--postfix ".r841072252411_C01.hifireads.bam" \
--separate_by_chrom \
--merge_by_chrom_bams \
--known_sites_bed_file D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
> D.pealeii/Alignment/AdditionalUMILongReads/align.19.3.2026.out &
```

### All-read pileup, isoform selection and expression

The March 2026 pileup also used **samtools 1.22.1**. The canonical All-reads input is [Code/Alignment/DataTable.Squid.PB3.AllReads.csv](Code/Alignment/DataTable.Squid.PB3.AllReads.csv). It selects the original merged BAMs under `D.pealeii/Alignment/AdditionalUMILongReads/ByChrom/`, before barcode recognition or deduplication, with the verified `GRIA`, `ADAR1` and `IQEC1` labels, Jan2025 coordinates/strands and no problematic-region mask. Native headers, execution logs, output identities and the current/historical routing agree on those inputs.

The historical Sample labels are unknown. The replacement uses explicit `PB3_AllReads_*` labels: in this directed route, with `filter_col=None`, Sample is read but its values are never used for BAM selection, output names, read/protein labels, replicate membership or calculations. Changing those labels leaves every downstream argument unchanged. This is a scientifically equivalent input table, not a claim of historical CSV byte identity. The existing `DataTable.Squid.MergedAdditionalUMILongReads.csv` is preserved for its later recognizable-barcode paths; it must not replace the All-reads table. Recognizable-barcode and deduplicated analyses remain separate populations.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads


nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads \
--processes 3 \
--threads 15 \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.PB3.AllReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/pileup.19.3.2026.out &
```

```bash
tmux new -s COMB18

COMB

INFILES=$(echo D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 80 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove ".merged.MinRQ998.unique_proteins.csv.gz" \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/DistinctProteins.Regular.19.3.2026.log
```

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads

nohup \
julia \
--project=. \
--threads 80 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/comp141693_c0_seq1.DistinctUniqueProteins.19.03.2026-13:06:03.csv $UMI_DIR/comp134400_c0_seq1_extended.DistinctUniqueProteins.19.03.2026-13:03:20.csv $UMI_DIR/comp141565_c6_seq3.DistinctUniqueProteins.19.03.2026-13:22:23.csv \
--allprotsfiles $UMI_DIR/comp141693_c0_seq1.merged.MinRQ998.unique_proteins.csv.gz $UMI_DIR/comp134400_c0_seq1_extended.merged.MinRQ998.unique_proteins.csv.gz $UMI_DIR/comp141565_c6_seq3.merged.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/comp141693_c0_seq1.merged.MinRQ998.reads.csv.gz $UMI_DIR/comp134400_c0_seq1_extended.merged.MinRQ998.reads.csv.gz $UMI_DIR/comp141565_c6_seq3.merged.MinRQ998.reads.csv.gz \
--samplenames GRIA2 ADAR1 IQEC1 \
--firstcolpos 15 \
--onlymaxdistinct \
--innerthreadedassignment \
--considerentropy \
--outdir $UMI_DIR \
--postfix_to_add .EntropyConsidered \
> $UMI_DIR/expressionlevels.EntropyConsidered.19.3.2026.out &
```

### UMI processing and biological-replicate analysis

[additional_umi_long_reads_2b.py](Code/Notebooks/additional_umi_long_reads_2b.py) contains the native UMI/read preparation: barcode recognition, replicate-aware read mapping, “Save reads with recognized barcodes, before de-duplication”, and “Save reads after de-duplication”. It depends on the delivered reads, aligned per-transcript BAMs and the original all-read tables. Its current settings select the deduplicated output branch; reproducing the earlier recognizable-barcode handoff also requires the matching branch settings.

The main scientific notebook is [squids_pacbio_5.additional_UMIs.MergedSamples.py](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.py). It produces per-site replicate comparisons, selected-isoform recurrence and expression summaries for Figure 4b, Figure S4 and supplementary tables. Original direct read counts and reassigned expression mass have different meanings; recurrence/relative abundance is calculated within each gene and biological replicate.

### Reads with recognizable barcodes — Figure S10

Use the barcode-filtered BAMs and data table produced by the UMI notebook before the following stages.

```bash
OUT_DIR=D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir $OUT_DIR \
--processes 3 \
--threads 20 \
--gz_compression \
directed_sequencing_data \
--data_table $OUT_DIR/DataTable.Squid.MergedAdditionalUMILongReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> $OUT_DIR/pileup.16.4.2026.out &
```

```bash
tmux new -s COMB17

COMB

INFILES=$(echo D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 80 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove ".merged.MinRQ998.unique_proteins.csv.gz" \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes/DistinctProteins.Regular.16.4.2026.log
```

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes


nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/comp141693_c0_seq1.DistinctUniqueProteins.*.csv \
$UMI_DIR/comp134400_c0_seq1_extended.DistinctUniqueProteins.*.csv \
$UMI_DIR/comp141565_c6_seq3.DistinctUniqueProteins.*.csv \
--allprotsfiles $UMI_DIR/comp141693_c0_seq1.merged.MinRQ998.unique_proteins.csv.gz \
$UMI_DIR/comp134400_c0_seq1_extended.merged.MinRQ998.unique_proteins.csv.gz \
$UMI_DIR/comp141565_c6_seq3.merged.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/comp141693_c0_seq1.merged.MinRQ998.reads.csv.gz \
$UMI_DIR/comp134400_c0_seq1_extended.merged.MinRQ998.reads.csv.gz \
$UMI_DIR/comp141565_c6_seq3.merged.MinRQ998.reads.csv.gz \
--samplenames GRIA2 ADAR1 IQEC1 \
--firstcolpos 15 \
--onlymaxdistinct \
--innerthreadedassignment \
--considerentropy \
--outdir $UMI_DIR \
--postfix_to_add .EntropyConsidered \
> $UMI_DIR/expressionlevels.5.7.2026.out &
```

The corresponding scientific notebook retains the author's filename spelling: [squids_pacbio_5.additional_UMIs.MergedSamples.ReadsWithRecognizableBardoces.py](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.ReadsWithRecognizableBardoces.py).

### UMI-deduplicated reads — Figure S10

This branch uses the deduplicated BAM/data-table handoff and remains separate from both all reads and recognizable-barcode reads.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads

# cp \
# D.pealeii/Alignment/AdditionalUMILongReads/DataTable.Squid.MergedAdditionalUMILongReads.csv \
# D.pealeii/Alignment/AdditionalUMILongReads.DeduppedReads/DataTable.csv

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads \
--processes 3 \
--threads 15 \
--gz_compression \
directed_sequencing_data \
--data_table D.pealeii/Alignment/AdditionalUMILongReads.DeduppedReads/DataTable.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads/pileup.2.7.2026.out &
```

```bash
tmux new -s COMB17

COMB

INFILES=$(echo D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 80 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove ".merged.MinRQ998.unique_proteins.csv.gz" \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads/DistinctProteins.Regular.6.7.2026.log
```

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads

nohup \
julia \
--project=. \
--threads 40 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/comp141693_c0_seq1.DistinctUniqueProteins.*.csv \
$UMI_DIR/comp134400_c0_seq1_extended.DistinctUniqueProteins.*.csv \
$UMI_DIR/comp141565_c6_seq3.DistinctUniqueProteins.*.csv \
--allprotsfiles $UMI_DIR/comp141693_c0_seq1.merged.MinRQ998.unique_proteins.csv.gz \
$UMI_DIR/comp134400_c0_seq1_extended.merged.MinRQ998.unique_proteins.csv.gz \
$UMI_DIR/comp141565_c6_seq3.merged.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/comp141693_c0_seq1.merged.MinRQ998.reads.csv.gz \
$UMI_DIR/comp134400_c0_seq1_extended.merged.MinRQ998.reads.csv.gz \
$UMI_DIR/comp141565_c6_seq3.merged.MinRQ998.reads.csv.gz \
--samplenames GRIA2 ADAR1 IQEC1 \
--firstcolpos 15 \
--onlymaxdistinct \
--innerthreadedassignment \
--considerentropy \
--outdir $UMI_DIR \
--postfix_to_add .EntropyConsidered \
> $UMI_DIR/expressionlevels.6.7.2026.out &
```

Continue with [squids_pacbio_5.additional_UMIs.MergedSamples.DeduppedReads.py](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.DeduppedReads.py). The comparison is assembled in [various_plots.BQ30.py](Code/Notebooks/various_plots.BQ30.py) and rendered as Figure S10.

## Octopus Iso-Seq

**BioProject PRJNA791920** supplies the octopus long-read dataset. The polished CCS files were provided by the authors of *MicroRNAs are deeply linked to the emergence of the complex octopus brain*. The local input directory is `O.vulgaris/Data/PRJNA791920/IsoSeqPolished/CCS`; the supplied polished-read handoff must be preserved rather than assumed identical to an arbitrary reprocessing of public raw data.

The native `samples.csv` maps the seven runs as follows:

| Run | Tissue label in the native sample table |
|---|---|
| SRR17321895 | Axial nerve cord |
| SRR17321896 | Non-neuronal tissues mix |
| SRR17321897 | Frontal & Vertical lobe |
| SRR17321898 | Pedunculate & Olfactory lobe |
| SRR17321899 | Stellate g. & Visceral g. |
| SRR17321900 | Sucker |
| SRR17321901 | Retina & Optic lobe |

### Prepare and align

Sampled native records identify supplied **CCS 4.0.0**, local **lima 2.6.0**, **Iso-Seq refine 3.8.2** and **pbmm2 1.7.0**. Preparation uses the tracked `Code/Data/ClontechSMARTer-NEBcDNA.primers.fasta` primer set.

```bash
mkdir -p O.vulgaris/Data/PRJNA791920/IsoSeqPolished/CCS

nohup \
python Code/prepare_data.py \
--in_dir O.vulgaris/Data/PRJNA791920/IsoSeqPolished/CCS \
--out_dir O.vulgaris/Data/PRJNA791920/IsoSeqPolished \
--processes 7 \
--threads 6 \
pacbio_polished_ccs_isoseq \
> O.vulgaris/Data/PRJNA791920/IsoSeqPolished/prepare_data.out &
```

The current alignment parser accepts this command; preserve the recorded external-tool versions and input population.

```bash
mkdir -p O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered

nohup \
python Code/align.py \
--genome O.vulgaris/Annotations/orfs_oct.fa \
--in_dir O.vulgaris/Data/PRJNA791920/IsoSeqPolished/Refined \
--out_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered \
--processes 7 \
--threads 6 \
whole_transcriptome_isoseq \
--known_sites_bed_file O.vulgaris/Annotations/O.vul.EditingSites.bed \
--postfix .flnc.bam \
> O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/align.out &
```

### Pooled samples, minimum 50 mapped reads

The retained pooled analysis uses the sample table, per-transcript/per-sample alignment summary, explicit SNP/noise criteria and final editing scheme below. Preserve these criteria and the distinction between total mapped-read eligibility and per-position coverage.

```bash
OUTDIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples

mkdir -p $OUTDIR

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--snp_noise_level 0.05 \
--out_dir $OUTDIR \
--processes 40 \
--threads 5 \
--gz_compression \
--keep_pileup_files \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--total_mapped_reads 50 \
--pooled_transcript_noise_threshold 1 \
--max_snps_per_gene_to_allow_editing_detection 3 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--keep_bam_files \
> $OUTDIR/pileup.9.2.2026.out 2>&1 &
```

```bash
INDIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples
OUTDIR=$INDIR/DistinctProteins
mkdir -p $OUTDIR
INFILES=$(echo $INDIR/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir $OUTDIR \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee $INDIR/DistinctProteins.regular.9.2.2026.log
```

```bash
DIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples

python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir $DIR/ProteinsFiles \
--reads_dir $DIR/ReadsFiles \
--distinct_proteins_dir $DIR/DistinctProteins \
--proteins_postfix ".gz" \
--reads_postfix ".gz" \
--out_dir $DIR

OUTDIR=$DIR/ExpressionLevels

mkdir $OUTDIR

nohup \
julia \
--project=. \
--threads 40 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $DIR/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $DIR/UniqueProteinsForExpressionLevels.txt \
--allreadsfilesfofn $DIR/ReadsForExpressionLevels.txt \
--samplenamesfile $DIR/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--onlymaxdistinct \
--considerentropy \
--postfix_to_add .EntropyConsidered \
--innerthreadedassignment \
--outdir $OUTDIR \
> $DIR/expressionlevels.10.2.2026.out &
```

### Pooled samples, minimum 1,000 mapped reads

This branch supports the coverage comparison. It is a separate eligible transcript population.

```bash
OUTDIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.PooledSamples

mkdir -p $OUTDIR

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--snp_noise_level 0.05 \
--out_dir $OUTDIR \
--processes 40 \
--threads 5 \
--gz_compression \
--keep_pileup_files \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--total_mapped_reads 1000 \
--pooled_transcript_noise_threshold 1 \
--max_snps_per_gene_to_allow_editing_detection 3 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--keep_bam_files \
> $OUTDIR/pileup.11.2.2026.out 2>&1 &
```

```bash
INDIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.PooledSamples
OUTDIR=$INDIR/DistinctProteins
mkdir -p $OUTDIR
INFILES=$(echo $INDIR/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir $OUTDIR \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee $INDIR/DistinctProteins.regular.11.2.2026.log
```

The scientific producer is [o.vul_tmr50_vs_tmr1000.pooled.py](Code/Notebooks/o.vul_tmr50_vs_tmr1000.pooled.py). It supplies the pooled Figure 6 summaries and supporting Figures S5, S13 and S14. [noise_positions_and_reads_octopus.py](Code/Notebooks/noise_positions_and_reads_octopus.py) supplies native noise summaries. Preserve the scientific notebook's final annotation flags, eligible transcripts and one-isoform entries added under its stated rules.

## Octopus single-cell

**BioProject PRJNA796958 / GEO GSE193622** supplies the octopus single-cell data. Upstream preprocessing was carried out by **David Gorelik** and is documented in [Kobi_octupus_RNA_editing](https://github.com/zbidav/Kobi_octupus_RNA_editing). That URL is verified from the local upstream repository's Git remote.

This project begins from the upstream per-transcript BAMs, alignment summaries and read/cell annotations. Original commands retain the lab handoff paths:

| Population | Upstream handoff |
|---|---|
| All retained cells | `/private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample` and `Whole_collapsed_sample.regatedByChromBySampleSummary.tsv` |
| Neuronal | `splitSites/neuronal_reads/` and `neuronal_summary.tsv` |
| Non-neuronal | `splitSites/non_neuronal_reads/` and `non_neuronal_summary.tsv` |
| Read/cell metadata | `splitSites/Whole_collapsed_sample_extracted_read_data.withAnnotations.csv`; cell-class lookup [Code/Data/NeuronalityOfAnnotation.tsv](Code/Data/NeuronalityOfAnnotation.tsv) |

The `regatedByChromBySampleSummary` spelling above is the actual name used by the downstream commands. Consult David's repository for the upstream preprocessing and annotation handoff. Downstream here, reads are analyzed within the specified cell population to detect editing, translate possible proteins, select distinct isoforms and summarize gene-level diversity. A pooled cell-class result is not an individual-cell measurement.


Expression dependencies are distinct from isoform-count producers:

| Original heading | Expression role in the current workflow |
|---|---|
| All-cell 50 | The revised `ExpressionLevels.EntropyConsidered.MaxDistinct` assignments are required for the full-cohort association. The earlier `ExpressionLevels` cache is used only by historical per-cell diagnostics, with no retained manuscript consumer. |
| All-cell 1,000 | Historical non-entropy run; expression values are unused. The former existence filter was replaced with the identical 3,969-gene ordered distinct-solution cohort. |
| Neuronal 50 | Expression values are unused for Figure 6. The identical 3,593-gene solution cohort feeds 6,241 FixedTMR50 output genes after the established one-isoform fill. |
| Neuronal 1,000 | The original heading duplicates all-cell commands; it is not an independent neuronal expression workflow. |
| Non-neuronal 50 | Expression values are unused for Figure 6. The identical 931-gene solution cohort feeds 3,138 FixedTMR50 output genes after the established one-isoform fill. |
| Non-neuronal 1,000 | The original heading duplicates all-cell commands; it is not an independent non-neuronal expression workflow. |

Both FixedTMR50 producer tables were reproduced exactly from native solution/read summaries, including selected solution identities, gene order and counts. Final figure rendering was not rerun. No entropy-considered 1,000-read expression run is assumed.

### All cells: 50-read branch

The native downstream notebook loads the position, protein, selection and expression files from this branch.

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_bq 30 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table /private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample.regatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir /private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample \
--interfix_start ".comp" \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/pileup.26.1.25.out &
```

```bash
tmux new -s COMB18

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 80 --proc 10 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins.regular.29.1.25.log
```

The older all-cell 50-read expression recipe is omitted from this workflow. Its `ConcatAggGeneCellCouplesXPlusCopies.tsv` export and per-cell diagnostics have no retained figure, supplementary-table or numerical-text consumer. Their source remains in the notebook, and the historical recipe remains in [README.md](README.md). Figure 6 uses the separate native cell-class isoform-count producers.

### Revised all-cell assigned expression for the cell-type association analysis

This run prepares the entropy-considered maximum-selection assignments used by the cell-type association analysis below. Directory creation uses the standard `mkdir` command; all scientific arguments are unchanged. Expanding those assignments into annotated original-support rows remains a separate prerequisite.

```bash
OUT_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ExpressionLevels.EntropyConsidered.MaxDistinct

mkdir "$OUT_DIR"

# mkdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ExpressionLevels.EntropyConsidered.MaxDistinct

python \
Code/Simulations/prepare_fofns_for_expression.py \
--reads_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ReadsFiles \
--proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ProteinsFiles \
--distinct_proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins \
--proteins_postfix ".gz" \
--reads_postfix ".gz" \
--out_dir $OUT_DIR

# DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/DistinctProteinsForExpressionLevels.txt)
# ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/UniqueProteinsForExpressionLevels.txt)
# SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/ChromsNamesForExpressionLevels.txt)

nohup \
julia \
--project=. \
--threads 40 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $OUT_DIR/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $OUT_DIR/UniqueProteinsForExpressionLevels.txt \
--allreadsfilesfofn $OUT_DIR/ReadsForExpressionLevels.txt \
--samplenamesfile $OUT_DIR/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--considerentropy \
--innerthreadedassignment \
--onlymaxdistinct \
--outdir $OUT_DIR \
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/expressionlevels.EntropyConsidered.MaxDistinct.27.7.2026.out &
```

### All cells: 1,000-read branch

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_bq 30 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table /private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample.regatedByChromBySampleSummary.tsv \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir /private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample \
--interfix_start ".comp" \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/pileup.27.1.25.out &
```

```bash
tmux new -s COMB13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 80 --proc 10 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/DistinctProteins.regular.29.1.2025.log
```

Assigned-expression values are not required for this branch. The notebook now selects the same ordered genes using the native distinct-solution files; retain the pileup and isoform-selection producers above. The historical expression recipe remains in [README.md](README.md).

### Neuronal cells: 50-read branch

```bash

OUT_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.Neuronal
ALIGNMENT_DIR=/private10/Projects/David/Kobi_octupus/splitSites/neuronal_reads/

mkdir -p $OUT_DIR

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_bq 30 \
--out_dir $OUT_DIR \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table $ALIGNMENT_DIR/neuronal_summary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir $ALIGNMENT_DIR \
--interfix_start ".comp" \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--min_mapped_reads_per_position 0 \
> $OUT_DIR/pileup.10.4.25.out &
```

```bash
tmux new -s COMB18

COMB

BASE_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.Neuronal
OUT_DIR=$BASE_DIR/DistinctProteins

INFILES=$(echo $BASE_DIR/ProteinsFiles/*.unique_proteins.csv.gz)

mkdir $OUT_DIR

julia \
--project=. \
--threads 60 --proc 8 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir $OUT_DIR \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee $OUT_DIR/DistinctProteins.regular.10.4.25.log
```

Assigned-expression values are not required for this branch. The notebook now selects the same ordered genes using the native distinct-solution files; retain the pileup and isoform-selection producers above. The historical expression recipe remains in [README.md](README.md).

### Non-neuronal cells: 50-read branch

```bash

OUT_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.NonNeuronal
ALIGNMENT_DIR=/private10/Projects/David/Kobi_octupus/splitSites/non_neuronal_reads

mkdir -p $OUT_DIR

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_bq 30 \
--out_dir $OUT_DIR \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table $ALIGNMENT_DIR/non_neuronal_summary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir $ALIGNMENT_DIR \
--interfix_start ".comp" \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--min_mapped_reads_per_position 0 \
> $OUT_DIR/pileup.10.4.25.out &
```

```bash
tmux new -s COMB16

COMB

BASE_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.NonNeuronal
OUT_DIR=$BASE_DIR/DistinctProteins

INFILES=$(echo $BASE_DIR/ProteinsFiles/*.unique_proteins.csv.gz)

mkdir $OUT_DIR

julia \
--project=. \
--threads 60 --proc 8 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir $OUT_DIR \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee $OUT_DIR/DistinctProteins.regular.10.4.25.log
```

Assigned-expression values are not required for this branch. The notebook now selects the same ordered genes using the native distinct-solution files; retain the pileup and isoform-selection producers above. The historical expression recipe remains in [README.md](README.md).

The main Figure 6 cell-class producer is [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py), including its `FixedTMR50` export. The original README's “1,000 — neuronal” and “1,000 — non-neuronal” sections repeat the all-cell commands without cell-class-specific paths; they are not additional cell-class workflows and are not repeated here.

The separate expanded assigned-expression analysis in [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py) supplies the isoform–cell-type association analysis described under Specific analyses. It requires the revised assignments above and the original pseudo-read/cell-annotation mappings.


The manuscript’s 856-gene association is produced by the full-cohort Fisher/BH section of the [single-cell notebook](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py). Its execution table identifies selection/annotation preparation, original-support expansion and the common filters. Top100 analyses are separate diagnostics.

## Neural-expression annotations

**BioProject PRJNA285380** supplies the O. bimaculoides RNA-seq data used for neural-expression annotation transferred to O. vulgaris through orthology. Original FASTQs belong under `O.bim/Data/PRJNA285380/Raw`, and the native sample/tissue mapping is [Code/NeuralTranscripts/SraRunTable.txt](Code/NeuralTranscripts/SraRunTable.txt).

That table contains **18 runs**. The original README says “12 RNA-seq samples”, while the local manuscript describes 12 tissues; the relationship between that wording and the 18 run records is **TO VERIFY**. Use the actual run/tissue mapping and retained quantifications, rather than imposing an undocumented 12-run subset.

### Prepare short reads

```bash
nohup \
python Code/prepare_data.py \
--in_dir O.bim/Data/PRJNA285380/Raw/ \
--out_dir O.bim/Data/PRJNA285380 \
illumina \
--min_qual_mean 30 \
> O.bim/Data/PRJNA285380/prepare_data.out &
```

### Quantify expression and transfer neural annotations

The May 2023 processing records identify **Salmon 1.7.0** and **OrthoFinder 2.5.5**. The pipeline creates the Salmon index, quantifies O. bimaculoides transcripts, assigns neural-expression classes and maps them to O. vulgaris orthologs.

```bash
nohup \
python Code/neural_transcripts.py \
--obim_trinity_file O.bim/Annotations/orfs_bim.fa \
--ovul_trinity_file O.vulgaris/Annotations/orfs_oct.fa \
--obim_salmon_out_dir O.bim/Expression \
--obim_in_fastq_dir O.bim/Data/PRJNA285380/TrimmedWoDup \
--orthofinder_out_dir O.vulgaris/Annotations/OrthoFinderAgainstObim \
--ovul_neural_transcripts_file O.vulgaris/Annotations/NeuralVsNonNeuralExpression.BySalmonAndOrthoFinder.tsv \
> O.vulgaris/Annotations/neural_transcripts.out &
```

The generated output is `O.vulgaris/Annotations/NeuralVsNonNeuralExpression.BySalmonAndOrthoFinder.tsv`. Current [Code/neural_transcripts.py](Code/neural_transcripts.py) includes the corrected tissue labels and ortholog identifier handling. Figure S14 uses the corrected eligible-gene population and its final table handoff; an older similarly named annotation CSV is not an interchangeable replacement.

## Simulations and other processing


The July 2026 PB1 and Illumina complete-noise processing explicitly used **samtools 1.18 / HTSlib 1.17**. These commands use the **Jan2025** reference and annotations, whereas the original PB1/Illumina analyses above use the original reference. Preserve both the version pin and the reference selection; numerical differences cannot be attributed to software version alone from these records.

#### PB1 complete-noise positions

```bash
SAMTOOLS_118=/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools
OUTDIR=D.pealeii/MpileupAndTranscripts/CompleteNoisePositions

"$SAMTOOLS_118" --version

mkdir -p "$OUTDIR"

nohup \
python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--samtools_path "$SAMTOOLS_118" \
--out_dir "$OUTDIR" \
--processes 2 \
--threads 15 \
--gz_compression \
--keep_pileup_files \
--override_existing_pileup_files \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.LongReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
--keep_non_refbase_noisy_positions \
--stop_after_making_noise_positions \
> "$OUTDIR/pileup_with_subparsers.OnlyNoisePositions.PB1.samtools118.21.7.2026.out" 2>&1 &
```

#### Illumina complete-noise positions

```bash
SAMTOOLS_118=/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools
OUTDIR=D.pealeii/MpileupAndTranscripts/CompleteNoisePositions

"$SAMTOOLS_118" --version

mkdir -p "$OUTDIR"

nohup \
python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--include_flags 3 \
--exclude_flags 2304 \
--parity PE \
--min_bq 30 \
--samtools_path "$SAMTOOLS_118" \
--out_dir "$OUTDIR" \
--processes 10 \
--threads 10 \
--gz_compression \
--keep_pileup_files \
--override_existing_pileup_files \
directed_sequencing_data \
--data_table D.pealeii/Alignment/Illumina/reads.ByChrom/data_table.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
--keep_non_refbase_noisy_positions \
--stop_after_making_noise_positions \
> "$OUTDIR/pileup_with_subparsers.OnlyNoisePositions.Illumina.samtools118..21.7.2026.out" 2>&1 &
```

The explicit executable path and the doubled dot in the Illumina log filename are retained from the original command. Both blocks use `--override_existing_pileup_files`; their output directory is the named complete-noise branch.

[noise_positions_and_reads_squid.py](Code/Notebooks/noise_positions_and_reads_squid.py) joins the retained noise positions to original read mappings. The squid and pooled-octopus notebook sections produce the 12-mismatch position tables under `Code/Notebooks/12mm/v2/positions/`. The “12-mismatch analysis” section of [various_plots.BQ30.py](Code/Notebooks/various_plots.BQ30.py) produces sample/panel counts consumed by Figures S2 and S15. The relevant unit is a mismatch site within its sample/transcript population; preserve each denominator and the treatment of undefined proportions.



### Empirical graph assessment — Figure S16

Generate complete and partially observed simulated protein populations with the original input settings and random seeds in [Code/empirical_graph_assessment.py](Code/empirical_graph_assessment.py).

```bash
COMB

mkdir -p Simulations/GraphAssessment

# python \
# Code/parse_known_sites.py \
# --editing_sites_excel_file O.vulgaris/Annotations/editing_sites.xlsx \
# --sheet_name "O.vul Editing sites" \
# --out_csv_file Simulations/GraphAssessment/O.vul.EditingSites.csv \
# --out_bed_file Simulations/GraphAssessment/O.vul.EditingSites.bed

# nohup python Code/pileup_with_subparsers.py \
# --transcriptome D.pealeii/Annotations/orfs_squ.fa \
# --known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
# --parity SE \
# --out_dir Simulations/GraphAssessment \
# simulate_complete_and_corresponding_partially_unkown_data \
# directed_sequencing_data \
# --data_table Code/Alignment/DataTable.Squid.LongReads.csv \
# > D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/pileup.out &

nohup \
python Code/empirical_graph_assessment.py \
> Simulations/GraphAssessment/log.11.12.24.out &
```

The retained notebook reads the consistent-fraction selections from `Simulations/GraphAssessment/SameFrac`:

```bash
tmux a -t comb18

mkdir Simulations/GraphAssessment/SameFrac

# sed 's|GraphAssessment|GraphAssessment/SameFrac|g' Simulations/GraphAssessment/FalsePositivesOutFiles.fofn


INFILES=$(echo Simulations/GraphAssessment/*.UniqueProteins.tsv.gz)

julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .UniqueProteins.tsv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir Simulations/GraphAssessment/SameFrac \
--fracstep 0.2 \
--fracrepetitions 1 \
--consistentfracsampling \
--algrepetitions 8 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee Simulations/GraphAssessment/SameFrac/DistinctProteins.regular.21.1.25.log
```

The false-positive comparison uses the generated file-of-filenames inputs:

```bash
tmux a -t comb18

COMPLETE_INFILES=$(cat Simulations/GraphAssessment/CompleteInfiles.fofn)
ERRORED_NA_INFILES=$(cat Simulations/GraphAssessment/ErroredNAFiles.fofn)
OUT_FILES=$(cat Simulations/GraphAssessment/FalsePositivesOutFiles.fofn)

julia \
--project=. \
--threads 40 \
Code/Simulations/mis_5_assessment.jl \
--complete_infiles $COMPLETE_INFILES \
--errored_na_infiles $ERRORED_NA_INFILES \
--out_files $OUT_FILES \
2>&1 | tee Simulations/GraphAssessment/mis_5_assessment.24.12.2024.log
```

Continue with [Code/Notebooks/empirical_graph_assessment.py](Code/Notebooks/empirical_graph_assessment.py), which summarizes recovered and false-positive isoforms for Figure S16. These simulations assess the heuristic on specified populations; they do not establish global optimality of every manuscript solution.

### Editing simulation — Figure S19

The original source is [Code/Simulations/unifinished_editing_isoforms.jl](Code/Simulations/unifinished_editing_isoforms.jl), retaining the author's filename spelling. The original README names this script but supplies no execution command. The retained simulation used Julia 1.7.2, one thread, seed 1892, N = 10,000 reads, K = 10 sites and the nine ordered editing/removal probability pairs from {0.001, 0.002, 0.003}. Each panel contains 45 unique off-diagonal correlations and 15 right-open histogram bins; the upper edge uses `nextfloat(maximum)` to include the maximum. The scoped producer [manuscript_s19.jl](Code/Simulations/manuscript_s19.jl) preserves the original loop and random-draw order. With Julia 1.7.2 and StatsBase 0.34.3 it reproduced all 900 matrix values, 405 displayed values and 135 histogram rows exactly. The final renderer reads `s19_correlations.tsv` and `s19_histograms.tsv`; it does not rerun the simulation. Run from the repository root with a new destination:

```bash
OPENBLAS_NUM_THREADS=1 julia --project=. --startup-file=no --compiled-modules=no --threads=1 \
  Code/Simulations/manuscript_s19.jl Code/Notebooks/figure_sources/new_run/FigureS19 \
  Code/Notebooks/manuscript_intermediates/new_run/s19
```

### Chimeric-read analysis — Figure S20

The retained squid analysis uses original unique RNA patterns and selected assigned-expression tables. The original command for the unique-read implementation is:

```bash
nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/chimeric_reads_origin_squid_2.jl \
> /private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis2/chimeric_reads_origin_squid_2.18.5.2026.out &
```

The scientific notebook is [visualize_chimeric_reads_origin_squid_2.py](Code/Notebooks/visualize_chimeric_reads_origin_squid_2.py). The retained comparison uses common/rare fractions X = 0.01 and Y = 0.30, with original parent/read identities and the selected assignment order. The corrected producer [manuscript_s20.py](Code/Notebooks/manuscript_s20.py) joins original unique RNA patterns to the fixed selected proteins, preserving original multiplicities, ties-to-even rounding, prefix/suffix witnesses and the rule that any eligible original RNA can support a rare-protein chimera call. All 21 native RNA matrices and classification/witness tables matched exactly; the 63 expression rows, two box summaries and 42 chimera rows also matched exactly. Expression normalization occurs in native source-row order **before sorting**, then aligns to the canonical selected IDs; moving the normalization after sorting changes last-bit floating-point values. The original Julia command is historical provenance for this analysis; use the corrected producer below for the accepted final handoff:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s20.py \
  --assignments "$MANUSCRIPT_ASSIGNMENTS" \
  --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS20" --work "$MANUSCRIPT_WORK/s20"
```


Annotation-based recoding potential is produced by [recoding_sites_in_squid.py](Code/Notebooks/recoding_sites_in_squid.py), using the Jan2025 annotation. Its input → classification → triplet-summary → export sequence is documented in that notebook.

# Scientific notebooks and workflow

The scientific progression is **editing positions → reads and unique RNA patterns →
proteins and unique proteins → distinct isoform solutions → assigned expression →
scientific summaries and exports**. Original direct support counts the original reads
supporting an isoform; reassigned support is additional weighted expression mass.
Selection identifies distinct isoforms, while assigned expression supplies their
relative abundance. Historical exploratory figures are retained in the notebooks;
`final_manuscript_figures` is the publication renderer.

Each relevant notebook begins with **Reproducing manuscript results**. Its execution
table identifies setup, nonadjacent prerequisites, branch-specific inputs and exports.
Do not use Run All. Retained native exports are valid inputs only for the documented
cohort/version/schema; a matching filename alone does not establish validity.

## Reconstruct the notebook representations

`.ipynb` files and their saved outputs are ignored by Git. The trackable sources are
`.py:percent`, paired by [jupytext.toml](Code/Notebooks/jupytext.toml). In a fresh clone:

```bash
cd Code/Notebooks
jupytext --to ipynb \
  squids_pacbio_5.BQ30.py \
  squids_illumina_5.py \
  squids_pacbio_5.additional_UMIs.MergedSamples.py \
  squids_pacbio_5.additional_UMIs.MergedSamples.ReadsWithRecognizableBardoces.py \
  squids_pacbio_5.additional_UMIs.MergedSamples.DeduppedReads.py \
  o.vul_tmr50_vs_tmr1000.pooled.py \
  o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py \
  o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py \
  various_plots.BQ30.py \
  noise_positions_and_reads_squid.py \
  noise_positions_and_reads_octopus.py \
  empirical_graph_assessment.py \
  visualize_chimeric_reads_origin_squid_2.py \
  recoding_sites_in_squid.py \
  final_manuscript_figures.py
```

This reconstructs all manuscript-relevant notebook representations without execution. Resolve local edits before synchronizing an existing pair. A percent notebook is not necessarily an executable top-to-bottom Python script. Original legacy exploratory syntax and pre-existing PB1 entropy cell-boundary differences are documented validation limits, not a reason to execute the full notebook. Guide IDs refer to the saved author notebooks; reconstructed notebooks may assign new IDs, so use the accompanying headings and source comments.

# Generating manuscript figures and tables

Choose one empty table root and a separate empty figure-output root. All producers
and the renderer use `MANUSCRIPT_FIGURE_SOURCES` through `manuscript_paths.py`.
Assignment caches and RNA matrices belong in `MANUSCRIPT_WORK`, outside figure_sources.
The default table location is `Code/Notebooks/figure_sources/new_run`; existing outputs
are not a cache-validity signal. A new run must not overwrite the approved collection.

```bash
# Run from the repository root; choose a unique RUN_ID.
RUN_ID=your_run_id
export MANUSCRIPT_DATA_ROOT="$PWD"  # or the repository containing retained native data
export MANUSCRIPT_FIGURE_SOURCES="$PWD/Code/Notebooks/figure_sources/runs/$RUN_ID"
export MANUSCRIPT_FINAL_FIGURES="$PWD/Code/Notebooks/final_figures/runs/$RUN_ID"
export MANUSCRIPT_WORK="$PWD/Code/Notebooks/manuscript_intermediates/$RUN_ID"
export MANUSCRIPT_ASSIGNMENTS="$MANUSCRIPT_WORK/assignments"
export MANUSCRIPT_FONTS=/path/to/permitted/OpenSans/fonts
export MANUSCRIPT_ARIAL=/path/to/permitted/Arial.TTF
export NUMBA_CACHE_DIR="$MANUSCRIPT_WORK/numba_cache"
export MANUSCRIPT_FIGURE2_ASSETS=/path/to/approved/Figure2/assets
# Do not precreate producer output subdirectories; the producers guard creation.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_inputs.py \
  --project-root "$MANUSCRIPT_DATA_ROOT" --output "$MANUSCRIPT_ASSIGNMENTS"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_expression_tables.py \
  --assignments "$MANUSCRIPT_ASSIGNMENTS" --output "$MANUSCRIPT_FIGURE_SOURCES/shared_expression"
```

Start Jupyter from the configured shell, or set the same variables in the notebook
before importing the manuscript modules. A running kernel does not inherit later
shell-variable changes. Select the intended Python environment explicitly.

The fixed selection files under `manuscript_input_data/` are source configuration,
not generated plotting data. The input module authenticates the native files,
selected parameters, cardinalities and ordered protein identities. It reconstructs
assigned-expression tables without selecting new solutions or reassigning expression.

`various_plots.BQ30` is an optional cross-platform **table** producer. Pure table
adapters in `manuscript_plot_tables.py` accept explicit native input and new output
paths; their notebook execution table links the scientific producers. Use the same
root for their destinations. From `Code/Notebooks`, the retained-native-table
adapters can be called directly (use one export route per table; each output must be absent):

```python
from pathlib import Path
from manuscript_paths import FIGURE_SOURCES, ASSIGNMENTS_DIR
from manuscript_inputs import PROJECT
from manuscript_plot_tables import (
    figure2_table, figure3_tables, figure6_tables, export_native_s3, figure_s3_tables,
    figure_s1_tables, figure_s6_tables, figure_s10_tables, figure_s12_tables,
    figure_s13_tables, figure_s16_table, dependency_matrix_tables, assemble_fisher_tables,
)
native_notebooks = PROJECT / 'Code/Notebooks'
work = ASSIGNMENTS_DIR.parent
figure2_table(PROJECT, FIGURE_SOURCES / 'Figure2/PCLO_top100_unique_reads.tsv')
figure3_tables(native_notebooks, ASSIGNMENTS_DIR, FIGURE_SOURCES / 'Figure3', PROJECT)
figure6_tables(native_notebooks, FIGURE_SOURCES / 'Figure6')
figure_s1_tables(native_notebooks, FIGURE_SOURCES / 'FigureS1')
export_native_s3(PROJECT, work / 's3')
figure_s3_tables(work / 's3', FIGURE_SOURCES / 'FigureS3')
figure_s6_tables(native_notebooks, ASSIGNMENTS_DIR, FIGURE_SOURCES / 'FigureS6', PROJECT)
figure_s10_tables(native_notebooks, ASSIGNMENTS_DIR, FIGURE_SOURCES / 'FigureS10')
figure_s12_tables(native_notebooks, FIGURE_SOURCES / 'FigureS12')
figure_s13_tables(native_notebooks, PROJECT, FIGURE_SOURCES / 'FigureS13')
# Run the documented BY, S18 and S16 notebook branches before these handoffs:
dependency_matrix_tables(native_notebooks, work / 'by', FIGURE_SOURCES)
assemble_fisher_tables(work / 's18', FIGURE_SOURCES / 'FigureS18')
(FIGURE_SOURCES / 'FigureS16').mkdir(parents=True, exist_ok=True)
figure_s16_table(work / 's16', FIGURE_SOURCES / 'FigureS16/figure_s16_summaries.tsv')
```

## Producer-to-renderer navigation

`F` below denotes the selected `MANUSCRIPT_FIGURE_SOURCES` root. Each computational
row is consumed by the corresponding named section of `final_manuscript_figures.py`.
See each scientific notebook's execution table for prerequisites and native producers.

| Figure/table | Scientific producer | Generated handoff under F | Final consumer |
|---|---|---|---|
| 1 | Author manual artwork | External asset, not generated data | Author assembly |
| 2 | PB1 unique RNA; figure2_table | Figure2/PCLO_top100_unique_reads.tsv | figure2-asset / figure2-assembly |
| 3 | PB1/Illumina; figure3_tables | Figure3/squid_original_figure3{a,b,c,d}_plot.tsv | figure3-* |
| 4 | expression_tables; PB3 batch-export-figure4b | shared_expression/figure4a_expression.tsv; Figure4/figure4b_recurrence.tsv | figure4-* |
| 5 | PB1/Illumina All-reads entropy; BY cells + dependency_matrix_tables | Figure5/{figure5c_entropy,figure5_pclo_matrices.BY}.tsv | figure5-* |
| 6 | Pooled/cell-class notebooks; figure6_tables | Figure6/figure6_{cohorts,distributions,mean_markers}.tsv | figure6-* |
| S1 | Native read/noise exports; figure_s1_tables; violin_inputs | FigureS1/figure_s1_* | figure-s1-* |
| S2 | manuscript_s2; mismatch producers | FigureS2/s2_*; shared_mismatches/* | figure-s2-* |
| S3 | Original-read recoding; export_native_s3 → figure_s3_tables | FigureS3/figure_s3_{counts,statistics}.tsv | figure-s3-* |
| S4 | PB3 batch-export-figure-s4 | FigureS4/figure_s4_plot.tsv | figure-s4-* |
| S5 | Pooled batch-export-figure-s5 | FigureS5/figure_s5_plot.tsv | figure-s5-* |
| S6 | Native depth/AA/dispersion; figure_s6_tables | FigureS6/figure_s6* | figure-s6-* |
| S7 | manuscript_s7 matrices → scores → combine | FigureS7/s7_scores.tsv | s7-complete-updated-20261116-* |
| S8/S9 | manuscript_expression_tables | shared_expression/figure_s{8,9}* | figure-s8 / figure-s9 |
| S10 | PB1 fractions; PB3 separate barcode/UMI exports; figure_s10_tables | FigureS10/figure_s10_* | figure-s10-* |
| S11 | PB1 BY + native MI; dependency_matrix_tables | FigureS11/figure_s11_matrices.BY.tsv | figure-s11-* |
| S12 | Native Corrs.PB1/IL1; figure_s12_tables | FigureS12/figure_s12_{pairs,boxes,outliers}.tsv | figure-s12-* |
| S13 | Pooled selection/read membership; figure_s13_tables | FigureS13/figure_s13_* | figure-s13-* |
| S14 | Pooled s14-native-handoff / batch-export-figure-s14 | FigureS14/figure_s14_{cohorts,curves,test}.tsv | figure-s14-* |
| S15 | Noise and VP mm12 cells | shared_mismatches/combined_s2_s15_mm12_* | s15-* |
| S16 | empirical_graph_assessment; figure_s16_table | FigureS16/figure_s16_summaries.tsv | figure-s16-* |
| S17 | expression_tables; violin_inputs.s17_statistics_geometry | shared_expression/figure_s17_observations.tsv; FigureS17/* | figure-s17-* |
| S18 | PB1/Illumina native S18 Fisher/BH cells; assemble_fisher_tables | FigureS18/figure_s18_fisher.tsv and summaries | figure-s18-* |
| S19 | Code/Simulations/manuscript_s19.jl | FigureS19/s19_{correlations,histograms}.tsv | figure-s19-* |
| S20 | manuscript_s20, original RNA and fixed assignments | FigureS20/figure_s20_{expression,boxes,chimeras}.tsv | figure-s20-* |
| Tables S1/S3 | manuscript_metrics; original PB1/PB3 native stages | New metrics/workbooks outside F | manuscript_workbooks.export_pacbio_metrics_workbook |
| Table S2 | Illumina original-read sample SD cell | New TableS2_recoding_sample_SD.tsv outside F | export_sample_sd_workbook |
| Table S4 | New S18 gene Fisher/BH summary | FigureS18/s18_fisher_gene_summary.tsv | export_fisher_workbook |

## Exact S7, S19 and S20 commands

Keep S7's seven initialized OpenBLAS threads and its internal one-thread fitting;
changing the distance reduction order can change the final score bits.

```bash
OPENBLAS_NUM_THREADS=7 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s7.py prepare \
  --project-root "$MANUSCRIPT_DATA_ROOT" --assignments "$MANUSCRIPT_ASSIGNMENTS" --output "$MANUSCRIPT_WORK/s7_matrices"
for branch in PacBio_ordinary Illumina_ordinary PacBio_dissimilar; do
  OPENBLAS_NUM_THREADS=7 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s7.py scores \
    --matrices "$MANUSCRIPT_WORK/s7_matrices" --branch "$branch" --output "$MANUSCRIPT_WORK/$branch.tsv"
done
mkdir -p "$MANUSCRIPT_FIGURE_SOURCES/FigureS7"
python Code/Notebooks/manuscript_s7.py combine --branches \
  "$MANUSCRIPT_WORK/PacBio_ordinary.tsv" "$MANUSCRIPT_WORK/Illumina_ordinary.tsv" "$MANUSCRIPT_WORK/PacBio_dissimilar.tsv" \
  --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS7/s7_scores.tsv"
OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 Code/Simulations/manuscript_s19.jl \
  "$MANUSCRIPT_FIGURE_SOURCES/FigureS19" "$MANUSCRIPT_WORK/s19"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s20.py \
  --project-root "$MANUSCRIPT_DATA_ROOT" --assignments "$MANUSCRIPT_ASSIGNMENTS" \
  --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS20" --work "$MANUSCRIPT_WORK/s20"
```

## Rendering and supplementary workbooks

Open the final notebook with `Code/Notebooks` as the working directory. Its execution
table specifies shared setup, selected section order and final crop. Set
`FIGURES_TO_BUILD` explicitly; Figure 3 additionally uses `RENDER_FIGURE3`. Run only
the selected input/draw/export cells. F5 and S7 have dedicated crops; F2 has explicit
assembly. Do not run archived comparisons or collection-refresh code. New canvases
and evidence remain under the selected final directory's `render_work`; final outputs
are in `svg`, `pdf` and `previews`. Existing publication files are guarded. S11/S18 retain their approved Matplotlib physical crop and Poppler 600-dpi rasterization; the other general page exports use the documented MuPDF route.

For S1/S2/S17, exact retained SVG contours/jitter are external graphical inputs from
the original plotting cells; numerical values come from native tables. Do not infer
numerical data from pixels or assume a different plotting runtime recreates those
contours. `manuscript_s2.py sites` projects the native position tables separately;
`tables` consumes those projections, the reference/known-site annotations and the
matching distribution SVG. The three projection branches are squid, octopus and
octopus_motif; their populations are intentionally distinct.

Supply the four permitted historical SVG files in an external `GEOMETRY` directory,
using the basenames below. Their SHA-256 identities are part of the display contract:

| Figure | SVG basename | SHA-256 |
|---|---|---|
| S1 | Per chrom noise levels - PacBio1 vs. Illumina vs. octopus pooled.svg | `93d6706671f4782ef3647132d2862fd8728e7c4698916519a5d484c9dc32c744` |
| S2 | Distribution of pooled editing levels according to their editing status in current vs. previous study - Squid.svg | `a4bfa94291863f51e0e5a234b22fe8b56a2c3c4c9f59ad50e9bbc0103fc49b6c` |
| S17 | Mean edited positions vs expression levels - PacBio1.svg | `0998b5b1538251b007aba8915805f5e25ce3bd9a44fb94ea47a8ba2b2d7f4dbe` |
| S17 | Mean edited positions vs expression levels - Illumina.svg | `50e4e39430b76cfb2867759b495e09b686fa75a47a35fdc932684abd991e4b15` |

From `Code/Notebooks`, after the native noise and expression exports:

```python
from manuscript_figures.violin_inputs import s1_noise_geometry, s17_statistics_geometry
geometry = Path('/path/to/permitted/historical/SVGs')
s1_noise_geometry(native_notebooks,
    geometry / 'Per chrom noise levels - PacBio1 vs. Illumina vs. octopus pooled.svg',
    FIGURE_SOURCES / 'FigureS1')
s17_statistics_geometry(FIGURE_SOURCES / 'shared_expression/figure_s17_observations.tsv',
    [geometry / f'Mean edited positions vs expression levels - {platform}.svg'
     for platform in ['PacBio1', 'Illumina']], FIGURE_SOURCES / 'FigureS17')
```

For S2, run from the repository root. The links collect this run's three distinct
native projections; they do not copy approved tables:

```bash
GEOMETRY=/path/to/permitted/historical/SVGs
for branch in squid octopus octopus_motif; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s2.py sites \
    --branch "$branch" --project-root "$MANUSCRIPT_DATA_ROOT" \
    --output "$MANUSCRIPT_WORK/s2_sites/$branch"
done
mkdir "$MANUSCRIPT_WORK/s2_sites/combined"
ln -s "$MANUSCRIPT_WORK"/s2_sites/{squid,octopus,octopus_motif}/*.tsv "$MANUSCRIPT_WORK/s2_sites/combined/"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s2.py tables \
  --sites "$MANUSCRIPT_WORK/s2_sites/combined" --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS2" \
  --project-root "$MANUSCRIPT_DATA_ROOT" \
  --distribution-svg "$GEOMETRY/Distribution of pooled editing levels according to their editing status in current vs. previous study - Squid.svg"
```

The approved workbook/template is read-only. `manuscript_workbooks.py` writes new
workbooks while preserving unrelated XML members, author metadata, labels, primers
and manual fields. S1/S3 use `manuscript_metrics.py --project-root ... --output NEW_DIR`;
S2 uses original-read sample SD (`ddof=1`); S4 uses the current S18 Fisher/BH summary.
Supply `review_inputs/9.11.2026/Supplementary_tables_12_revised_completed.xlsx` and the
explicit S2 template where required. No workbook is reconstructed by guessing manual
fields. Six accepted last-decimal S1/S3 differences must remain full precision; do not
round them to force binary equality.

# External assets and reproduction limitations

Figure 1 requires the author's final manual/BioRender PDF. Figure 2 requires the
hash-identified `BioRender_base.pdf` and `BioRender_reference.pdf` documented in
`manuscript_figures/docker_export.py`, plus matching permitted
Arial files. Set `MANUSCRIPT_ARIAL` to the permitted Arial TTF file and
`MANUSCRIPT_FIGURE2_ASSETS` to the two-PDF manual input directory. The `figure2-asset-*` cells regenerate panel b
from this run’s PCLO rank table before assembly. Keep these outside figure_sources. Source availability does not establish
redistribution permission; no proprietary font or manual-artwork binaries are shipped.

The existing rendering image is `combinatorics-figure2:2026-10-13`, ID
`sha256:e1fcbeb34fbe87a26e76d553e8987277b29e3792eba52b7555f96fc1d5a47c97`.
Its [Dockerfile](Code/Notebooks/manuscript_figures/Dockerfile) depends on the external
`python-stats:1.0` base. A registry-manifest request returned access denied /
authentication required; a publicly obtainable exact base has not been established.
Obtain the matching image from its maintainer or validate a rebuild before claiming
identical artwork. Containers use no network, read-only inputs, two CPUs and 4 GiB RAM.

Native sequencing/processed inputs, the Jan2025 reference extension, original
single-cell metadata and retained solution files are external scientific inputs.
The accepted ScanProsite response remains at its existing documented external path
`Code/Notebooks/figure_sources/scanprosite_accepted.json` (SHA-256
`c02bc0b3e0563aab6590a06efc405e1a6fae0f013885ec5e2d8cd82513ad4ed5`);
its original query date/database release is unknown. This historical annotation is
not part of the new generated-table layout and is not regenerated or moved.

## Validation boundary and known issues

The 23 November 2026 integration test reproduced 24 of 26 publication figures
from documented retained native inputs into new table and figure directories.
All 24 match the approved PNG pixels; scientific table comparisons retain only
the explicitly justified floating-point differences described in the producers.
This is not a fresh raw-sequencing pipeline run. Figure 1 requires the missing
manual source; S16 panel A is exact, but its complete native summary exceeded the
30-minute limit in both the one-worker and four-worker attempts.

Fresh isolated `pacbiocomb` installation and all six CLI version checks passed.
Two fresh `combinatorics2` attempts reached the five-minute dependency-solver cap;
fresh scientific imports and plotting remain unverified. Julia 1.7.2 instantiation
remains incomplete after download/time and host temporary-storage limits; the final
workspace-temporary-directory retry stopped when its disk-budget monitor timed out.
The successful scientific and rendering tests used the documented installed runtimes.
The exact Docker base also requires authorized acquisition. These results establish
retained-input reproduction in those runtimes, not a complete clean-environment build.

Tables S1–S4 were reconstructed in new workbooks. S1/S3 retain the six accepted
last-bit differences; S2 reproduces original-read sample SD, with the archived
histogram-order reduction difference documented; S4 matches every member of the
approved Fisher/BH workbook. Existing scientific artifacts and saved notebook
displays were preserved.

- **Confirmed-open:** S16 needs a longer bounded native-summary run or a separately
  validated aggregation improvement; its full table and figure remain unverified.
- **Confirmed-open:** acquisition/version provenance for the named external inputs,
  manual assets and exact Docker base remains necessary for independent reproduction.
- **Confirmed-open:** six historical workbook last-bit differences are accepted
  when native inputs/formulas agree; their original serialization operation is not
  independently established. They are not scientific corrections.
- **Confirmed-open:** saved heuristic selections are not proven globally maximum or
  necessarily maximal. Preserve canonical solutions; selection improvements require
  a separate scientific task.
- **Suspected/unverified:** representative-row protein metadata may require a
  field-specific aggregation review; an effect on manuscript statistics is unproven.
- **Fixed in accepted source:** categorical All-reads entropy, binary recoding,
  Pearson BY/Fisher families and corrected neural labels. The observed selected
  assigned-expression entropy and All-reads marginal entropy use different empirical
  constructions; their difference is not automatically exact total correlation.

For historical workflows, see [README.md](README.md). It is preserved for provenance.
