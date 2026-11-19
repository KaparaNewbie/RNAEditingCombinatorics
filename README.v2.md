> **A-to-I editing generates unparalleled complexity in the neural transcriptome of cephalopods**
>
> _Kobi Shapira, Ruti Balter, David Gorelik, Joshua J C Rosenthal, Erez Y. Levanon & Eli Eisenberg_

**Revised-manuscript workflow — 19 November 2026.** This document describes the manuscript analyses, their processing commands and the producer–consumer chain for figures and tables. Commands target the current tracked source where equivalent behavior is established; qualified historical steps and unresolved reproduction requirements are identified explicitly.

For additional historical processing commands, alternative analyses, and extended development records, see the original [README.md](README.md). That document preserves provenance; it is not an independently validated reproduction path.

The workflow is:

**Original reads and annotations → preprocessing and alignment → editing detection, translation and isoform selection → scientific analysis notebooks → generated scientific tables → `final_manuscript_figures` → publication figures.**

The primary squid analyses use the original targeted PacBio data (PB1), the original pooled Illumina data, and the additional PacBio biological replicates (PB3). Octopus Iso-Seq and single-cell data support the transcriptome-wide comparisons. The earlier PacBio Data 2 dataset is outside the revised manuscript workflow.

## Contents

- [Software and setup](#software-and-setup)
- [References and annotations](#references-and-annotations)
- [Original squid PacBio data — PB1](#original-squid-pacbio-data--pb1)
- [Original squid Illumina data](#original-squid-illumina-data)
- [Additional squid biological replicates — PB3](#additional-squid-biological-replicates--pb3)
- [Octopus Iso-Seq data](#octopus-iso-seq-data)
- [Octopus neural-expression annotations](#octopus-neural-expression-annotations)
- [Octopus single-cell data](#octopus-single-cell-data)
- [Complete-noise positions and mismatch summaries](#complete-noise-positions-and-mismatch-summaries)
- [Simulations and chimeric-read analysis](#simulations-and-chimeric-read-analysis)
- [Specific analyses](#specific-analyses)
- [Analysis notebooks and figure production](#analysis-notebooks-and-figure-production)
- [Known issues and reproduction requirements](#known-issues-and-reproduction-requirements)

## Software and setup

The study used Conda environments that evolved during the original analysis and subsequent revision. The latest specification, [combinatorics2_simplified.yml](combinatorics2_simplified.yml), is an appropriate starting point for a modern analysis environment. It is not an immutable historical lockfile, and exact numerical equivalence to every historical analysis has not been established. The `combinatorics2` environment was created from [combinatorics2.yml](combinatorics2.yml) in February 2025 and later updated using the simplified specification. Processing-stage versions are identified alongside the relevant commands below.

The original environment specifications are [combinatorics.yml](combinatorics.yml) and [pacbiocomb.yml](pacbiocomb.yml). Historical package/build lists are available in [combinatorics.explicit.txt](combinatorics.explicit.txt) and [pacbiocomb.explicit.txt](pacbiocomb.explicit.txt); these filename lists are not complete Conda URL lockfiles.

### Clone and choose the environment

Choose an existing parent directory before cloning:

```bash
PROJECT_DIR="/path/to/some/dir" # replace with your existing parent directory
cd "$PROJECT_DIR"
git clone https://github.com/KaparaNewbie/RNAEditingCombinatorics
```

For a new modern installation, the following setup command uses the latest specification. This is an installation starting point, not a replacement for the stage-specific historical environments.

```bash
cd RNAEditingCombinatorics
conda env create -f combinatorics2_simplified.yml
conda activate combinatorics2
```

The alternative setup commands below define the original environments and aliases for the current interactive shell. They are alternatives to the modern installation above. The original work used CentOS 7; the historical environment files and paths may require a compatible system.

```bash
cd RNAEditingCombinatorics

conda env create -f combinatorics.yml
conda activate combinatorics  # test successful creation

conda env create -f pacbiocomb.yml
conda activate pacbiocomb  # test successful creation
```

```bash
alias COMB="cd \"${PROJECT_DIR}/RNAEditingCombinatorics\" && conda activate combinatorics"
alias PBCOMB="cd \"${PROJECT_DIR}/RNAEditingCombinatorics\" && conda activate pacbiocomb"
```

Unless a block changes directory, processing commands are run from the repository root. `COMB` returns there and activates `combinatorics`; it does not select a particular historical version by itself. In particular, PB3 March 2026 processing requires the recorded samtools 1.22.1 environment even where an original block retains the `COMB` alias. Verify the executable selected for that stage before execution.

Commands retain scientific inputs, filters, numerical settings, selected solution filenames and log destinations. Interface corrections are used only where source comparison supports equivalent behavior. These are separate processing stages, not a script to run top to bottom. Many target large lab servers; original Illumina isoform/expression calculations can require approximately 2 TB of RAM. Timestamped files identify retained selections: a new inference run need not choose the same isoform IDs or tie order. Resolve wildcard inputs to the intended files before execution.

### Command compatibility

The command interfaces and relevant defaults were inspected against tracked source at `01ef67b90f30638aa195e81f307dcbaa422b0ea0`. This is static inspection, not a completed reproduction run.

| Stage | Current-source compatibility and required action |
|---|---|
| Read preparation | PB1 now supplies `pacbio` with the established CCS defaults explicitly. Illumina and polished Iso-Seq subcommands match the parser. Current setup resolves FastQC, MultiQC, PRINSEQ-lite, Trimmomatic and its adapter file while constructing the parser, including for PacBio; provide a complete environment. Preserve the recorded sequencing-tool versions and adapter identity. |
| Alignment | **Blocked in current source:** `align.py` registers `--separate_by_chrom` twice in the PacBio parser, preventing parser construction for every subcommand. `--merge_by_chrom_bams` is absent although `pacbio_main` requires it; PB1 additionally lacks the subcommand and now-required known-sites argument. The four alignment recipes below remain explicitly historical until a compatible source revision or reviewed source repair is available. |
| Original Illumina pileup | **Blocked in current source:** the legacy `pileup.py` command parses, but passes nine arguments to `Pileup.mpileup.mpileup`, which now requires ten (`override_existing_pileup_files` is missing). A migration to the directed subparser must establish the same input population, default behavior and selected outputs; the historical recipe is retained with this qualification. |
| Directed and pooled pileup | Listed options match `pileup_with_subparsers.py`. Its default samtools path is `~/anaconda3/envs/combinatorics2/bin/samtools`, regardless of the active shell environment. PB1 explicitly selects the recorded 1.18 executable below; PB3 requires 1.22.1; July complete-noise commands explicitly select 1.18. For other stages, establish the matching executable/version before a fresh run. |
| Isoform selection and expression | Julia option names match the tracked parsers. Four older single-cell expression recipes omit the now-required read-table handoff (all-cell 50/1,000 and neuronal/non-neuronal 50); they remain qualified historical recipes. Use the tracked project, preserve `firstcolpos` (15 squid, 16 octopus), fractions, algorithms, grouping/matrix criteria and selected files. MIS sampling defaults to seed 1892; expression defaults to fraction 1.0 and Ascending/Descending. These defaults do not guarantee identical stochastic selections. |
| Annotation and auxiliary programs | ORF/BED conversion, known-site parsing and O. bim neural annotation match current interfaces. File-list preparation now requires `--reads_dir`; the revised maximum-selection single-cell recipe supplies it, while the four older expression recipes do not. Cloud transfers and simulation/chimera scripts contain site-specific settings; their input and handoff requirements are described in the corresponding sections. |

Current pipeline source uses Python 3.10+ syntax/type annotations; a historical environment name alone does not establish compatibility. The modern specification is the starting point for current Python source, while external tools still require the stage-specific versions.

### Julia

The Julia workflows use **Julia 1.7.2** and the tracked [Project.toml](Project.toml) / [Manifest.toml](Manifest.toml). The manifest specifies the dependency sources, including the custom BioAlignments source. Preserve that project environment when running the Julia commands.

```bash
COMB
mkdir Julia
cd Julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz
tar zxvf julia-1.7.2-linux-x86_64.tar.gz
```

Add the extracted Julia executable to `PATH`, then instantiate the project using the original interactive instructions. The `pkg>` lines are entered in Julia's package mode.

```bash
COMB

julia
pkg> activate .
pkg> instantiate
exit()
```

### Small demonstration

[Demo/DEMO.md](Demo/DEMO.md) describes a small demonstration after setup and acquisition of the original squid annotation. It illustrates pipeline capabilities; it is not a complete manuscript reproduction.

## References and annotations

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

## Original squid PacBio data — PB1

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
> D.pealeii/Alignment/BestN1/align.out &
```

**Historical recipe; currently blocked.** Adding `pacbio` alone is insufficient: the duplicate parser option, missing `merge_by_chrom_bams` destination and required known-sites argument must be resolved. The old and current core PB1 alignment use the CCS preset, best-n = 1, two processes and 25 threads, but the current wrapper adds separation/merging interfaces. No source repair or altered alignment command is supplied here.

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

## Original squid Illumina data

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

**Historical alignment recipe; the current `align.py` parser is blocked as described under Command compatibility.**

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

**Historical recipe; current `pileup.py` has the required-argument mismatch described under Command compatibility.**

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

## Additional squid biological replicates — PB3

| Public accession | Biological replicate | Original local input |
|---|---|---|
| SRR41297225 | JR1 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR1.r841072252411_C01.hifireads.bam` |
| SRR41297224 | JR2 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR2.r841072252411_C01.hifireads.bam` |
| SRR41297223 | JR3 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR3.r841072252411_C01.hifireads.bam` |

These are delivered **HiFi/CCS BAMs with UMIs**, not raw subreads. The accession-to-replicate assignments are author supplied and the exact filenames are verified locally. Provider records identify CCS **8.2.0**, lima **2.11.0**, pbtrim **1.1.0** and jasmine **2.3.0**. The local pipeline starts from those delivered files.

PB3 targets GRIA2 (`comp141693_c0_seq1`), ADAR1 (`comp134400_c0_seq1_extended`) and IQEC1 (`comp141565_c6_seq3`). Keep biological replicate identity through UMI processing, per-site comparisons and assigned-expression summaries.

### Align the delivered HiFi reads

March 2026 processing used local **pbmm2 26.1.99** and **samtools 1.22.1** with the extended Jan2025 reference. The original alignment command also separates and merges BAMs by transcript.

**Historical alignment recipe; the current `align.py` parser is blocked as described under Command compatibility.**

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

The March 2026 pileup also used **samtools 1.22.1**. Its data table must identify the all-read merged BAMs from that run. The current file at `D.pealeii/Alignment/AdditionalUMILongReads/DataTable.Squid.MergedAdditionalUMILongReads.csv` instead points into the recognizable-barcode output directory. Recovering the stage-specific table contents is necessary before reproducing this run; the current shared filename alone does not establish the March inputs.

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
--data_table D.pealeii/Alignment/AdditionalUMILongReads/DataTable.Squid.MergedAdditionalUMILongReads.csv \
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

## Octopus Iso-Seq data

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

**Historical alignment recipe; the current `align.py` parser is blocked as described under Command compatibility.**

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

## Octopus neural-expression annotations

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

## Octopus single-cell data

**BioProject PRJNA796958 / GEO GSE193622** supplies the octopus single-cell data. Upstream preprocessing was carried out by **David Gorelik** and is documented in [Kobi_octupus_RNA_editing](https://github.com/zbidav/Kobi_octupus_RNA_editing). That URL is verified from the local upstream repository's Git remote.

This project begins from the upstream per-transcript BAMs, alignment summaries and read/cell annotations. Original commands retain the lab handoff paths:

| Population | Upstream handoff |
|---|---|
| All retained cells | `/private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample` and `Whole_collapsed_sample.regatedByChromBySampleSummary.tsv` |
| Neuronal | `splitSites/neuronal_reads/` and `neuronal_summary.tsv` |
| Non-neuronal | `splitSites/non_neuronal_reads/` and `non_neuronal_summary.tsv` |
| Read/cell metadata | `splitSites/Whole_collapsed_sample_extracted_read_data.withAnnotations.csv`; cell-class lookup [Code/Data/NeuronalityOfAnnotation.tsv](Code/Data/NeuronalityOfAnnotation.tsv) |

The `regatedByChromBySampleSummary` spelling above is the actual name used by the downstream commands. Consult David's repository for the upstream preprocessing and annotation handoff. Downstream here, reads are analyzed within the specified cell population to detect editing, translate possible proteins, select distinct isoforms and summarize gene-level diversity. A pooled cell-class result is not an individual-cell measurement.

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

**Historical expression recipe; incomplete for current source.** File-list preparation now requires `--reads_dir` and the corresponding read suffix; expression assignment requires `--allreadsfilesfofn`. The matching native `ReadsFiles` population and equivalence to these retained assignments must be established before modernization. This recipe does not reproduce the manuscript's 856-gene result.

```bash
python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ProteinsFiles \
--distinct_proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins \
--proteins_postfix ".gz"

DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/ChromsNamesForExpressionLevels.txt)

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ExpressionLevels


nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/UniqueProteinsForExpressionLevels.txt \
--samplenamesfile O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ExpressionLevels \
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/expressionlevels.30.1.25.out &
```

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

**Historical expression recipe; incomplete for current source.** File-list preparation now requires `--reads_dir` and the corresponding read suffix; expression assignment requires `--allreadsfilesfofn`. The matching native `ReadsFiles` population and equivalence to these retained assignments must be established before modernization. This recipe does not reproduce the manuscript's 856-gene result.

```bash
SC_1000_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000

python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir $SC_1000_DIR/ProteinsFiles \
--distinct_proteins_dir $SC_1000_DIR/DistinctProteins \
--proteins_postfix ".gz" \
--out_dir $SC_1000_DIR

mkdir $SC_1000_DIR/ExpressionLevels

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $SC_1000_DIR/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $SC_1000_DIR/UniqueProteinsForExpressionLevels.txt \
--samplenamesfile $SC_1000_DIR/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir $SC_1000_DIR/ExpressionLevels \
> $SC_1000_DIR/expressionlevels.2.2.25.out &
```

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

**Historical expression recipe; incomplete for current source.** File-list preparation now requires `--reads_dir` and the corresponding read suffix; expression assignment requires `--allreadsfilesfofn`. The matching native `ReadsFiles` population and equivalence to these retained assignments must be established before modernization. This recipe does not reproduce the manuscript's 856-gene result.

```bash
tmux a -t COMB18

python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir  $BASE_DIR/ProteinsFiles \
--distinct_proteins_dir $BASE_DIR/DistinctProteins \
--proteins_postfix ".gz"

mkdir $BASE_DIR/ExpressionLevels

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $BASE_DIR/DistinctProteins/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $BASE_DIR/DistinctProteins/UniqueProteinsForExpressionLevels.txt \
--samplenamesfile $BASE_DIR/DistinctProteins/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir $BASE_DIR/ExpressionLevels \
> $BASE_DIR/expressionlevels.30.1.25.out &
```

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

**Historical expression recipe; incomplete for current source.** File-list preparation now requires `--reads_dir` and the corresponding read suffix; expression assignment requires `--allreadsfilesfofn`. The matching native `ReadsFiles` population and equivalence to these retained assignments must be established before modernization. This recipe does not reproduce the manuscript's 856-gene result.

```bash
tmux a -t COMB16

python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir  $BASE_DIR/ProteinsFiles \
--distinct_proteins_dir $BASE_DIR/DistinctProteins \
--proteins_postfix ".gz"

mkdir $BASE_DIR/ExpressionLevels

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $BASE_DIR/DistinctProteins/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $BASE_DIR/DistinctProteins/UniqueProteinsForExpressionLevels.txt \
--samplenamesfile $BASE_DIR/DistinctProteins/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir $BASE_DIR/ExpressionLevels \
> $BASE_DIR/expressionlevels.30.1.25.out &
```

The main Figure 6 cell-class producer is [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py), including its `FixedTMR50` export. The original README's “1,000 — neuronal” and “1,000 — non-neuronal” sections repeat the all-cell commands without cell-class-specific paths; they are not additional cell-class workflows and are not repeated here.

The separate expanded assigned-expression analysis in [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py) supplies the isoform–cell-type association analysis described under Specific analyses. It requires the revised assignments above and the original pseudo-read/cell-annotation mappings.

## Complete-noise positions and mismatch summaries

The July 2026 PB1 and Illumina complete-noise processing explicitly used **samtools 1.18 / HTSlib 1.17**. These commands use the **Jan2025** reference and annotations, whereas the original PB1/Illumina analyses above use the original reference. Preserve both the version pin and the reference selection; numerical differences cannot be attributed to software version alone from these records.

### PB1 complete-noise positions

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

### Illumina complete-noise positions

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

## Simulations and chimeric-read analysis

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

### Unfinished editing simulation — Figure S19

The original source is [Code/Simulations/unifinished_editing_isoforms.jl](Code/Simulations/unifinished_editing_isoforms.jl), retaining the author's filename spelling. The original README names this script but supplies no execution command. The exact retained run settings/seed and the producer of the final `s19_correlations.tsv` and `s19_histograms.tsv` handoff remain to be documented. The final renderer reads those tables; it does not rerun the simulation.

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

The scientific notebook is [visualize_chimeric_reads_origin_squid_2.py](Code/Notebooks/visualize_chimeric_reads_origin_squid_2.py). The retained comparison uses common/rare fractions X = 0.01 and Y = 0.30, with original parent/read identities and the selected assignment order. The current Figure S20 handoff also depends on the local corrected classification and table preparation code under `figure_rebuild_work/correction_2026-10-12/preparation/`. That handoff is not yet available through a complete tracked reproduction route; rerunning the old Julia command alone is not evidence that the current final plotting tables have been reproduced.

## Specific analyses

These analyses support manuscript text or a focused result. The named sections identify the intended calculation within notebooks that also contain exploratory work.

### Recoding potential at annotated squid editing sites

[recoding_sites_in_squid.py](Code/Notebooks/recoding_sites_in_squid.py), “Recoding potential at annotated squid RNA-editing sites”, uses the Jan2025 squid FASTA, ORF BED and editing-site BED. It enumerates nonempty subsets of annotated coding-strand A-to-G edits within each CDS triplet and classifies whether they preserve or change translation. The unit is a **CDS triplet within a transcript**; shared homology labels do not make transcripts one genomic gene. This is annotation-based potential, not an estimate of editing probabilities or evidence that the combinations co-occur.

Run its input setup, classification, annotated-site checks, per-triplet construction and export sequence. Generated outputs under `D.pealeii/Annotations/Jan2025/RecodingPotential/` include `squid_cds_codon_recoding.tsv.gz`, `squid_transcript_recoding_summary.tsv` and `squid_codon_editing_patterns.tsv`, together with the site inclusion/exclusion table. This analysis supports the manuscript text; it has no required final-figure rendering step.

### Association of isoform composition with cell type

The **manuscript result is the full-cohort analysis of 856 testable genes**, with 13 nominal P < 0.05 and one significant after BH correction (adjusted P ≈ 0.00856). The historical implementation is in commit `a8b75499c7c5386a79b289ac05d06789a4655ee5`, [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py), under “Distinct isoforms per annotation”. Read that version without changing the working tree using `git show a8b7549:Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py`.

The analysis pools totalRNA and nuclearRNA **originally supporting pseudo-reads**, retains genes with at least five inferred pooled isoforms, requires nonmissing SeuratCluster/Annotation, and excludes “not separated”, “unstable” and TBA1–TBA8. It counts isoform × cell-type observations within each gene, removes zero margins, and tests genes with at least two observed isoforms and two observed annotations. Reassigned expression mass is not an original-support count.

In that historical source, Jupytext cell ordinals 595, 605 and 610 filter/count/pivot; 622 defines the generalized Fisher test; 624 tests the complete cohort and applies BH over **all 856 testable genes**; 638 summarizes it. These are source ordinals including Markdown, not notebook execution counts or preserved cell IDs. The Monte Carlo settings are 100,000 resamples, batch 10,000 and seed `(20260730 + CRC32(Chrom)) mod 2^32`; BH uses alpha 0.05. Recorded medians are six original pseudo-reads, three observed isoforms and five cell-type annotations **per testable gene**.

The **additional analyses** currently in that notebook are `Top100ByCellTypes` and `Top100ByReads`. They rank genes by annotated cell-type count or all loaded read-table rows, break ties by Chrom, test each gene in their union once, and apply BH separately within each selection. They are not the manuscript's 856-gene result and are not independent validations of it.

**Required reproduction path:** native `SC.TotalCoverage50` read/protein/selection tables and retained entropy-considered assignments → annotated expanded original-support rows → the full-cohort historical cells above → per-gene contingency tables and Fisher/BH results → manuscript text. The active loader expects `ExpressionLevels.EntropyConsidered.MaxDistinct.MaxExpressionLevels/<chrom>.gz`; [parallel_max_sc_expression_worker.py](Code/Notebooks/parallel_max_sc_expression_worker.py) contains expansion functions. A documented, tracked route must preserve the retained solution keys, raw-read annotation joins, cache/input identities and full-cohort test family. Historical source is available in Git and an existing isolated reconstruction supports the November manuscript result, but its complete preparation/export route is not yet integrated into the current tracked notebook workflow. Running the current Top100 section does not fill this gap. No association analysis was implemented or rerun for this README revision.

### Domain enrichment

The domain sections of [squids_pacbio_5.BQ30.py](Code/Notebooks/squids_pacbio_5.BQ30.py) use the accepted union of protein-domain annotations and eligible reference-A positions. The `revision-20261109-domains-*` cells generate domain-site and Fisher summary tables without double-counting overlapping domain intervals. Editing-level comparisons are a separate analysis.

The required ScanProsite annotation values are currently read from saved notebook output. Source-only notebook reconstruction cannot recover those values; the accepted annotation input and acquisition procedure remain required for reproduction.

## Analysis notebooks and figure production

### Reconstruct notebooks with Jupytext

Large `.ipynb` notebooks are intentionally omitted from Git. Their tracked `.py:percent` representations are maintained with Jupytext using [Code/Notebooks/jupytext.toml](Code/Notebooks/jupytext.toml), which specifies `formats = "ipynb,py:percent"`.

For a fresh clone without local notebooks, reconstruct the required notebooks from their Python files. For example:

```bash
cd Code/Notebooks
jupytext --to ipynb squids_pacbio_5.BQ30.py squids_illumina_5.py final_manuscript_figures.py
```

Apply the same conversion to the other selected producer sources listed below. This creates notebook representations without running cells. Jupytext configuration is discovered from the notebook directory. Once paired files exist, normal Jupytext synchronization maintains the pair; resolve which representation contains the intended changes before syncing. Reconstruction does not supply saved outputs, external data or stable historical cell IDs that some existing consumers expect.

### Execution order and table contracts

1. Acquire the native reads, matching references and sample/cell metadata.
2. Complete the corresponding preparation, alignment, pileup, translation, heuristic selection and assigned-expression stages above.
3. Run the required sections of the original scientific notebooks, preserving their selected solution identities, populations and input versions. Generate their scientific tables first.
4. Assemble the final plotting-table handoffs from those native exports. Use the specified current branch where a later correction supersedes an earlier export.
5. Open `final_manuscript_figures.ipynb` in **combinatorics2**, with `Code/Notebooks` as its working directory, and render the required figures from verified tables and assets.

Source reconstruction is available, but the unresolved dependencies below prevent a claim of a complete tested clean-clone run. A saved table's existence is not sufficient evidence that it matches the current inputs or scientific method.

| Producer under `Code/Notebooks/` | Native inputs and generated scientific tables | Consumer / retained result |
|---|---|---|
| [squids_pacbio_5.BQ30.py](Code/Notebooks/squids_pacbio_5.BQ30.py) | PB1 native read/protein matrices, selected solutions and assigned expression; unique-RNA rank table, expression/recoding/correlation summaries | Figure 2b; squid Figures 3–5 and supporting panels |
| [squids_illumina_5.py](Code/Notebooks/squids_illumina_5.py) | Original 19-gene population and retained fixed expression assignments; expression/recoding/correlation summaries and Table S2 statistics | Figures 3–5, Table S2 and supporting panels |
| PB1/Illumina native recoding summaries and Table S2 export | Original-read-weighted recoding observations → `figure_s3_*` summaries and `TableS2_recoding_sample_SD.tsv` | Figure S3 and Table S2 |
| PB1/Illumina selected-expression sections and `various_plots.BQ30` | Isoform ranks and assigned-expression summaries → `figure4a_expression`, `figure_s8_s9_expression`, fitted-curve and marker tables | Figure 4a and Figures S8/S9 |
| PB1/Illumina original correlation-distance sections | Defined site-pair correlations grouped by distance → `figure_s12_{boxes,outliers,pairs}.tsv` | Figure S12; preserve its existing BY analysis |
| PB1 and Illumina, current categorical entropy sections | Native original read multiplicities per informative codon; selected isoform assigned expression → `figure_sources/figure5c_entropy.tsv` and method metadata | `figure5-load` / `figure5-entropy` / `figure5-export` → Figure 5c |
| PB1 and Illumina, `revision-20261109-by-*` | Defined unordered site pairs within gene/platform → BY-adjusted Pearson pair/matrix tables | Figure 5 and Figure S11; retain their MI values and site order |
| PB1 and Illumina, `revision-20261109-s18-*` | Selected common/rare isoforms → original directly supporting read A/G/NA counts → Fisher site/gene summaries | Figure S18 and Table S4; Figure S17 uses separate original-support summaries |
| [squids_pacbio_5.additional_UMIs.MergedSamples.py](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.py) | PB3 positions/read maps/replicate contributions → `figure4b_recurrence.tsv`, `source_s4.tsv` and native table metrics | Figure 4b, Figure S4, Tables S1/S3 |
| PB3 recognizable-barcode and deduplicated notebooks; PB1 fraction-0.1 branch | Branch-specific selections and assigned expression → paired/ranked expression summaries | `various_plots.BQ30` → Figure S10 |
| [o.vul_tmr50_vs_tmr1000.pooled.py](Code/Notebooks/o.vul_tmr50_vs_tmr1000.pooled.py) | Native pooled 50/1,000-read populations, read/protein tables and neural annotations → `MaxDistinctProtsForFig6.TMR50.Octopus.Pooled.csv`, 1,000-read counterpart and haplotype/coverage summaries | Figure 6 and Figures S5/S13/S14 |
| [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py) | Native cell-class position/protein/selection tables → `MaxDistinctProtsForFig6.{Neuronal,NonNeuronal}.FixedTMR50.Octopus.SC.csv` | `various_plots.BQ30` → Figure 6 |
| [noise_positions_and_reads_squid.py](Code/Notebooks/noise_positions_and_reads_squid.py), [noise_positions_and_reads_octopus.py](Code/Notebooks/noise_positions_and_reads_octopus.py), pooled octopus mismatch sections | Complete-noise/original read mappings → noise summaries and `12mm/v2/positions` transfers | `various_plots.BQ30` → Figures S1/S2/S15 |
| [various_plots.BQ30.py](Code/Notebooks/various_plots.BQ30.py) | Scientific exports from the platform-specific producers, fixed 80K-pair comparison and mismatch tables → cross-platform summaries and panel counts | Figures 3/6 and supporting Figures S1/S2/S6/S8/S9/S10/S15 |
| [empirical_graph_assessment.py](Code/Notebooks/empirical_graph_assessment.py) | Native simulated populations, SameFrac solutions and false-positive tables → S16 summaries | Figure S16 |
| [visualize_chimeric_reads_origin_squid_2.py](Code/Notebooks/visualize_chimeric_reads_origin_squid_2.py) | Native unique RNA patterns, selected assignments and corrected classification handoff → S20 expression/box/chimera tables | Figure S20 |
| [recoding_sites_in_squid.py](Code/Notebooks/recoding_sites_in_squid.py) and the single-cell association sections above | Annotated CDS triplet outcomes; per-gene original-support contingency tables | Specific manuscript text results |
| [final_manuscript_figures.py](Code/Notebooks/final_manuscript_figures.py) | Verified generated tables, existing plotting configuration, source PDFs and fonts | `final_figures/svg`, `final_figures/pdf`, `final_figures/previews` |

Important scientific definitions are preserved across these handoffs:

- **Figure 5c observed entropy:** Shannon entropy of the selected distinct isoforms' normalized assigned-expression probabilities, including eligible reassigned support.
- **Figure 5c hypothetical entropy:** sum of categorical amino-acid entropies per codon from **All reads** in the original PB1/Illumina populations. Count original informative-read multiplicities; exclude an ambiguous AA state at that codon, not automatically the whole read. The originally-supporting-read version is a diagnostic. The observed/hypothetical difference is not automatically exact total correlation or a strict maximum-entropy gap for the selected distribution.
- **Binary recoding for clustering:** [Code/EditingUtils/recoding.py](Code/EditingUtils/recoding.py) maps reference-only states to 0, nonempty AA sets excluding the reference to 1, and mixed-reference or missing states to NA. Preserve the existing imputation, retained columns, selected isoforms, order and random settings. This is distinct from categorical AA entropy.
- **Figure S18/Table S4:** two-sided Fisher tests use integer original-support A/G counts; per-group editing denominators are A+G, and BH correction is within gene over testable sites. **Pearson BY** applies to the defined unordered-pair family within gene/platform.
- **Replicate recurrence:** assigned mass is normalized within gene and biological replicate. An assigned-expression weight is not an original-read count. Keep the original read → protein → unique protein → selected isoform relationships.

### Final figures and supplementary tables

Recorded October/November figure work used Python **3.13.14**, NumPy **2.4.6**, pandas **3.0.3** and Matplotlib **3.11.0** in combinatorics2; notebook pairing records Jupytext **1.19.4**. These are dated execution/source records, not a promise about versions resolved by a new installation.

The final notebook consumes persisted scientific summaries; expensive inference belongs in the original scientific producers. Figure-specific plotting parameters belong in the figure notebook or its existing configuration. Generated TSV/CSV tables and SVG/PDF/PNG figures are **not required tracked source files**. Reproduction depends on their producer code, input definitions and documented handoffs.

The current renderer uses generated tables under `Code/Notebooks/figure_rebuild_work/`, particularly `correction_2026-10-12/inputs` and `phase4_batch/inputs`, plus the current correction-specific inputs. Examples include `figure5_pclo_matrices.BY.tsv`, `figure_s11_matrices.BY.tsv`, `combined_s2_s15_mm12_panel_counts.tsv`, `figure_s18_fisher.tsv` and the later corrected S14 tables. These locations remain active dependencies.

For Figure 5c, run the original PB1 and Illumina entropy/export sections sequentially after their selected-expression and original-read prerequisites. The combined table must contain the two PB1 and 19 Illumina cohorts and the matching method/input metadata before rendering. The required cohort and layout definitions in `figure_sources` are not currently tracked.

The canonical Figure S7 renderer consumes the corrected score handoff `s7_old_vs_updated_summary.tsv`. Its full score-generation route still depends on local preparation scripts, even though binary recoding itself is now in tracked source. Preserve the exact matrix, clustering and random-number settings.

Figure 2 combines approved manual/BioRender PDFs and the original unique-RNA rank panel. Two source PDFs are tracked under `Fig 2 inputs/`; the assembly handoff, configuration and fonts are additional prerequisites. The recorded assembly container is `combinatorics-figure2:2026-10-13` (Python 3.11.16, PyMuPDF 1.26.7). Its local instructions are in `figure_rebuild_work/figure2_assembly_2026-10-13/reports/REPRODUCE.txt`. Arial and Open Sans must be available through permitted acquisition; proprietary font files must not be redistributed. Figure 1 remains external/manual artwork.

Supplementary workbook generation also requires the author's approved workbook template and the native metric producers. The Table S2 producer uses the sample standard deviation (`ddof=1`); Table S4 uses the revised Fisher results. Later PB3 count-check tables are not substitutes for the complete Tables S1/S3 producers.

Select the intended figure sections explicitly: current configuration, render gates and old introductory notebook prose do not describe one universal execution path. If the legacy collection-refresh section is used, run the corrected S7 export **after** it. Figure 5 has a preservation guard, but a collection refresh is not a general guarantee against replacing other current artwork.

## Known issues and reproduction requirements

Static source inspection has not established a complete clean-clone reproduction. Generated tables and graphics need not be tracked; their required inputs and reproducible producers must be available.

| Category / status | Practical limitation |
|---|---|
| Execution blockers — confirmed-open | All four `align.py` recipes are blocked by the duplicate parser option and incomplete PacBio interface. Original Illumina `pileup.py` omits a now-required internal argument. Four older single-cell expression recipes lack the required read-table/file-list arguments. Resolve the exact source/interface compatibility before execution. |
| Required inputs — confirmed-open | Exact Jan2025 reference/annotation preparation, supplied polished octopus CCS, native sample/cell metadata, retained solution identities, the Illumina raw-directory handoff and PB3 stage-specific data tables remain required. Cloud scripts depend on the author's storage/configuration. |
| Required producers — confirmed-open | Integrate the full-cohort 856-gene preparation/test/export route. Complete tracked scientific-table/figure adapters, including `canonical_io`, Fisher/BY helpers, S7/S14/S19/S20 handoffs, supplementary-workbook assembly and Figure 2 assembly/assets. Domain annotation values currently held in saved notebook output also need an available input. See the figure dependencies above. |
| Versions and compatibility — confirmed-open | Preserve stage-specific samtools/PacBio/Illumina tools, adapter files, Julia 1.7.2 and the tracked custom BioAlignments project. Current Python source requires a compatible modern environment. Hard-coded lab paths and local assembly/font requirements need documented equivalents; a new environment solve is not a historical lockfile. |
| Scientific limitations — confirmed-open; effect on retained results unquantified | Descending heuristic solutions do not guarantee maximal extension or a global maximum. `expressionlevels.jl` also documents a limitation in choosing missing candidates. Preserve canonical selections; these observations alone do not quantify an error in a manuscript result. |
| Scientific concerns — suspected/unverified | Representative-row protein metadata may need field-specific aggregation review. Whether this affects a retained statistic has not been established. |
| Fixed in current source | Categorical All-reads entropy, binary recoding and the accepted Pearson BY/Fisher and neural-label corrections are present; the remaining production/assembly gaps above still apply. |

For help with the pipeline or file formats, open an issue in this repository.
