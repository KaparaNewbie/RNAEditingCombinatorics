> **A-to-I editing generates unparalleled complexity in the neural transcriptome of cephalopods**
>
> _Kobi Shapira, Ruti Balter, David Gorelik, Joshua J C Rosenthal, Erez Y. Levanon & Eli Eisenberg_

This repository contains the processing code, scientific notebooks and figure-generation workflow for the manuscript. The guide follows the analyses from sequencing inputs to manuscript figures and supplementary tables. Historical commands and development records are archived in [Logs.Final.md](Logs.Final.md).

The workflow is:

**Original reads and annotations → preprocessing and alignment → RNA editing detection, protein-sequence inference and isoform selection → scientific analysis notebooks → generated scientific tables → `final_manuscript_figures.ipynb` → publication figures.**

The primary squid analyses use the original targeted PacBio data (PB1), the original pooled Illumina data, and additional squid PacBio biological replicates (JR1–JR3). Octopus Iso-Seq and single-cell data support the transcriptome-wide comparisons.


# Contents
<br>
<table>
  <tr>
    <th>Section</th>
    <th>Subsections</th>
  </tr>
  <tr>
    <td><a href="#installation-and-software-requirements">Installation and software requirements</a></td>
    <td>—</td>
  </tr>
  <tr>
    <td><a href="#reference-data-and-annotations">Reference data and annotations</a></td>
    <td>
      <details>
        <summary>Show subsections</summary>
        <ul>
          <li><a href="#original-squid-reference-and-editing-annotations">Original squid reference and editing annotations</a></li>
          <li><a href="#extended-squid-reference-jan2025">Extended squid reference (Jan2025)</a></li>
          <li><a href="#octopus-reference-and-editing-annotations">Octopus reference and editing annotations</a></li>
          <li><a href="#octopus-bed-gff-and-gtf-for-the-single-cell-handoff">Octopus BED, GFF and GTF for the single-cell handoff</a></li>
          <li><a href="#o-bimaculoides-reference">O. bimaculoides reference</a></li>
        </ul>
      </details>
    </td>
  </tr>
  <tr>
    <td><a href="#analyses">Analyses</a></td>
    <td>
      <details>
        <summary>Show subsections</summary>
        <ul>
          <li><a href="#original-squid-pacbio-pb1">Original squid PacBio (PB1)</a></li>
          <li><a href="#squid-illumina">Squid Illumina</a></li>
          <li><a href="#additional-squid-pacbio-biological-replicates">Additional squid PacBio biological replicates</a></li>
          <li><a href="#octopus-iso-seq">Octopus Iso-Seq</a></li>
          <li><a href="#octopus-single-cell">Octopus single-cell</a></li>
          <li><a href="#neural-expression-annotations">Neural-expression annotations</a></li>
          <li><a href="#simulations-and-other-processing">Simulations and other processing</a></li>
        </ul>
      </details>
    </td>
  </tr>
  <tr>
    <td><a href="#scientific-notebooks-and-workflow">Scientific notebooks and workflow</a></td>
    <td>
      <details>
        <summary>Show subsections</summary>
        <ul>
          <li><a href="#reconstruct-the-notebook-representations">Reconstruct the notebook representations</a></li>
        </ul>
      </details>
    </td>
  </tr>
  <tr>
    <td><a href="#figures-and-supplementary-tables">Figures and supplementary tables</a></td>
    <td>—</td>
  </tr>
  <tr>
    <td><a href="#methodological-limitations-and-future-improvements">Methodological limitations and future improvements</a></td>
    <td>—</td>
  </tr>
</table>
<br>

# Installation and software requirements

Use separate environments for Python analysis, PacBio processing and Julia. Run commands from the repository root unless indicated otherwise. Install Julia separately using the [stage-specific version table](software_versions/historical_stages.2026-11-23.md); Docker is needed for the [figure-export steps](docs/figure_reproduction.md).

```bash
git clone https://github.com/KaparaNewbie/RNAEditingCombinatorics
cd RNAEditingCombinatorics
conda env create -f combinatorics2_simplified.yml
conda env create -f pacbiocomb.yml
conda activate combinatorics2
SAMTOOLS=$(command -v samtools)
MANUSCRIPT_REPO="$PWD"
alias COMB='cd "$MANUSCRIPT_REPO" && conda activate combinatorics2'
alias PBCOMB='cd "$MANUSCRIPT_REPO" && conda activate pacbiocomb'
# With the documented Julia version installed and on PATH:
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Use [combinatorics2_simplified.yml](combinatorics2_simplified.yml), [pacbiocomb.yml](pacbiocomb.yml), [Project.toml](Project.toml) and [Manifest.toml](Manifest.toml). Exact installed versions and historical requirements are documented in [software_versions/README.md](software_versions/README.md) and the [stage-specific version table](software_versions/historical_stages.2026-11-23.md); the current Conda environments do not necessarily supply historical tool versions.

| Environment | Purpose and requirements |
|---|---|
| combinatorics2 | Python processing, analyses, scientific notebooks and plotting. |
| pacbiocomb | PacBio command-line tools. |
| Julia | Independent-set selection, expression assignment and simulations. |

`SAMTOOLS` selects the active environment's executable; set it to an absolute executable path if needed. The `SAMTOOLS_118` and `SAMTOOLS_1221` overrides below select the versions documented in the [stage-specific table](software_versions/historical_stages.2026-11-23.md). For pileup stages whose historical version is not established, the generic `SAMTOOLS` path provides portability without claiming historical version equivalence; an alignment version alone does not establish the pileup version.

Figure export additionally requires Chrome for Kaleido, `rsvg-convert`, Poppler and permitted Open Sans/Arial fonts. Figure 2 assembly and PDF processing use the Docker image containing PyMuPDF; see [external assets](docs/figure_reproduction.md#special-external-assets-and-figure-specific-requirements).

Adjust paths to your system and use the reference and annotation versions specified for each analysis. Raw Illumina isoform selection and reassignment can require approximately 2 TB RAM; figure generation from processed tables is less demanding. Some simulation scripts also contain site-specific paths.

The complete demo was validated in the documented installed environments using four CPU cores and took approximately four minutes. See [Demo/DEMO.md](Demo/DEMO.md) for instructions and expected outputs.

Clean-install limitations are listed in [software_versions/README.md](software_versions/README.md); `combinatorics2` and Julia reproduction currently relies on the documented installed environments.

# Reference data and annotations

Use the reference and annotation versions named in each workflow. The original squid transcriptome and the extended Jan2025 transcriptome are distinct inputs.

| Reference or annotation | Local location and role |
|---|---|
| Original squid transcriptome | `D.pealeii/Annotations/orfs_squ.fa`; original PB1 and Illumina processing |
| Extended squid transcriptome | Supplied as [reference_data/Jan2025/orfs_squ.fa](reference_data/Jan2025/orfs_squ.fa); copy to `D.pealeii/Annotations/Jan2025/orfs_squ.fa` for additional squid PacBio biological replicates and July 2026 complete-noise processing |
| Squid editing sites | Generate the original `D.pea.EditingSites.csv` and `.bed` from the public Excel below; use these annotations with either squid FASTA. No validated known-site records are available for `comp134400_c0_seq1_extended`. |
| Octopus transcriptome and editing sites | `O.vulgaris/Annotations/orfs_oct.fa`, `O.vul.EditingSites.bed` / `.csv` |
| O. bimaculoides transcriptome | `O.bim/Annotations/orfs_bim.fa`; expression and orthology-based neural annotation |
| ORF intervals | Species-specific `orfs_*.bed`, generated by [Code/orfs_fasta_to_bed.py](Code/orfs_fasta_to_bed.py) from FASTA header coordinates |
| Single-cell transcript annotation | `O.vulgaris/Annotations/oct.bed` → `oct.gff` → `oct.gtf`, using the commands below; the exact reference/index handoff is also required from the upstream repository |
| O. bimaculoides GTF | `O.bim/Annotations/orfs_bim.gtf`; producer [prepare_salmon_gtf_input_from_transcriptome.py](Code/NeuralTranscripts/prepare_salmon_gtf_input_from_transcriptome.py). The current `neural_transcripts.py` route indexes the transcript FASTA directly. |

The original squid FASTA has SHA-256 `49f094d8f83aabb5961760353aca9d2d349c185a19a14451e78dd0c45f1c6e09` and is obtained using the download commands below.

## Original squid reference and editing annotations

Use the installed samtools executable selected by `SAMTOOLS` above to create the FASTA index; no additional Conda environment is needed for this step.

```bash
mkdir -p D.pealeii/Annotations
cd D.pealeii/Annotations

wget https://www.tau.ac.il/~elieis/squid/orfs_squ.fa
"$SAMTOOLS" faidx orfs_squ.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx
wget https://www.tau.ac.il/~elieis/squid/ST4.txt.gz
wget https://www.tau.ac.il/~elieis/squid/ST5.txt.gz
gunzip ST4.txt.gz
gunzip ST5.txt.gz

COMB  # return to the main project' dir
```

Generate the original squid CSV and six-column BED with the existing [known-sites parser](Code/parse_known_sites.py). The BED start is zero-based and its end is the original one-based editing location. `--ignore_score` retains the original `.` score convention. Preserve existing files: the command below generates outputs only when neither exists; otherwise, verify the saved files before reuse. Use a fresh output directory if an existing file fails verification.

```bash
if [ ! -e D.pealeii/Annotations/D.pea.EditingSites.csv ] &&
   [ ! -e D.pealeii/Annotations/D.pea.EditingSites.bed ]; then
  python Code/parse_known_sites.py \
    --editing_sites_excel_file D.pealeii/Annotations/editing_sites.xlsx \
    --sheet_name "D.pea Editing sites" \
    --out_csv_file D.pealeii/Annotations/D.pea.EditingSites.csv \
    --out_bed_file D.pealeii/Annotations/D.pea.EditingSites.bed \
    --ignore_score
fi
```

The public squid sheet and the original generated files contain 82,714 records, including 14 attached to the original partial `comp134400_c0_seq1`; none belongs to `comp134400_c0_seq1_extended`. Keep the original transcript identifiers and coordinates unchanged. Verify the files before using them:

```python
from pathlib import Path
import pandas as pd

annotations = Path("D.pealeii/Annotations")
sites = pd.read_excel(annotations / "editing_sites.xlsx", sheet_name="D.pea Editing sites")
csv = pd.read_csv(annotations / "D.pea.EditingSites.csv")
bed = pd.read_csv(annotations / "D.pea.EditingSites.bed", sep="\t")
extended_adar1 = "comp134400_c0_seq1_extended"
assert len(sites) == len(csv) == len(bed) == 82714
for frame, column in [(sites, "Trinity id"), (csv, "Trinity id"), (bed, "#TrinityID")]:
    assert not frame[column].eq(extended_adar1).any()
pd.testing.assert_frame_equal(sites, csv, check_dtype=False,
                              check_exact=False, rtol=1e-14, atol=1e-15)
assert bed["#TrinityID"].equals(sites["Trinity id"])
assert bed["Start"].eq(sites["Editing location (base1)"] - 1).all()
assert bed["End"].eq(sites["Editing location (base1)"]).all()
assert bed["SwissProtName"].fillna("").equals(sites["SwissProt name"].fillna(""))
assert bed["Score"].eq(".").all() and bed["Strand"].eq("+").all()
```

Extract ORF coordinates from the original transcriptome:

```bash
python Code/orfs_fasta_to_bed.py \
--in_fasta D.pealeii/Annotations/orfs_squ.fa \
--out_bed D.pealeii/Annotations/orfs_squ.bed
```

## Extended squid reference (Jan2025)

The supplied [Jan2025 FASTA](reference_data/Jan2025/orfs_squ.fa) includes the complete *D. pealeii* ADAR1 transcript `comp134400_c0_seq1_extended`, provided by Joshua Rosenthal. Its ORF is 3,741 nucleotides long. For the biological source, see [Vallecillo-Viejo et al. (2023), *Squid express conserved ADAR orthologs that possess novel features*](https://doi.org/10.3389/fgeed.2023.1181713). This repository supplies the exact frozen workflow FASTA, without requiring reconstruction from the publication.

Its SHA-256 is `b19dd2d1eb6aa9e8601fbc93b00145fe787b208bc7f5c7d119eaad08e97f1221`. From the repository root, copy it to the working location and create the FASTA index:

```bash
mkdir -p D.pealeii/Annotations/Jan2025
cp reference_data/Jan2025/orfs_squ.fa \
  D.pealeii/Annotations/Jan2025/orfs_squ.fa
"$SAMTOOLS" faidx D.pealeii/Annotations/Jan2025/orfs_squ.fa
```

Generate the ORF BED directly from the extended FASTA headers. Preserve the file if it already exists:

```bash
if [ ! -e D.pealeii/Annotations/Jan2025/orfs_squ.bed ]; then
  python Code/orfs_fasta_to_bed.py \
    --in_fasta D.pealeii/Annotations/Jan2025/orfs_squ.fa \
    --out_bed D.pealeii/Annotations/Jan2025/orfs_squ.bed
fi
```

Use the original public-Excel-derived `D.pealeii/Annotations/D.pea.EditingSites.bed` for known-site annotation with this FASTA, and the corresponding original `.csv` for known-site comparisons. The Jan2025 FASTA retains all 11,963 original transcript sequences unchanged and adds `comp134400_c0_seq1_extended`; retain the full extended sequence and the original annotations for all other identifiers. Do not transfer coordinates from the partial ADAR transcript to the extended transcript or reconstruct them by approximate mapping.

The historical `D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed` and `.csv` contain 35 manually constructed extended-ADAR1 records whose mapping and editing-site coordinates are unverified. These records are unsuitable for known-site annotation and must not enter newly recommended annotation files. Zero records for the extended identifier means that a validated known-site annotation is unavailable; keep ADAR1 previously reported editing sites as **N/A**, not zero, in Supplementary Table S3.

**Historical reproduction and corrected inputs.** Preserve historical annotations and generated outputs without overwriting them. For new processing, choose fresh output directories and update their downstream paths consistently. The processing examples below use the original known-site BED as the corrected annotation input; exact reproduction of retained historical outputs is a separate provenance task and does not validate their known-site labels. No scientific outputs have been regenerated with the corrected inputs. Legacy notebook paths still referring to the Jan2025 known-site CSV/BED require a separate input correction before new annotation-based analysis. In particular, the historical extended-ADAR1 results from [recoding_sites_in_squid.ipynb](Code/Notebooks/recoding_sites_in_squid.py) depend on the unverified records and must not be interpreted as validated recoding potential.

## Octopus reference and editing annotations

```bash
mkdir -p O.vulgaris/Annotations
cd O.vulgaris/Annotations
wget https://www.tau.ac.il/~elieis/squid/orfs_oct.fa
"$SAMTOOLS" faidx orfs_oct.fa

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

## Octopus BED, GFF and GTF for the single-cell handoff

Construct the transcript BED and convert it to GFF/GTF using the versioned AGAT container in the commands below; see the [software requirements](software_versions/README.md#additional-documented-tool-requirements).

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

## O. bimaculoides reference

```bash
mkdir -p O.bim/Annotations
cd O.bim/Annotations
wget https://www.tau.ac.il/~elieis/squid/orfs_bim.fa
"$SAMTOOLS" faidx orfs_bim.fa

COMB # back to the main directory of the project
```

# Analyses

## Original squid PacBio (PB1)

| Public accession | Manuscript target | Original local input |
|---|---|---|
| SRR28713506 | GRIA2 (called `GRIA` in older tables/commands) | `D.pealeii/Data/Raw/GRIA-CNS-RESUB.C0x1291.subreads.bam` |
| SRR28713505 | PCLO | `D.pealeii/Data/Raw/PCLO-CNS-RESUB.C0x1291.subreads.bam` |

These original PacBio subreads BAMs map to target regions `comp141693_c0_seq1` and `comp141882_c0_seq14` in [DataTable.Squid.LongReads.csv](Code/Alignment/DataTable.Squid.LongReads.csv).

### Prepare consensus reads

Generate consensus reads from the demultiplexed subreads with CCS, using the PB1 version in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md). The Python preparation wrapper requires the complete `combinatorics2` toolset, including Trimmomatic adapters, even for PacBio input.

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

Preparation creates BasicCCS, HiFi and deduplicated branches. The downstream alignment uses **BasicCCS**.

### Align

Align with pbmm2 using the PB1 version in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md). The resulting GRIA/PCLO `.aligned.sorted.bam` files are selected by the native long-read data table.

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

Alignment uses the CCS preset, best-n = 1, two processes and 25 threads, without chromosome separation or merging.

### Detect RNA editing and infer protein sequences

Use samtools as documented for PB1 in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md), the original squid reference, SE handling, minimum read quality 0.998 and base quality 30. Set `SAMTOOLS_118` to your executable; `--samtools_path` overrides the pileup wrapper's environment-specific default.

```bash
COMB
SAMTOOLS_118=/path/to/samtools-1.18/bin/samtools

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

Outputs in `RQ998.TopNoisyPositions3.BQ30` include position, read, unique-read, protein and unique-protein tables, with read memberships and original multiplicities used downstream.

### Select distinct isoforms

Julia generates repeated heuristic independent-set solutions. The largest saved solution is not guaranteed to be globally maximum or maximal.

The selection seed defaults to 1892. `--firstcolpos 15` is the 1-based index of the first amino-acid position column in the squid unique-protein tables, after 14 metadata columns. Octopus tables have 15 metadata columns and use `--firstcolpos 16`.

`maximal_independent_set_5.jl` automatically appends a timestamp to each generated `DistinctUniqueProteins` CSV filename, so different runs normally produce separate solution files. Expression assignment must select the intended generated solution files; the accepted examples below retain exact timestamped filenames where needed. Use a wildcard only when it matches exactly the intended file.

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

### Assign isoform expression

Assigned expression includes eligible support reassigned from compatible reads, in addition to original direct read support. Expression assignment defaults to fraction 1.0 and the Ascending/Descending algorithms. The entropy-considered run below uses the specified native solution files and reads:

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
> $WD/expressionlevels.regular.out &
```

### Ten-percent comparison for Figure S10

The 0.1-fraction selections and assigned-expression outputs form a separate comparison with the full-data selection.

#### Select distinct isoforms

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)

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

#### Assign isoform expression

```bash
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
> $WD/expressionlevels.regular.fraction01.entropy_considered.out &
```

### Amino-acid dissimilarity comparisons for Figure 3

These comparisons use Miyata amino-acid groups and the Grantham 100 threshold with the grouping, score cutoff and validator below.

#### Miyata amino-acid groups

##### Select distinct isoforms

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

##### Assign isoform expression

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
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.AAgroupsMiyata1979.EntropyConsidered.out &
```

#### Grantham 100 comparison

##### Select distinct isoforms

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

##### Assign isoform expression

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
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.GRANTHAM1974-100.EntropyConsidered.out &
```

## Squid Illumina

**SRR28713507 — PooledShortReads** provides paired-end FASTQ data for the original 19-transcript analysis. Place the input files at:

- `D.pealeii/Data/Illumina/Raw/reads_1.fastq.gz`
- `D.pealeii/Data/Illumina/Raw/reads_2.fastq.gz`

Prepared reads are written to `D.pealeii/Data/Illumina/TrimmedWoDup`.

[DataTable.Squid.ShortReads.csv](Code/Alignment/DataTable.Squid.ShortReads.csv) maps target regions to per-transcript BAMs. It contains 20 targets; `VGA_BPS13` is the target absent from the final 19-gene analysis and Table S2. The manuscript reports successful PCR amplification of 19 of the 20 targets.

### Prepare short reads

Preparation filters by mean quality and trims reads to 265 bases. Preparation-version evidence is recorded in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md).

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

### Align

Use bwa-mem2 and samtools with the Illumina alignment/filtering versions in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md).

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

### Detect RNA editing and infer protein sequences

This step retains positions, reads, unique-read, protein and unique-protein tables, together with the original-read mapping. These tables supply the downstream analyses; the raw samtools `.pileup` intermediates are not required after conversion.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/Illumina

nohup python Code/pileup.py \
--samtools_path "$SAMTOOLS" \
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
>> D.pealeii/MpileupAndTranscripts/Illumina/pileup.out &
```

### Select distinct isoforms

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

### Assign isoform expression

This illustrative command processes three of the 19 Illumina genes—TRIM2, CA2D3 and K0513—in that order, using the reads, unique-protein tables and distinct solutions generated above. It uses one Julia invocation with 60 threads and is not the complete 19-gene manuscript execution.

```bash
WD="D.pealeii/MpileupAndTranscripts/Illumina"

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $WD/comp141044_c0_seq2.DistinctUniqueProteins.*.csv \
$WD/comp140439_c0_seq1.DistinctUniqueProteins.*.csv \
$WD/comp141840_c0_seq2.DistinctUniqueProteins.*.csv \
--allprotsfiles $WD/reads.sorted.aligned.filtered.comp141044_c0_seq2.unique_proteins.csv \
$WD/reads.sorted.aligned.filtered.comp140439_c0_seq1.unique_proteins.csv \
$WD/reads.sorted.aligned.filtered.comp141840_c0_seq2.unique_proteins.csv \
--allreadsfiles $WD/reads.sorted.aligned.filtered.comp141044_c0_seq2.reads.csv \
$WD/reads.sorted.aligned.filtered.comp140439_c0_seq1.reads.csv \
$WD/reads.sorted.aligned.filtered.comp141840_c0_seq2.reads.csv \
--samplenames TRIM2 CA2D3 K0513 \
--firstcolpos 15 \
--fractions 1.0 \
--algs Ascending Descending \
--onlymaxdistinct \
--innerthreadedassignment \
--considerentropy \
--outdir $WD \
--postfix_to_add .EntropyConsidered \
> $WD/expressionlevels.example.out &
```

The output is `<gene>.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv` for each gene, plus selection metadata with the suffix `.Updated.EntropyConsidered.csv`.

### Matched-depth comparison for Figure S6

This branch samples **80,000 read pairs** with `--consider_parity_when_sampling` and seed 1892 for Figure S6c, separately from the full-data Illumina analysis.

#### Detect RNA editing and infer protein sequences

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/Illumina80KPairs

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
--data_table Code/Alignment/DataTable.Squid.ShortReads.csv \
> D.pealeii/MpileupAndTranscripts/Illumina80KPairs/pileup.out &
```

#### Select distinct isoforms

```bash
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/Illumina80KPairs/DistinctProteins.regular.log
```

## Additional squid PacBio biological replicates

| Public accession | Biological replicate | Original local input |
|---|---|---|
| SRR41297225 | JR1 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR1.r841072252411_C01.hifireads.bam` |
| SRR41297224 | JR2 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR2.r841072252411_C01.hifireads.bam` |
| SRR41297223 | JR3 | `D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR3.r841072252411_C01.hifireads.bam` |

The additional squid PacBio biological replicates (JR1–JR3) start from delivered **HiFi/CCS BAMs with UMIs**; provider software is listed in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md).

These replicates target GRIA2 (`comp141693_c0_seq1`), ADAR1 (`comp134400_c0_seq1_extended`) and IQEC1 (`comp141565_c6_seq3`). Biological replicate identity is used throughout UMI processing, per-site comparisons and assigned-expression summaries.

### Align the delivered HiFi reads

Use pbmm2 and samtools with the extended Jan2025 reference and the additional-replicate alignment versions in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md). The `--merge_by_chrom_bams` option merges JR1/JR2/JR3 alignments per transcript.

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
--known_sites_bed_file D.pealeii/Annotations/D.pea.EditingSites.bed \
> D.pealeii/Alignment/AdditionalUMILongReads/align.out &
```

### All reads

#### Detect RNA editing and infer protein sequences

[DataTable.Squid.PB3.AllReads.csv](Code/Alignment/DataTable.Squid.PB3.AllReads.csv) selects merged BAMs under `D.pealeii/Alignment/AdditionalUMILongReads/ByChrom/`, before barcode recognition or deduplication, with `GRIA`, `ADAR1` and `IQEC1` labels, Jan2025 coordinates/strands and no problematic-region mask. For this pileup branch, use the executable documented in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md); set `SAMTOOLS_1221` below and pass it with `--samtools_path`. Recognizable-barcode and deduplicated reads use separate input tables below.

```bash
SAMTOOLS_1221=/path/to/samtools-1.22.1/bin/samtools

mkdir -p D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--samtools_path "$SAMTOOLS_1221" \
--out_dir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads \
--processes 3 \
--threads 15 \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.PB3.AllReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/pileup.out &
```

#### Select distinct isoforms

```bash
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/DistinctProteins.Regular.log
```

#### Assign isoform expression

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
> $UMI_DIR/expressionlevels.EntropyConsidered.out &
```

### UMI-based analyses

#### Prepare barcode-filtered and deduplicated BAMs

[additional_umi_long_reads_2b.ipynb](Code/Notebooks/additional_umi_long_reads_2b.py) uses the delivered HiFi BAMs, per-transcript merged alignments and original all-read tables. Follow its barcode-recognition and replicate-mapping sections, then the two separate export sections:

- **“Save reads with recognized barcodes, before de-duplication”** keeps each gene's `OldRead` names from `best_gene_specific_pcr_amplified_concat_alignments_df` (recognized gene-specific and PCR barcodes). It filters the merged BAMs by these original read names and writes the same BAM basenames under `D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes/`, retaining reads before UMI deduplication.
- **“Save reads after de-duplication”** keeps `OldRead` names from `distinct_umis_gene_specific_pcr_amplified_concat_alignments_df` and writes to `D.pealeii/Alignment/AdditionalUMILongReads.DeduppedReads/`. The preceding deduplication groups indistinguishable UMIs within each gene and biological replicate (`min_shared_umi_sub_seq_len = 5`, `max_errors = 1`), selects a maximum-entropy read per connected component and resolves ties with seed 1892.

No population-switch setting is required: run each export after its own prerequisite calculations. The sections also save `ReadsWithRecognizableBarcodes.tsv` and `DedupedReadsByUMISeq.tsv` in the all-read table directory. Their BAM-indexing commands are commented out; create the `.bai` indexes before the downstream pileup steps. The two input tables linked below select these separate BAM populations. Their `Start` coordinates are zero-based and `End` coordinates are exclusive; empty `ProbRegionsBED` fields specify no problematic-region mask.

#### Reads with recognizable barcodes

##### Detect RNA editing and infer protein sequences

Use the UMI notebook's barcode-filtered BAMs and [DataTable.Squid.PB3.ReadsWithRecognizableBarcodes.csv](Code/Alignment/DataTable.Squid.PB3.ReadsWithRecognizableBarcodes.csv) for the following stages.

```bash
OUT_DIR=D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes

mkdir -p "$OUT_DIR"

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir $OUT_DIR \
--processes 3 \
--threads 20 \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.PB3.ReadsWithRecognizableBarcodes.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> $OUT_DIR/pileup.out &
```

##### Select distinct isoforms

```bash
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.ReadsWithRecognizableBarcodes/DistinctProteins.Regular.log
```

##### Assign isoform expression

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
> $UMI_DIR/expressionlevels.out &
```

#### UMI-deduplicated reads

##### Detect RNA editing and infer protein sequences

Use [DataTable.Squid.PB3.DeduppedReads.csv](Code/Alignment/DataTable.Squid.PB3.DeduppedReads.csv) with the BAMs exported by the UMI notebook's “Save reads after de-duplication” section and their `.bai` indexes.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads \
--processes 3 \
--threads 15 \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.PB3.DeduppedReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads/pileup.out &
```

##### Select distinct isoforms

```bash
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads.DeduppedReads/DistinctProteins.Regular.log
```

##### Assign isoform expression

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
> $UMI_DIR/expressionlevels.out &
```

## Octopus Iso-Seq

**BioProject PRJNA791920** supplies the octopus long-read dataset. This workflow starts from polished CCS files provided by the authors of [*MicroRNAs are deeply linked to the emergence of the complex octopus brain*](https://doi.org/10.1126/sciadv.add9938), under `O.vulgaris/Data/PRJNA791920/IsoSeqPolished/CCS`. Reproduction requires this polished-read handoff.

The seven sequencing runs correspond to the following tissue samples:

| Run | Tissue sample |
|---|---|
| SRR17321895 | Axial nerve cord |
| SRR17321896 | Non-neuronal tissues mix |
| SRR17321897 | Frontal & Vertical lobe |
| SRR17321898 | Pedunculate & Olfactory lobe |
| SRR17321899 | Stellate g. & Visceral g. |
| SRR17321900 | Sucker |
| SRR17321901 | Retina & Optic lobe |

### Prepare polished CCS reads

Use lima and Iso-Seq refine with [ClontechSMARTer-NEBcDNA.primers.fasta](Code/Data/ClontechSMARTer-NEBcDNA.primers.fasta). Follow the Octopus processing versions in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md); newer versions allowed by `pacbiocomb.yml` are not established substitutes for reproducing the original outputs.

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

### Align

Align the refined reads with pbmm2.

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

#### Detect RNA editing and infer protein sequences

The preceding alignment step also generates `AggregatedByChromBySampleSummary.tsv`, summarizing mapped reads and sample counts per transcript. The pooled analysis uses this summary to select eligible transcripts, then applies the SNP/noise criteria below. Total mapped reads determine transcript eligibility; per-position coverage is a separate filter.

```bash
OUTDIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples

mkdir -p $OUTDIR

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
--keep_bam_files \
> $OUTDIR/pileup.out 2>&1 &
```

#### Select distinct isoforms

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
2>&1 | tee $INDIR/DistinctProteins.regular.log
```

#### Assign isoform expression

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
> $DIR/expressionlevels.out &
```

### Pooled samples, minimum 1,000 mapped reads

This branch supports the coverage comparison. It is a separate eligible transcript population.

#### Detect RNA editing and infer protein sequences

```bash
OUTDIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.PooledSamples

mkdir -p $OUTDIR

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
--keep_bam_files \
> $OUTDIR/pileup.out 2>&1 &
```

#### Select distinct isoforms

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
2>&1 | tee $INDIR/DistinctProteins.regular.log
```

## Octopus single-cell

**BioProject PRJNA796958 / GEO GSE193622** supplies the octopus single-cell data. Upstream preprocessing by **David Gorelik** Z"L is documented in [OctopusSCRNAEditingCombinatorics](https://github.com/KaparaNewbie/OctopusSCRNAEditingCombinatorics) and its [pipeline guide](https://github.com/KaparaNewbie/OctopusSCRNAEditingCombinatorics/blob/main/pipeline/README.md). It produces per-transcript BAMs, alignment summaries, neuronal/non-neuronal read populations and read/cell annotations.

Set the input root to your copy of the upstream handoff:

```bash
SC_INPUT_ROOT=/path/to/octopus-single-cell-inputs
```

Under this root, `Whole_collapsed_sample/` contains the original all-cell BAMs, with `Whole_collapsed_sample.regatedByChromBySampleSummary.tsv` alongside it. The `neuronal_reads/` and `non_neuronal_reads/` subdirectories contain their respective BAMs and `neuronal_summary.tsv` / `non_neuronal_summary.tsv`. Read/cell annotations are in `Whole_collapsed_sample_extracted_read_data.withAnnotations.csv`; the cell-class lookup is [NeuronalityOfAnnotation.tsv](Code/Data/NeuronalityOfAnnotation.tsv). Preserve the original all-cell population: the upstream pipeline's `all_reads/` directory contains only reads with neuronal/non-neuronal annotations and is not an interchangeable input.

Reads are pooled within each specified cell population for RNA editing detection, protein-sequence inference and isoform selection. The resulting gene-level diversity measures describe cell-class populations, rather than individual cells.

| Cell population / minimum mapped reads per transcript | Analysis role |
|---|---|
| All cells, ≥50 | Main all-cell analysis; assigned expression supplies the isoform–cell-type association analysis. |
| All cells, ≥1,000 | Higher-coverage isoform-diversity comparison; expression reassignment is not required for this comparison. |
| Neuronal cells, ≥50 | Isoform diversity among neuronal cells. |
| Non-neuronal cells, ≥50 | Isoform diversity among non-neuronal cells. |

### All cells: 50-read branch

#### Detect RNA editing and infer protein sequences

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
--alignments_stats_table "$SC_INPUT_ROOT/Whole_collapsed_sample.regatedByChromBySampleSummary.tsv" \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir "$SC_INPUT_ROOT/Whole_collapsed_sample" \
--interfix_start ".comp" \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/pileup.out &
```

#### Select distinct isoforms

```bash
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
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins.regular.log
```

#### Assign isoform expression

Generate entropy-considered maximum-selection assignments for the cell-type association analysis.

```bash
OUT_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ExpressionLevels.EntropyConsidered.MaxDistinct

mkdir "$OUT_DIR"

python \
Code/Simulations/prepare_fofns_for_expression.py \
--reads_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ReadsFiles \
--proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ProteinsFiles \
--distinct_proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins \
--proteins_postfix ".gz" \
--reads_postfix ".gz" \
--out_dir $OUT_DIR

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
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/expressionlevels.EntropyConsidered.MaxDistinct.out &
```

### All cells: 1,000-read branch

#### Detect RNA editing and infer protein sequences

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
--alignments_stats_table "$SC_INPUT_ROOT/Whole_collapsed_sample.regatedByChromBySampleSummary.tsv" \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir "$SC_INPUT_ROOT/Whole_collapsed_sample" \
--interfix_start ".comp" \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/pileup.out &
```

#### Select distinct isoforms

```bash
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
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000/DistinctProteins.regular.log
```

### Neuronal cells: 50-read branch

#### Detect RNA editing and infer protein sequences

```bash
OUT_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.Neuronal
ALIGNMENT_DIR="$SC_INPUT_ROOT/neuronal_reads/"

mkdir -p $OUT_DIR

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
> $OUT_DIR/pileup.out &
```

#### Select distinct isoforms

```bash
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
2>&1 | tee $OUT_DIR/DistinctProteins.regular.log
```

### Non-neuronal cells: 50-read branch

#### Detect RNA editing and infer protein sequences

```bash
OUT_DIR=O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.NonNeuronal
ALIGNMENT_DIR="$SC_INPUT_ROOT/non_neuronal_reads"

mkdir -p $OUT_DIR

nohup python Code/pileup_with_subparsers.py \
--samtools_path "$SAMTOOLS" \
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
> $OUT_DIR/pileup.out &
```

#### Select distinct isoforms

```bash
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
2>&1 | tee $OUT_DIR/DistinctProteins.regular.log
```

## Neural-expression annotations

**BioProject PRJNA285380** supplies O. bimaculoides RNA-seq data for neural-expression annotation transferred to O. vulgaris through orthology. Place FASTQs under `O.bim/Data/PRJNA285380/Raw`; [Code/NeuralTranscripts/SraRunTable.txt](Code/NeuralTranscripts/SraRunTable.txt) gives the run/tissue mapping for the 18 runs.

### Prepare short reads

```bash
mkdir -p O.bim/Data/PRJNA285380/Raw

nohup \
python Code/prepare_data.py \
--in_dir O.bim/Data/PRJNA285380/Raw/ \
--out_dir O.bim/Data/PRJNA285380 \
illumina \
--min_qual_mean 30 \
> O.bim/Data/PRJNA285380/prepare_data.out &
```

### Quantify expression and transfer neural annotations

Use Salmon and OrthoFinder with the versions documented in the [software requirements](software_versions/README.md#additional-documented-tool-requirements). The pipeline builds the Salmon index, quantifies O. bimaculoides transcripts, assigns neural-expression classes and maps them to O. vulgaris orthologs.

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

[Code/neural_transcripts.py](Code/neural_transcripts.py) writes `O.vulgaris/Annotations/NeuralVsNonNeuralExpression.BySalmonAndOrthoFinder.tsv`, the annotation used by the Figure S14 producer.

## Simulations and other processing

Complete-noise processing uses the samtools/HTSlib versions in the [stage-specific version table](software_versions/historical_stages.2026-11-23.md) and the **Jan2025** FASTA/CDS intervals with the original public-Excel-derived known-site BED. The original PB1/Illumina workflows above use the original reference.

### PB1 complete-noise positions

```bash
SAMTOOLS_118=/path/to/samtools-1.18/bin/samtools
OUTDIR=D.pealeii/MpileupAndTranscripts/CompleteNoisePositions

mkdir -p "$OUTDIR"

nohup \
python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
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
> "$OUTDIR/pileup.OnlyNoisePositions.PB1.out" 2>&1 &
```

### Illumina complete-noise positions

```bash
SAMTOOLS_118=/path/to/samtools-1.18/bin/samtools
OUTDIR=D.pealeii/MpileupAndTranscripts/CompleteNoisePositions

mkdir -p "$OUTDIR"

nohup \
python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
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
--data_table Code/Alignment/DataTable.Squid.ShortReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
--keep_non_refbase_noisy_positions \
--stop_after_making_noise_positions \
> "$OUTDIR/pileup.OnlyNoisePositions.Illumina.out" 2>&1 &
```

These branches retain raw pileups for inspection or reruns and stop after exporting noise-position tables and original-read mappings. The downstream noise analysis uses those tables and mappings.

### Empirical graph assessment — Figure S16

#### Generate simulated proteins

Generate complete and partially observed simulated protein populations with the original input settings and random seeds in [Code/empirical_graph_assessment.py](Code/empirical_graph_assessment.py).

```bash
COMB

mkdir -p Simulations/GraphAssessment

nohup \
python Code/empirical_graph_assessment.py \
> Simulations/GraphAssessment/empirical_graph_assessment.out &
```

#### Select distinct isoforms

Write the consistent-fraction selections to `Simulations/GraphAssessment/SameFrac`:

```bash
mkdir Simulations/GraphAssessment/SameFrac

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
2>&1 | tee Simulations/GraphAssessment/SameFrac/DistinctProteins.regular.log
```

#### Assess false-positive isoforms

The false-positive comparison uses the generated file-of-filenames inputs:

```bash
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
2>&1 | tee Simulations/GraphAssessment/mis_5_assessment.log
```

### Editing simulation — Figure S19

The simulation models repeated RNA editing and read removal to examine how correlations between editing sites depend on editing and removal probabilities. The original analysis is implemented in [unifinished_editing_isoforms.jl](Code/Simulations/unifinished_editing_isoforms.jl).

For Figure S19, [manuscript_s19.jl](Code/Simulations/manuscript_s19.jl) runs the relevant simulation and exports the results as tables for figure generation.

```bash
julia --project=. --threads 1 \
  Code/Simulations/manuscript_s19.jl \
  Code/Notebooks/figure_sources/new_run/FigureS19 \
  Code/Notebooks/manuscript_intermediates/new_run/s19
```

Use new output directories for each run. Figure rendering is described in the [figure reproduction guide](docs/figure_reproduction.md).

### Chimeric-read analysis — Figure S20

#### Analyze chimeric-read origins

This analysis tests whether rare RNA editing patterns are compatible with a prefix from one common RNA pattern and a suffix from another. Such compatibility suggests a potential chimeric origin; it does not directly establish molecular chimerism.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis2

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/chimeric_reads_origin_squid_2.jl \
> D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis2/chimeric_reads_origin_squid_2.out &
```

# Scientific notebooks and workflow

The scientific sequence is **editing positions → reads and unique RNA patterns → predicted proteins and unique proteins → distinct isoform solutions → assigned expression → summaries and exports**. Direct support counts original reads; reassigned support contributes weighted expression mass. Protein sequences are inferred computationally from RNA-derived coding sequences; protein molecules were not directly measured.

The notebooks' **Reproducing manuscript results** tables list setup, prerequisites, branch-specific inputs and exports. Run the indicated sections in order; the notebooks also contain exploratory analyses. Cached exports must match the documented cohort, input version and schema.

Notebook links display the full reconstructed `.ipynb` filename and point to its tracked Jupytext `.py:percent` source; `.ipynb` files and saved outputs are not tracked. Ordinary Python scripts retain their `.py` names. Use the reconstruction command below to open the scientific sources in Jupyter.

| Notebook | Analysis / purpose |
|---|---|
| [squids_pacbio_5.BQ30.ipynb](Code/Notebooks/squids_pacbio_5.BQ30.py) | Original squid PacBio (PB1): expression, recoding, site correlations, All-reads entropy, original-read support and domain analysis. |
| [squids_illumina_5.ipynb](Code/Notebooks/squids_illumina_5.py) | Squid Illumina: isoform expression, recoding, site correlations and read-support analyses. |
| [squids_pacbio_5.additional_UMIs.MergedSamples.ipynb](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.py) | Additional squid PacBio all reads: per-site replicate comparisons, isoform recurrence and assigned expression within each gene and biological replicate. |
| [additional_umi_long_reads_2b.ipynb](Code/Notebooks/additional_umi_long_reads_2b.py) | Recognizes barcodes, maps biological replicates and exports barcode-filtered and UMI-deduplicated BAMs for the downstream processing above. |
| [squids_pacbio_5.additional_UMIs.MergedSamples.ReadsWithRecognizableBardoces.ipynb](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.ReadsWithRecognizableBardoces.py) | Assigned-expression summaries for reads with recognizable barcodes, before UMI deduplication. |
| [squids_pacbio_5.additional_UMIs.MergedSamples.DeduppedReads.ipynb](Code/Notebooks/squids_pacbio_5.additional_UMIs.MergedSamples.DeduppedReads.py) | Assigned-expression summaries for UMI-deduplicated reads. |
| [o.vul_tmr50_vs_tmr1000.pooled.ipynb](Code/Notebooks/o.vul_tmr50_vs_tmr1000.pooled.py) | Pooled Octopus Iso-Seq: isoform diversity, expression and neural-enrichment comparisons at the 50- and 1,000-read thresholds; exports 12-mismatch position tables. |
| [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.ipynb](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py) | Octopus single-cell: isoform distributions and full-cohort isoform–cell-type associations expanding original read support with cell annotations and applying Fisher/BH tests (856 genes). |
| [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.ipynb](Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py) | Neuronal/non-neuronal comparisons of isoform diversity using separate cell-class populations and the `FixedTMR50` export. |
| [noise_positions_and_reads_squid.ipynb](Code/Notebooks/noise_positions_and_reads_squid.py) | Squid noise and 12-mismatch analyses using native positions and original read mappings. |
| [noise_positions_and_reads_octopus.ipynb](Code/Notebooks/noise_positions_and_reads_octopus.py) | Octopus noise and haplotype analyses using native positions and original read mappings. |
| [empirical_graph_assessment.ipynb](Code/Notebooks/empirical_graph_assessment.py) | Summarizes recovered and false-positive isoforms from the complete/partially observed protein simulations. |
| [visualize_chimeric_reads_origin_squid_2.ipynb](Code/Notebooks/visualize_chimeric_reads_origin_squid_2.py) | Analyzes rare-protein chimeras using original RNA patterns and selected assignments through `manuscript_s20.py`. |
| [recoding_sites_in_squid.ipynb](Code/Notebooks/recoding_sites_in_squid.py) | Classifies annotation-based potential recoding per CDS triplet. Its historical extended-ADAR1 input is unverified; see the [reference correction](#extended-squid-reference-jan2025) before new use. |
| [various_plots.BQ30.ipynb](Code/Notebooks/various_plots.BQ30.py) | Produces cross-dataset scientific summaries and table exports, including depth/dispersion, correlations and 12-mismatch counts. |
| [final_manuscript_figures.ipynb](Code/Notebooks/final_manuscript_figures.py) | Renders final publication figures from prepared scientific tables and graphical assets. |

New runs of `maximal_independent_set_5.jl` and `expressionlevels.jl` can yield different heuristic solutions; [selected_solutions.tsv](Code/Notebooks/manuscript_input_data/selected_solutions.tsv) identifies the manuscript's fixed choices and source files, but a fresh clone does not include the required large native inputs. For exact reproduction, `manuscript_inputs.py` verifies and extracts retained assignments without rerunning selection or reassignment, `manuscript_expression_tables.py` prepares plot-ready expression tables, and the final notebook renders previously prepared scientific results, as described in the [reproduction guide](docs/figure_reproduction.md).

## Reconstruct the notebook representations

In a fresh clone, reconstruct the notebooks using the repository's [Jupytext configuration](Code/Notebooks/jupytext.toml):

```bash
cd Code/Notebooks
jupytext --to ipynb \
  additional_umi_long_reads_2b.py \
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

Conversion does not execute cells. Resolve local edits before synchronizing an existing pair, and use execution-table headings to locate cells whose reconstructed IDs differ.

# Figures and supplementary tables

[final_manuscript_figures.ipynb](Code/Notebooks/final_manuscript_figures.py) is the main notebook for final figures. It uses scientific results and tables prepared by the analysis notebooks and supporting scripts, without rerunning all upstream calculations. Some figures also include manual artwork or external graphical assets.

See [Publication figure and supplementary-table reproduction](docs/figure_reproduction.md) for retained-input requirements, table preparation, rendering and supplementary-workbook instructions.

# Methodological limitations and future improvements

Heuristic distinct-isoform selection does not guarantee a globally maximum or maximal solution.

**Aggregation of read-level statistics.** When multiple RNA patterns encode the same protein sequence, some metadata in the unique-protein table may reflect a representative input row rather than an aggregate of all supporting reads (see [Code/Pileup/proteins.py](Code/Pileup/proteins.py)). Where necessary, the scientific analysis notebooks recalculate these statistics from the original supporting reads, as demonstrated by the per-isoform editing summaries in [squids_pacbio_5.BQ30.ipynb](Code/Notebooks/squids_pacbio_5.BQ30.py). A systematic review of the remaining aggregated fields and their downstream uses would be a useful future improvement.

Historical workflows are archived in [Logs.Final.md](Logs.Final.md); environment requirements and validation records are linked from [software_versions/README.md](software_versions/README.md).
