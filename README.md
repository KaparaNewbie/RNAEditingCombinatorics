An accompanying repo to the paper:

<!-- <span style="font-size:2em;">A-to-I editing generates unparalleled complexity in the neural proteome of cephalopods</span> -->

_**A-to-I editing generates unparalleled complexity in the neural proteome of cephalopods**_

_Kobi Shapira, Ruti Balter, Joshua J C Rosenthal, Erez Y. Levanon & Eli Eisenberg_

<br>

<p align="left">
  <img src="ReadmeImages/readme_image_1.jpeg" width="30%">
&nbsp; &nbsp; 
  <img src="ReadmeImages/readme_image_2.jpeg" width="30%">
  &nbsp; &nbsp; 
  <img src="ReadmeImages/readme_image_3.jpeg" width="30%">
</p>

<br>

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<br>

# Getting help

If you need help of any kind (e.g., running one or more of the following scripts, understanding the various file formats, etc.), feel free to open a new issue.

# Setup

All code was run on CentOS 7.
<br>

1. Clone the repository

```bash
PROJECT_DIR=</path/to/some/dir> # choose your own path :)
cd $PROJECT_DIR
git clone https://github.com/KaparaNewbie/RNAEditingCombinatorics
```

2. Create two conda environments; the first is the main one, while the other strictly serves [PacBio packages](https://github.com/PacificBiosciences/pbbioconda) that require python 2.7

```bash
cd RNAEditingCombinatorics

conda env create -f combinatorics.yml
conda activate combinatorics  # test successful creation

conda env create -f pacbiocomb.yml
conda activate pacbiocomb  # test successful creation
```

Creating the first of these may take a few hours if you are using a (relatively) old conda installation, so make sure to update it and make use of the new mamba solver.

3. Easily change into the working path & activate conda on the fly

```bash
echo "alias COMB="cd ${PROJECT_DIR}/RNAEditingCombinatorics\;conda activate combinatorics"" >> ~/.bashrc
echo "alias PBCOMB="cd ${PROJECT_DIR}/RNAEditingCombinatorics\;conda activate pacbiocomb"" >> ~/.bashrc
source ~/.bashrc
```

From now on, we assume you are inside `${PROJECT_DIR}/RNAEditingCombinatorics`, and the `combinatorics` conda enviroment is activated.

4. Julia

- Install the Julia programing language

```bash
COMB
mkdir Julia
cd Julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz
tar zxvf julia-1.7.2-linux-x86_64.tar.gz
```

- Add Julia to `PATH` (see [here](https://julialang.org/downloads/platform/#running_julia) if you need any help)

- Instantiate the Julia environment

```bash
COMB

julia
pkg> activate .
pkg> instantiate
exit()
```

## scsnv

```bash
# mkdir scSNV
# cd scSNV

git clone https://github.com/GWW/scsnv # Mar 25, 2023 release

# WRITE DOCKERFILE

cd scsnv
docker build -t getting-started .

```

## umicollapse

```bash
cd Code/umicollapse

docker build -t umicollapse:latest .


mkdir test
cd test
curl -O -L https://github.com/CGATOxford/UMI-tools/releases/download/1.0.0/example.bam
samtools index example.bam
cd ..

docker run --rm -v $(pwd):/data umicollapse -i /data/example.bam -o /data/dedup_example.bam

./umicollapse bam -i test/example.bam -o test/dedup_example.bam



conda create --name umicollapse umicollapse
conda activate umicollapse

umicollapse bam -i example.bam -o dedup_example.bam
```

# Demo

To make sure thing are set up correctly, you can try running the demo which captures many of the capabilities of the code in the repository.
After completting the setup section, you'll also need to get the latest squid annotation (see below).
Then, follow the instructions at `Demo/DEMO.MD`.

# Data & annotations

## Latest squid annotations

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

Convert the known editing sites excel file into a 6-col bed file:

```python
import pandas as pd
df = pd.read_excel("D.pealeii/Annotations/editing_sites.xlsx", sheet_name="D.pea Editing sites")
# df.to_csv("D.pealeii/Annotations/D.pea.EditingSites.csv", index=False)
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

Also extract the ORFs coordinates into a 6-bed file:

```bash
python Code/orfs_fasta_to_bed.py \
--in_fasta D.pealeii/Annotations/orfs_squ.fa \
--out_bed D.pealeii/Annotations/orfs_squ.bed
```

## 2017 squid annotations

```bash
mkdir D.pealeii/Annotations/December2017
cd D.pealeii/Annotations/December2017

wget https://www.tau.ac.il/~elieis/squid/Archive/v2/orfs_squ.fa
samtools faidx orfs_squ.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx

COMB
```

## Preparing squid long reads data

```bash
mkdir -p D.pealeii/Data/Raw
```

Download the raw demultiplexed data into `D.pealeii/Data/Raw`, and then process it:

```bash
mkdir D.pealeii/Data/CCS

nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Raw \
--out_dir D.pealeii/Data/CCS \
>> D.pealeii/Data/CCS/prepare_data.out &
```

## Preparing squid short reads data

Getting Ruti's reads

```bash
mkdir -p D.pealeii/Data/Illumina/Raw
```

Download the raw data into `D.pealeii/Data/Illumina/Raw` s.t. you have two files: `D.pealeii/Data/Illumina/Raw/reads_1.fastq.gz` and `D.pealeii/Data/Illumina/Raw/reads_1.fastq.gz`, and then process them so:

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

## O.vulgaris annotations

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

### O.vul annotations for single-cell analysis

Creating a GTF from the ORFs bed file to allow for single-cell analysis:

<!-- ```bash
# get the chosen AGAT container version
docker pull quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0



# use an AGAT's tool e.g. agat_convert_sp_gxf2gxf.pl
docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v O.vulgaris/Annotations:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_bed2gff.pl \
--bed /dir/orfs_oct.bed \
-o O.vulgaris/Annotations/orfs_oct.gff



cd O.vulgaris/Annotations

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_bed2gff.pl \
--bed /dir/orfs_oct.bed \
-o /dir/orfs_oct.gff


# docker run \
# -u $(id -u ${USER}):$(id -g ${USER}) \
# -v $PWD:/dir \
# -i \
# "quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0"
``` -->

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

<!--
```python
import pandas as pd
import numpy as np

gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_file = "/private7/projects/Combinatorics/O.vulgaris/Annotations/oct.gff"
gff_df = pd.read_csv(gff_file, sep='\t', comment='#', names=gff_cols)

# gff_df["attributes"] = gff_df["attributes"].apply(lambda x: dict([y.split("=") for y in x.split(";")]))
# gff_df["ID"] = gff_df["attributes"].apply(lambda x: x.get("ID"))
# gff_df["Parent"] = gff_df["attributes"].apply(lambda x: x.get("Parent", np.nan))

mrna_df = gff_df.loc[gff_df["type"] == "mRNA"]
cds_df = gff_df.loc[gff_df["type"] == "CDS"]

mrna_with_no_cds_df = mrna_df.loc[~mrna_df["seqid"].isin(cds_df["seqid"])]
``` -->

<!-- ```bash
cd O.vulgaris/Annotations

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_bed2gff.pl \
--bed /dir/oct.bed \
--primary_tag mRNA \
-o /dir/oct.mRNA.gff

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_bed2gff.pl \
--bed /dir/orfs_oct.bed \
--primary_tag CDS \
-o /dir/orfs_oct.CDS.gff


docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_sp_merge_annotations.pl \
-h

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_sp_merge_annotations.pl \
-gff /dir/oct.mRNA.gff \
-gff /dir/orfs_oct.CDS.gff \
-o /dir/orfs_oct.gff

# verify 1 feature per transcript
tail -n +2 oct.mRNA.gff | cut -f 1 | sort | uniq -c | awk '{$1=$1;print}' | cut -f 1 -d " " | sort | uniq
tail -n +2 orfs_oct.CDS.gff | cut -f 1 | sort | uniq -c | awk '{$1=$1;print}' | cut -f 1 -d " " | sort | uniq

bedtools sort -i orfs_oct.gff > orfs_oct.sorted.gff

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_gff3sort.pl \
-gff /dir/oct.mRNA.gff \
-gff /dir/orfs_oct.CDS.gff \
-o /dir/orfs_oct.gff
 --precise test.gff
``` -->

## O.vulgaris Iso-Seq data

Here we are using a polished version of the CCS reads which were kindly given to us by the authors of the paper [MicroRNAs are deeply linked to the emergence of the complex octopus brain](https://www.science.org/doi/10.1126/sciadv.add9938#T1).

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

## O.bim annotations

```bash
mkdir O.bim/Annotations
cd O.bim/Annotations
wget https://www.tau.ac.il/~elieis/squid/orfs_bim.fa
samtools faidx orfs_bim.fa

COMB # back to the main directory of the project
```

## Preparing O.bim short reads' data

Download the 12 RNA-seq samples from BioProject PRJNA285380 into `O.bim/Data/PRJNA285380/Raw`, and then:

```bash
nohup \
python Code/prepare_data.py \
--in_dir O.bim/Data/PRJNA285380/Raw/ \
--out_dir O.bim/Data/PRJNA285380 \
illumina \
--min_qual_mean 30 \
> O.bim/Data/PRJNA285380/prepare_data.out &
```

## O.bim expression quantification -> O.bim neural transcripts -> orthofinder orthologs -> O.vul neural transcripts

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

## New squid long-reads w/ UMIs

```bash
tmux new -s umi_data

COMB

RAW_DATA_DIR=D.pealeii/Data/RawWithUMIs
mkdir $RAW_DATA_DIR
cd $RAW_DATA_DIR
```

CCS and lima commands executed by Azenta:

```bash
@PG     ID:ccs  PN:ccs  VN:6.3.0        DS:Generate circular consensus sequences (ccs) from subreads.   CL:/opt/pacbio/pa-ccs/current/bin/ccs --streamed --suppress-reports --num-threads 200 --log-level INFO --log-file m64296e_241222_071206.ccs.log --report-json m64296e_241222_071206.ccs_reports.json --report-file m64296e_241222_071206.ccs_reports.txt --metrics-json m64296e_241222_071206.zmw_metrics.json.gz --hifi-summary-json m64296e_241222_071206.hifi_summary.json --stderr-json-log --all --bam m64296e_241222_071206.reads.bam --all-kinetics /data/pa/m64296e_241222_071206.consensusreadset.xml

@PG     ID:lima VN:2.7.1 (commit v2.7.1-1-gf067520)     CL:/gwngsfs/gwngs/tools/conda/envs/smrtlink12.0/share/smrtlink/install/smrtlink-release_12.0.0.177059/bundles/smrttools/install/smrttools-release_12.0.0.177059/private/pacbio/barcoding/binwrap/../../../../private/pacbio/barcoding/bin/lima --ccs /gwngsfs/gwngs/data/projects/30-1097162729/r64296e_20241219_203404D01_tmp//rawdata/r64296e_20241219_203404D01/r64296e_20241219_203404D01.hifireads.bam /gwngsfs/gwngs/tools/pipelines/revio_pacbio_workflow/barcodes/Sequel_RSII_96_barcodes_v3/Sequel_RSII_96_barcodes_v3.fasta /gwngsfs/gwngs/data/p
rojects/30-1097162729/r64296e_20241219_203404D01_tmp//rawdata/r64296e_20241219_203404D01/r64296e_20241219_203404D01.bam
 --split-bam --num-threads 8 --peek-guess --same

     \
--ccs \
hifireads.bam \
Sequel_RSII_96_barcodes_v3.fasta \
r64296e_20241219_203404D01.bam \
--split-bam \
--num-threads 8 \
--peek-guess \ # Try to infer the used barcodes subset.
--same # Only keep same barcodes in a pair in output
```

### Merge per-gene samples as "bulk"

```bash
mkdir D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk

COMB

samtools merge \
--threads 10 \
-o D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk/ADAR1.Merged.r64296e203404D01.hifireads.bam \
D.pealeii/Data/RawWithUMIs/30-1097162729/CCS/*ADAR1*.bam

samtools merge \
--threads 10 \
-o D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk/IQEC.Merged.r64296e203404D01.hifireads.bam \
D.pealeii/Data/RawWithUMIs/30-1097162729/CCS/*IQEC*.bam

PBCOMB

pbindex D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk/ADAR1.Merged.r64296e203404D01.hifireads.bam

pbindex D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk/IQEC.Merged.r64296e203404D01.hifireads.bam
```

# Alignment

### Squid long reads

```bash
mkdir -p D.pealeii/Alignment/BestN1

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir D.pealeii/Data/CCS/BasicCCS \
--out_dir D.pealeii/Alignment/BestN1 \
> D.pealeii/Alignment/BestN1/align.out &
```

### Squid long reads w/ UMIs

Using the new transcriptome with long ADAR1 version.

```bash
mkdir -p D.pealeii/Alignment/UMILongReads

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--in_dir D.pealeii/Data/RawWithUMIs/30-1097162729/CCS \
--out_dir D.pealeii/Alignment/UMILongReads \
--processes 6 \
--threads 6 \
whole_transcriptome_isoseq \
--postfix ".hifireads.bam" \
--known_sites_bed_file D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
> D.pealeii/Alignment/UMILongReads/align.25.3.2025.out &
```

- alu 18
- 16:07

### Squid long reads w/ UMIs - merged samples

Using the new transcriptome with long ADAR1 version.

```bash
mkdir -p D.pealeii/Alignment/UMILongReads.MergedSamples

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--in_dir D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk \
--out_dir D.pealeii/Alignment/UMILongReads.MergedSamples \
--processes 2 \
--threads 15 \
pacbio \
--postfix ".hifireads.bam" \
> D.pealeii/Alignment/UMILongReads.MergedSamples/align.25.3.2025.out &
```

- alu 17
- 16:12

```bash
mkdir -p D.pealeii/Alignment/UMILongReads.MergedSamples

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--in_dir D.pealeii/Data/RawWithUMIs/30-1097162729/CCSAsBulk \
--out_dir D.pealeii/Alignment/UMILongReads.MergedSamples \
--processes 2 \
--threads 15 \
pacbio \
--postfix ".hifireads.bam" \
> D.pealeii/Alignment/UMILongReads.MergedSamples/align.25.3.2025.out &
```

- alu 17
- 16:12

### Squid short reads

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

### Squid short reads - With duplicates

```bash
mkdir -p D.pealeii/Alignment/Illumina.WithDup

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir D.pealeii/Data/RutisReads/Trimmed \
--out_dir D.pealeii/Alignment/Illumina.WithDup \
--threads 20 \
illumina \
--postfix ".good.fastq.gz" \
--require_flags 3 \
--exclude_flags 2304 \
--separate_by_chrom \
> D.pealeii/Alignment/Illumina.WithDup/align.30.5.24.out &

# nohup \
# python Code/align.py \
# --genome D.pealeii/Annotations/orfs_squ.fa \
# --in_dir D.pealeii/Data/Illumina/Trimmed \
# --out_dir D.pealeii/Alignment/Illumina.WithDup \
# --threads 40 \
# illumina \
# --require_flags 3 \
# --exclude_flags 2304 \
# --separate_by_chrom \
# > D.pealeii/Alignment/Illumina.WithDup/align.out &
```

- alu 15

### O. vulgaris

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

# Editing detection & distinct proteins finding

## PacBio

### Pileup

Before running the following script, you should create a csv file like this:

| Sample                 | Gene | Region              | Start | End  | Strand | Path                                                                 | ProbRegionsBED |
| ---------------------- | ---- | ------------------- | ----- | ---- | ------ | -------------------------------------------------------------------- | -------------- |
| GRIA-CNS-RESUB.C0x1291 | GRIA | comp141693_c0_seq1  | 170   | 2999 | +      | D.pealeii/Alignment/BestN1/GRIA-CNS-RESUB.C0x1291.aligned.sorted.bam |                |
| PCLO-CNS-RESUB.C0x1291 | PCLO | comp141882_c0_seq14 | 0     | 6294 | +      | D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam |                |

For your convinience, a file named `Code/Alignment/DataTable.Squid.LongReads.csv` is included for this purpose. Modify its content (e.g., if you use different paths) as needed.

```bash
COMB

mkdir -p D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--processes 2 \
--threads 10 \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.LongReads.csv \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/pileup.out &
```

### Distinct proteins

#### Regular

By regular we mean fractions {0.2, 0.4, ..., 1.0}.

##### Finding isoforms

It may be beneficial for you to run the next step inside a `tmux` session.

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

Due to the random nature of the greedy algorithms we use and the long, resources-heavy computational time, we output each distinct proteins file with a `<time-stamp>` to enable running the simulation multiple times with the same output-dir without worrying that the new results will override the previous ones.

##### Calculating expression levels

```bash
# DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.<time-stamp>.csv \
# D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.<time-stamp>.csv"
# ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
# D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# SAMPLESNAMES="GRIA PCLO"

# nohup \
# julia \
# --project=. \
# --threads 30 \
# Code/Simulations/expressionlevels.jl \
# --distinctfiles $DISTINCTFILES \
# --allprotsfiles $ALLROTSFILES \
# --samplenames $SAMPLESNAMES \
# --outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
# > D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.regular.out &

WD="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30"


nohup \
julia \
--project=. \
--threads 50 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $WD/GRIA-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:29:20.csv $WD/PCLO-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:46:24.csv \
--allprotsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames GRIA PCLO \
--firstcolpos 15 \
--fractions 1.0 \
--outdir $WD \
> $WD/expressionlevels.regular.28.7.2025.out &
```
* alu 18
* 4049123



Calculating exp levels considering entropy of each unchosen prot:

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
* alu 18
* 3644062
  

#### Regular - fraction 0.1 only

##### Finding isoforms

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

##### Calculating expression levels

```bash
# DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.<time-stamp>.csv \
# D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.<time-stamp>.csv"
# ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
# D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# SAMPLESNAMES="GRIA PCLO"

# nohup \
# julia \
# --project=. \
# --threads 40 \
# Code/Simulations/expressionlevels.jl \
# --distinctfiles $DISTINCTFILES \
# --allprotsfiles $ALLROTSFILES \
# --samplenames $SAMPLESNAMES \
# --postfix_to_add .Fraction0_1 \
# --fractions 0.1 \
# --outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
# > D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.regular.fraction01.out &


WD="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30"

nohup \
julia \
--project=. \
--threads 40 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $WD/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.*.csv $WD/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.*.csv \
--allprotsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $WD/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz $WD/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames GRIA PCLO \
--postfix_to_add .Fraction0_1 \
--firstcolpos 15 \
--fractions 0.1 \
--outdir $WD \
> $WD/expressionlevels.regular.fraction01.28.7.2025.out &
```
* alu 17
* 1405645

#### Distinct dissimilar: AA_groups_Miyata1979

##### Finding isoforms

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

##### Calculating expression levels

```bash

DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.<time-stamp>.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.<time-stamp>.csv
"""
ALLROTSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz
"""
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 20 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--aagroups AA_groups_Miyata1979 \
--postfix_to_add .AAgroupsMiyata1979 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.AAgroupsMiyata1979.out &
```

#### Distinct dissimilar: GRANTHAM1974

##### Finding isoforms

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

##### Calculating expression levels

```bash
DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.<time-stamp>.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.<time-stamp>.csv
"""
ALLROTSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz
"""
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-100 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.GRANTHAM1974-100.out &
```

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

## Illumina - 80k sampled reads per transcript

### Pileup

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/Illumina80K

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
--out_dir D.pealeii/MpileupAndTranscripts/Illumina80K \
--processes 10 \
--threads 10 \
--keep_pileup_files \
--gz_compression \
--sample_reads \
--num_sampled_reads 80000 \
--seed 1892 \
directed_sequencing_data \
--data_table D.pealeii/Alignment/Illumina/reads.ByChrom/data_table.csv \
> D.pealeii/MpileupAndTranscripts/Illumina80K/pileup.out &
```

### Distinct proteins

#### Finding isoforms

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina80K/*.unique_proteins.csv.gz)


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
--outdir D.pealeii/MpileupAndTranscripts/Illumina80K \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/Illumina80K/DistinctProteins.regular.log
```

## Illumina - with duplicates

### Pileup

```bash
mkdir D.pealeii/MpileupAndTranscripts/Illumina.WithDup

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
--out_dir D.pealeii/MpileupAndTranscripts/Illumina.WithDup \
--processes 5 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
directed_sequencing_data \
--data_table Code/Alignment/DataTable.Squid.ShortReads.WithDuplicates.csv \
>> D.pealeii/MpileupAndTranscripts/Illumina.WithDup/pileup_with_subparsers.30.5.24.out &
```

- alu 15

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

### total_mapped_reads 50

#### Pileup

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

#### Distinct proteins

##### Finding isoforms

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

##### Expression levels

```bash
DIR=O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3

# python \
# Code/Simulations/prepare_fofns_for_expression.py \
# --proteins_dir $DIR/ProteinsFiles \
# --distinct_proteins_dir $DIR/DistinctProteins \
# --proteins_postfix ".gz" \
# --out_dir $DIR

# mkdir q

# nohup \
# julia \
# --project=. \
# --threads 30 \
# Code/Simulations/expressionlevels.jl \
# --distinctfilesfofn $DIR/DistinctProteinsForExpressionLevels.txt \
# --allprotsfilesfofn $DIR/UniqueProteinsForExpressionLevels.txt \
# --samplenamesfile $DIR/ChromsNamesForExpressionLevels.txt \
# --firstcolpos 16 \
# --fractions 0.2 0.4 0.6 0.8 1.0 \
# --outdir $DIR/ExpressionLevels \
# > $DIR/expressionlevels.8.5.25.out &

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
--threads 30 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $DIR/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $DIR/UniqueProteinsForExpressionLevels.txt \
--allreadsfilesfofn $DIR/ReadsForExpressionLevels.txt \
--samplenamesfile $DIR/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--onlymaxdistinct \
--outdir $OUTDIR \
> $DIR/expressionlevels.30.7.25.out &
```

- alu18
- 1269660

```bash
find O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/DistinctProteins -type f -name '*Updated.csv' | while read -r file; do
    tail -n 16 "$file" | awk '{print $3, $10}'
done | sort | uniq -c | sort -nr
```
  21608 Ascending 
  21597 Descending 
     11 Descending 1

So only in very few selected cases, the descending algorithm missess a distinct protein.

### total_mapped_reads 1000

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3
```

#### Pileup

```bash
nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3/pileup.out &
```

#### Distinct proteins

##### Finding isoforms

```bash
mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 24 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3/DistinctProteins.regular.log
```

## Squid long reads w/ UMIs

Using the new transcriptome with long ADAR1 version.

### Pileup

I'm using `--total_mapped_reads 2000` rather than `--total_mapped_reads 50` only as a crude way to use the real sequenced transcripts.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/UMILongReads


nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/UMILongReads \
--processes 20 \
--threads 5 \
--gz_compression \
--keep_pileup_files \
undirected_sequencing_data \
--alignments_stats_table D.pealeii/Alignment/UMILongReads/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 2000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--disregard_alt_base_freq_1 \
--main_by_chrom_dir D.pealeii/Alignment/UMILongReads/ByChrom \
--interfix_start ".r64296e203404D01" \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
--min_mapped_reads_per_position 0 \
> D.pealeii/MpileupAndTranscripts/UMILongReads/pileup.25.3.25.out &
```

- alu 18
- 14:38

Checking number of unique proteins:

```bash
zcat /private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads/ProteinsFiles/comp134400_c0_seq1_extended.unique_proteins.csv.gz | wc -l
# 137039
zcat /private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads/ProteinsFiles/comp134400_c0_seq1_extended.proteins.csv.gz | wc -l
#  148926
```

### Distinct proteins

```bash
tmux new -s COMB17

COMB

mkdir -p D.pealeii/MpileupAndTranscripts/UMILongReads/DistinctProteins

INFILES=$(echo D.pealeii/MpileupAndTranscripts/UMILongReads/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 60 --proc 10 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/UMILongReads/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/UMILongReads/DistinctProteins/DistinctProteins.Regular.25.3.2025.log
```

#### Expression levels

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/UMILongReads

python \
Code/Simulations/prepare_fofns_for_expression.py \
--proteins_dir $UMI_DIR/ProteinsFiles \
--distinct_proteins_dir $UMI_DIR/DistinctProteins \
--proteins_postfix ".gz" \
--out_dir $UMI_DIR

mkdir $UMI_DIR/ExpressionLevels

nohup \
julia \
--project=. \
--threads 30 \
Code/Simulations/expressionlevels.jl \
--distinctfilesfofn $UMI_DIR/DistinctProteinsForExpressionLevels.txt \
--allprotsfilesfofn $UMI_DIR/UniqueProteinsForExpressionLevels.txt \
--samplenamesfile $UMI_DIR/ChromsNamesForExpressionLevels.txt \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir $UMI_DIR/ExpressionLevels \
> $UMI_DIR/expressionlevels.26.3.25.out &
```

- alu18
- 12:45

## Squid long reads w/ UMIs - merged samples

Using the new transcriptome with long ADAR1 version.

### Pileup

I'm using `--total_mapped_reads 2000` rather than `--total_mapped_reads 50` only as a crude way to use the real sequenced transcripts.

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples

# cp Code/Alignment/DataTable.Squid.LongReads.csv D.pealeii/Alignment/UMILongReads.MergedSamples/DataTable.Squid.MergedUMILongReads.csv

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples \
--processes 2 \
--threads 15 \
--gz_compression \
directed_sequencing_data \
--data_table D.pealeii/Alignment/UMILongReads.MergedSamples/DataTable.Squid.MergedUMILongReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/pileup.25.3.25.out &
```

- alu 17
- 19:26

> (combinatorics) [kobish@alu17 Combinatorics]$ zcat /private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz | wc -l
> 147741
> (combinatorics) [kobish@alu17 Combinatorics]$ zcat /private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz | wc -l
> 153954

### Distinct proteins

#### Regular

##### Distinct isoforms

```bash
tmux new -s COMB18

COMB

# mkdir -p D.pealeii/MpileupAndTranscripts/UMILongReads.TotalCoverage50/DistinctProteins

INFILES=$(echo D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 60 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove ".r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz" \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/DistinctProteins.Regular.25.3.2025.log
```

- alu 13
- 15:07

```bash
tail -n 15 /private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.DistinctUniq
ueProteins.29.01.2025-22:57:33.csv | cut -f 1-5
# ~29K

tail -n 15 /private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.DistinctUniqueProteins.29.01.2025-18:40:32.csv | cut -f 1-5
# ~23K
```

##### Expression levels

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/ADAR1.Merged.DistinctUniqueProteins.*.csv $UMI_DIR/IQEC.Merged.DistinctUniqueProteins.*.csv \
--allprotsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames ADAR1 IQEC1 \
--firstcolpos 15 \
--onlymaxdistinct \
--outdir $UMI_DIR \
--minloglevel -1000 \
> $UMI_DIR/expressionlevels.7.8.25.out &
```
* alu 18
* 3661886


Considering entropy:

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/ADAR1.Merged.DistinctUniqueProteins.26.03.2025-04:39:41.csv $UMI_DIR/IQEC.Merged.DistinctUniqueProteins.26.03.2025-00:27.csv \
--allprotsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames ADAR1 IQEC1 \
--firstcolpos 15 \
--onlymaxdistinct \
--innerthreadedassignment \
--considerentropy \
--outdir $UMI_DIR \
--postfix_to_add .EntropyConsidered \
> $UMI_DIR/expressionlevels.EntropyConsidered.14.8.25.out &
```
* alu 17
* 2508766





## Squid long reads w/ UMIs - unique reads - merged samples

Using the unique reads found by the notebook `umi_long_read.ipynb` by exact UMI seq of 12 bases.

### Pileup

I'm using `--total_mapped_reads 2000` rather than `--total_mapped_reads 50` only as a crude way to use the real sequenced transcripts.

```bash
OUT_DIR=D.pealeii/MpileupAndTranscripts/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples

mkdir -p $OUT_DIR

nohup python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/Jan2025/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/Jan2025/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir $OUT_DIR \
--processes 2 \
--threads 30 \
--gz_compression \
directed_sequencing_data \
--data_table D.pealeii/Alignment/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples/DataTable.Squid.MergedUMILongReads.csv \
--cds_regions D.pealeii/Annotations/Jan2025/orfs_squ.bed \
> $OUT_DIR/pileup.4.7.25.out &
```

- alu 18
- 15:15

### Distinct proteins

#### Regular

##### Distinct isoforms

```bash
tmux new -s COMB18

COMB

# mkdir -p D.pealeii/MpileupAndTranscripts/UMILongReads.TotalCoverage50/DistinctProteins

INFILES=$(echo D.pealeii/MpileupAndTranscripts/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 80 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove ".r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz" \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples/DistinctProteins.Regular.4.7.2025.log
```

- alu 16
- 15:07

##### Expression levels

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples


nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/ADAR1.Merged.DistinctUniqueProteins.*.csv $UMI_DIR/IQEC.Merged.DistinctUniqueProteins.*.csv \
--allprotsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames ADAR1 IQEC1 \
--firstcolpos 15 \
--fractions 1.0 \
--outdir $UMI_DIR \
> $UMI_DIR/expressionlevels.28.7.2025.out &
```
- alu 14
- 287396

Considering entropy:

```bash
UMI_DIR=D.pealeii/MpileupAndTranscripts/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples

nohup \
julia \
--project=. \
--threads 60 \
Code/Simulations/expressionlevels.jl \
--distinctfiles $UMI_DIR/ADAR1.Merged.DistinctUniqueProteins.04.07.2025-15:59:37.csv $UMI_DIR/IQEC.Merged.DistinctUniqueProteins.04.07.2025-15:34:47.csv \
--allprotsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--allreadsfiles $UMI_DIR/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz $UMI_DIR/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz \
--samplenames ADAR1 IQEC1 \
--firstcolpos 15 \
--onlymaxdistinct \
--innerthreadedassignment \
--considerentropy \
--outdir $UMI_DIR \
--postfix_to_add .EntropyConsidered \
> $UMI_DIR/expressionlevels.EntropyConsidered.14.8.25.out &
```
* alu 17
* 1220256


## O.vul single-cell data

### total_mapped_reads 50

#### Pileup

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

- alu 18
- 26.1.25
- 14:03

#### Distinct proteins

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

#### Expression levels

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

- alu18
- 19:51

### total_mapped_reads 1000

#### Pileup

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

- alu 18
- 27.1.25
- 11:09

#### Distinct proteins

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

#### Expression levels

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

- alu18
- 12:45

### total_mapped_reads 50 - neuronal

#### Pileup

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

- alu 18
- 11:45

#### Distinct proteins

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

#### Expression levels

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

- alu18

### total_mapped_reads 1000 - neuronal

#### Pileup

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

- alu 18
- 27.1.25
- 11:09

#### Distinct proteins

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

#### Expression levels

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

- alu18
- 12:45

### total_mapped_reads 50 - non-neuronal

#### Pileup

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

- alu 16
- 11:50

#### Distinct proteins

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

#### Expression levels

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

- alu18
- 19:51

### total_mapped_reads 1000 - non-neuronal

#### Pileup

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

- alu 18
- 27.1.25
- 11:09

#### Distinct proteins

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

#### Expression levels

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

- alu18
- 12:45

# Simulations

## Unfinished editing isoforms

`Code/Simulations/unifinished_editing_isoforms.jl`

## Empirical graph assessment

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

```bash
tmux new -s comb18
COMB

INFILES=$(echo Simulations/GraphAssessment/*.UniqueProteins.tsv.gz)

julia \
--project=. \
--threads 40 --proc 8 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .UniqueProteins.tsv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir Simulations/GraphAssessment \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee Simulations/GraphAssessment/DistinctProteins.regular.12.12.log
```

- alu 18
- 12.12.2024
- 23:16

ls -l Simulations/GraphAssessment/_.UniqueProteins.tsv.gz | wc -l # 180 files
ls -l Simulations/GraphAssessment/_.DistinctUniqueProteins.\* | wc -l # 180 files

```bash
tmux a -t comb18

julia \
--project=. \
--threads 40 \
Code/Simulations/mis_5_assessment.jl \
--complete_infiles Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_09.Rep1.Complete.UniqueProteins.tsv.gz \
--errored_na_files Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_09.Rep1.Errored.PartiallyUnknown.UniqueProteins.tsv.gz \
--out_files Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_09.Rep1.FalsePositives.tsv.gz \
--testfraction 0.1 \
2>&1 | tee Simulations/GraphAssessment/mis_5_assessment.test.23.12.2024.log
```

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

Doing the distinct proteins analysis but with a single sampling per fraction.

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

- alu 13
- 21.1.2024
- 10:42

Not sure I got it right. Trying to print the names of sampled reads from 2/3 samples of corresponding files:

```bash

tmux a -t comb18

mkdir Simulations/GraphAssessment/SameFracTest

x=1
complete_files=$(echo Simulations/GraphAssessment/*Complete.UniqueProteins.tsv.gz | cut -d' ' -f1-$x)
errored_files=$(echo Simulations/GraphAssessment/*Errored.PartiallyUnknown.UniqueProteins.tsv.gz | cut -d' ' -f1-$x)
INFILES="${complete_files} ${errored_files}"
echo $INFILES


julia \
--project=. \
--threads 40 --proc 6 \
Code/Simulations/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .UniqueProteins.tsv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir Simulations/GraphAssessment/SameFracTest \
--fracstep 0.2 \
--fracrepetitions 1 \
--consistentfracsampling \
--algrepetitions 8 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee Simulations/GraphAssessment/SameFracTest/DistinctProteins.regular.20.1.25.log
```

- alu 18
- 20.1.2024
- 11:22

```bash
fractions=("0.2" "0.4" "0.6" "0.8" "1.0")
for fraction in "${fractions[@]}"; do
    f1=/private7/projects/Combinatorics/Simulations/GraphAssessment/SameFracTest/comp140826_c0_seq1.VPS8_HUMAN.UP0_09.Rep10.Errored.PartiallyUnknown.SampledReads.$fraction.csv
    f2=/private7/projects/Combinatorics/Simulations/GraphAssessment/SameFracTest/comp140826_c0_seq1.VPS8_HUMAN.UP0_09.Rep10.Complete.SampledReads.$fraction.csv
    # diff --brief <(tail -n +2 $f1 | sort) <(tail -n +2 $f2 | sort)
    diff --brief $f1 $f2
done

```

# L-GIREMI

https://github.com/gxiaolab/L-GIREMI

```bash
conda env create -f l_giremi.yml
```

## Annotations

We use the 2025 extended squid annotations for L-GIREMI, as we want to be able to compare it to all squid long-read datasets.

```bash
mkdir -p D.pealeii/L-GIREMI/Annotations
cp D.pealeii/Annotations/Jan2025/* D.pealeii/L-GIREMI/Annotations
```

### GTF

Creating a GTF from the ORFs bed file:

<!-- ```bash
cd D.pealeii/L-GIREMI/Annotations
``` -->

```python
from pathlib import Path

import pandas as pd
import Bio
from Bio import SeqIO

in_fasta = Path("D.pealeii/L-GIREMI/Annotations/orfs_squ.fa")
out_bed = Path("D.pealeii/L-GIREMI/Annotations/squ.bed") # a bed file of the complete transcriptome

records = list(SeqIO.parse(in_fasta, "fasta"))

transcripts = []
transcripts_starts = []
transcripts_ends = []
names = []
strands = []
orfs_starts = []
orfs_ends = []

# record = records[0]

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

exit()
```

```bash
cd D.pealeii/L-GIREMI/Annotations

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_bed2gff.pl \
--bed /dir/squ.bed \
--primary_tag mRNA \
-o /dir/squ.gff

docker run \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $PWD:/dir \
"quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0" \
agat_convert_sp_gff2gtf.pl \
--gff /dir/squ.gff \
-o /dir/squ.gtf
```

### Repeats - with EDTA

Test:

```bash
cd /home/alu/kobish/anaconda3/envs/l_giremi/share/EDTA/test/

# perl EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --threads 10

EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --threads 10
```

Actual analysis:

```bash
L_GIREMI

cd D.pealeii/L-GIREMI/Annotations
```

Shorten transcripts IDs - needed to be <= 13 for EDTA:

```python
from pathlib import Path

import pandas as pd
import Bio
from Bio import SeqIO, SeqRecord

in_fasta = Path("orfs_squ.fa")
out_fasta = Path("orfs_squ.short_ids.fa")
long_to_short_ids_file = Path("long_to_short_ids.csv")

records = list(SeqIO.parse(in_fasta, "fasta"))

original_ids = []
shortened_ids = []
new_records = []

# i = 0
# record = records[i]

# description = record.description.split("\t")


for i, record in enumerate(records):
    original_id = record.id
    seq = record.seq
    transcript_end = len(record.seq)
    shortened_id = f"chr_{i}"
    new_record = SeqRecord.SeqRecord(
        seq, id=shortened_id, description=""
    )
    original_ids.append(original_id)
    shortened_ids.append(shortened_id)
    new_records.append(new_record)

with open(out_fasta, "w") as output_file:
    SeqIO.write(new_records, output_file, "fasta-2line")

long_to_short_ids_df = pd.DataFrame(
    {
        "LongID": original_ids,
        "ShortID": shortened_ids
    }
)

long_to_short_ids_df.to_csv(long_to_short_ids_file, index=False)

exit()
```

```bash
nohup EDTA.pl --genome orfs_squ.short_ids.fa --threads 30 > out.19.5.2025 &
```

- alu 18
- 3677597

### Repeats - with EDTA - git cloned

```bash
COMB

cd Code

git clone https://github.com/oushujun/EDTA.git

conda env create -f EDTA/EDTA_2.2.x.yml
```

### Variant calling

Download the single squid DNA sample in the SRA archive:

```bash
mkdir -p D.pealeii/L-GIREMI/DNA/Raw

# prefetch --type sra --check-all yes --max-size 250GB -O D.pealeii/L-GIREMI/DNA/Raw SRX659582 > prefetech_log.19.5.2025.out 2>&1;

prefetch --type sra --check-all yes --max-size 250GB -O D.pealeii/L-GIREMI/DNA/Raw SRR1522984

fasterq-dump D.pealeii/L-GIREMI/DNA/Raw/SRR1522984/SRR1522984.sra --split-files --skip-technical --outdir D.pealeii/L-GIREMI/DNA/Raw
```

Preparing raw data

```bash
nohup \
python Code/prepare_data.py \
--in_dir O.bim/Data/PRJNA285380/Raw/ \
--out_dir O.bim/Data/PRJNA285380 \
illumina \
--min_qual_mean 30 \
> O.bim/Data/PRJNA285380/prepare_data.out &

# Notebooks

Finally, to analyze and visualize the data you'll need the following notebooks:

1. `Code/Notebooks/test_mapped_bam_parsing.ipynb`
2. `Code/Notebooks/squids_pacbio_5.BQ30.ipynb`
3. `Code/Notebooks/squids_illumina_5.ipynb`
4. `Code/Notebooks/o.vul_6_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.ipynb`
5. `Code/Notebooks/various_plots.BQ30.ipynb` (this last notebook depends on the output of notebooks 2-4)

Do note that some of paths appearing in the notebooks might be slightly different compared to those in this readme file.
```
