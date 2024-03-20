[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

TODO add a brief explanation about the project.


<!-- ![](octopus_with_tentacles_made_from_dna_strands_mixing_test_tubes_in_the_lab.png) -->

<br>

<p align="center">
    <img src="**octopus_with_tentacles_made_from_dna_strands_mixing_test_tubes_in_the_lab**.png"  width="300" height="300">
</p>

# Setup

All code was run on CentOS 7.
<br>     

Assume the project resides in `PROJECT_DIR`.

1. Clone the repository
    ```bash
    cd PROJECT_DIR

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

3. Julia
   1. Install the Julia programing language  
        ```bash
        cd ~
        mkdir Julia
        cd Julia
        wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz
        tar zxvf julia-1.7.2-linux-x86_64.tar.gz
        ```

    2. Add Julia to PATH (see [here](https://julialang.org/downloads/platform/#running_julia) if you need any help)

    3.  Instantiate the Julia environment
        ```bash
        cd RNAEditingCombinatorics

        julia
        pkg> activate .
        pkg> instantiate
        exit()
        ```



# Data & annotations

## Getting data & annotations from previous works

Files Yoav copied: 
`/private7/projects/yoav/squid_combinations`

Ruti's project:
`/private7/projects/ruti/new_data/30-461999153`

```
mkdir /private7/projects/Combinatorics

cd /private7/projects/Combinatorics

mkdir -p D.pealeii/Data/Raw

cp \
/private7/projects/ruti/new_data/30-461999153/rawdata/demultiplexed/* \
D.pealeii/Data/Raw

mkdir D.pealeii/Annotations
cd D.pealeii/Annotations

wget https://www.tau.ac.il/~elieis/squid/orfs_squ.fa
samtools faidx orfs_squ.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx
wget https://www.tau.ac.il/~elieis/squid/ST4.txt.gz
wget https://www.tau.ac.il/~elieis/squid/ST5.txt.gz
gunzip ST4.txt.gz
gunzip ST5.txt.gz

cd -  # return to the main dir of the project
```

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

## Getting annotations used by Ruti for the short reads analysis

```bash
COMB

mkdir D.pealeii/Annotations/December2017
cd D.pealeii/Annotations/December2017

wget https://www.tau.ac.il/~elieis/squid/Archive/v2/orfs_squ.fa
samtools faidx orfs_squ.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx
COMB  # return to the December2017 dir of the project
```
The orfs thing still needs work, right now the coordinates doesn't yield an acceptable reading frames.

<!-- Extracting reading frames from the genome:
```python
from Bio import SeqIO
genome = "/private7/projects/Combinatorics/D.pealeii/Annotations/December2017/orfs_squ.fa"
bed = "/private7/projects/Combination/D.pealeii/Annotations/December2017/orfs_squ.bed"
records = SeqIO.parse(Path(genome).absolute(), "fasta")
orfs = []
for record in records:
    description = record.description.split("\t")
    start_index = description.index("OrfStart") + 1
    end_index = description.index("OrfEnd") + 1
    strand_index = description.index("Strand") + 1
    orf_start = str(int(description[start_index]) - 1)
    orf_end = description[end_index]
    strand = description[strand_index]
    orf = [record.id, orf_start, orf_end, ".", ".", strand]
    orfs.append(orf)
    # seq = record.seq[orf_start:orf_end]
with Path(bed).absolute().open("w") as b:
    for orf in orfs:
        b.write("\t".join(orf) + "\n")
```

```python
import pandas as pd
df = pd.read_excel("D.pealeii/Annotations/December2017/editing_sites.xlsx", sheet_name="D.pea Editing sites")
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
``` -->

```bash
python Code/orfs_fasta_to_bed.py \
--in_fasta D.pealeii/Annotations/orfs_squ.fa \
--out_bed D.pealeii/Annotations/orfs_squ.bed
```
* alu 13
* 11.3.2024

## Getting Ruti's matrices

```bash
mkdir D.pealeii/Data/RutisMatrices


cp -r \
/private7/projects/ruti/Squid/MATRIXes_only_paired_full_comp_base_quality_filtered/MATRIXes_from_A_only_paired_full_comp_base_quality_filtered/. \
D.pealeii/Data/RutisMatrices/
```

## Getting Ruti's reads

```bash
mkdir D.pealeii/Data/RutisReads

mkdir D.pealeii/Data/RutisReads/Raw

# cp \
# /private7/projects/ruti/Squid/fastq/left.fq.gz \
# D.pealeii/Data/RutisReads/Raw
# cp \
# /private7/projects/ruti/Squid/fastq/right.fq.gz \
# D.pealeii/Data/RutisReads/Raw
# mv \
# D.pealeii/Data/RutisReads/Raw/left.fq.gz \
# D.pealeii/Data/RutisReads/Raw/reads_1.fq.gz
# mv \
# D.pealeii/Data/RutisReads/Raw/right.fq.gz \
# D.pealeii/Data/RutisReads/Raw/reads_2.fq.gz


cp -T \
/private7/projects/ruti/Squid/fastq/left.fq.gz \
D.pealeii/Data/RutisReads/Raw/reads_1.fastq.gz
cp -T \
/private7/projects/ruti/Squid/fastq/right.fq.gz \
D.pealeii/Data/RutisReads/Raw/reads_2.fastq.gz

```


## Preparing pabbio reads: from raw bam to high quality css reads (including dupliactes removal)

Note: this step is not mandatory. The number of reads in the existing ccs files match the report provided by GENEWIZ:

```bash
$ echo $(zcat /private7/projects/ruti/new_data/30-461999153/CCS/PCLO-CNS-RESUB.ccs.fastq.gz|wc -l)/4|bc
>>> 93450
$ echo $(zcat /private7/projects/ruti/new_data/30-461999153/CCS/GRIA-CNS-RESUB.ccs.fastq.gz|wc -l)/4|bc
>>> 104363
```

### Manual tests


#### _CCS without chunks_

```bash
cd /private7/projects/Combinatorics

conda activate pacbiocomb

# ccs

nohup \
ccs \
--num-threads 30 \
D.pealeii/Data/Raw/GRIA-CNS-RESUB.C0x1291.subreads.bam \
D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.bam \
> D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.out &

# hifi

nohup \
extracthifi \
D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.bam \
D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.bam \
> D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.out &

# index hifi 

pbindex \
D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.bam \

# remove duplicates

nohup \
pbmarkdup \
--rmdup \
--clobber \
--num-threads 30 \
D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.bam \
D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.WoDup.bam \
> D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.WoDup.out &
```

Now, let's check the number of reads in the hifi and WoDup bams.
The first one should be around 104363.

```bash
$ conda activate combinatorics

$ samtools view -c D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.bam
>>> 92098

$ samtools view -c D.pealeii/Data/CCS.Test/GRIA-CNS-RESUB.C0x1291.ccs.WoChunks.hifi.WoDup.bam
>>> 8828
```

The first one is looking good.
On the other hand, it seems that removing duplicates (at least with parameters I used) doesn't fit here. It's probably because all the reads originate from the same area and thus are considered similar to each other.


#### _Optimization of NaN-depleted positions_

```
samtools view /private7/projects/Combinatorics/D.pealeii/Data/CCS/HiFi/GRIA-CNS-RESUB.C0x1291.ccs.bam | head -n 1
```

```
m54278_210226_220806/4194397/ccs        4       *       0       255     *       *       0       0       ATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCCTGATCACAACGATGTGTTGGTCGAAAAACAGTCAATTTCATGTCTACGTCTTGACAATCATCGTACTGGCTTCCATCGTCGGATTACAACCGACTATGTCGGCATCATGGAATATTCAGATTGGTGCCATTTTTGACACACAAACTCTGACGTCAGCCATAAGATACGCCATCAAGACACAAAACAACAACCATGATGATCATCCACAAATGAAGCTACATGAGACTGTCGACTATTTAGACACATTAGACGCCTACAACCTGAGCAATGCAATCTGTGCTCAGTCGGCTCGGGGAATTTTCGCTCTATATATCCTAAACAACCCCGTGTCTCATCACGTTGTCCAGGCCTACAGTAACGCCCTGCACCTACCCGTCTTTGTCCTGGACGTCGCTGCTAGACATCAGCTGCCTAGTTACCCGTTCGAGATAAACATGAGCCCGTCCTTTGTACCTGCCATAACAGAGCTCATTCATTACTTCCGATGGAGAAGAATATATTACATTTTTGATTCAGATGACGGTTTAATTCACTTACAAGATATTTACGATACTTTTAATGATTCCAACAAACCATTGGATGTTGTTGTTCGTCGAGTCGAAGATTTAGACAACTCACACAAAATGCTACGACAACTTGACAGAATCGATCGGCCATTGGACCCACGCAATATTGTTGTGGTTATGTCCAGTACATCTGCTTACCGACGACTATTAAACCAAATTATGGACGTCGGGATGAACAAAGACGGATATCATTATATACTTGGCGGGCTGGGAATCGCCGAATTAGATTTAAGGAACTTCACCTATGGCGGCGTCACCATTACAGGTTTTCAAATCGTGGATAGAAATAACCCGATTATAAAGAAATTTCGGCGGAGGTGGTCTAGCTTGGAGACAATTGTCTGGCCGGGTGCAGGAAAAAAGATTTCATACGATGAAGCAGCTATTTTCGACAGCATCAGAACTTTAGTGAAAGCCTTGGAAAAAATGCACCGTGCTGATCCAAACGTATTTACAGGCGTAATTCGCAGGGGTGGGCTGTGGAATAACAACACCCGCGGTGTGCACTGCACAGCCTACCCTTTGGTGCCGTGGGTGCTGGGGGACAGAATCGTAGATGCCATCAAAAAGGTCAAGTTTAACGGGTTCACGGGCAATGTACAATTTGATAAAAATCTTCAACGCAAAAATTTTAAACTCGATGTTCTTCAATTAAGCTATGGGACACCGTTACAAAAGGTTGGAGAATGGTCCGATGAAAATAAGTTGATTACAGATATGTACCGCCCGCCAGAACACGACAGGTTTCTCCCGGTGGTAAACGACACAAAAATAATGTCTACACACCACCCAGTCAAATGATTACATTTGATGGCTACACAGACGGTTGCCAGACGGAGATTTAGGGATCCAATGAGGAAGAAGAGGCAGCACCCACGACAACAGAAGACAAGGCCTTTTATTCGCAAGAAATGCCTTATTGGAGCTCAGCGATTTAGCTTGAAGCAGGAATGGCTGTTACAAGGTTTTT    ~~~~~|~~~~~x~~r~z~~i~p~~~~{Oc~~~~~~~~[~w~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~y~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  RG:Z:c5f97c6b/1--1      bc:B:S,1,1      bq:i:51 ec:f:35.7105    np:i:34 rq:f:1  sn:B:f,6.70609,12.5494,6.61481,11.4512    zm:i:4194397
```
What do the tags represent?
* `RG:Z:c5f97c6b/1--1` - read group - doesn't matter
* `bc:B:S,1,1` - Barcode sequence identifying the sample - doesn't matter
* `bq:i:51` - Offset to base alignment quality (BAQ)
* `ec:f:35.7105`
* `np:i:34` - NumPasses (1 for subreads, variable for CCS---encodes number of complete passes of the insert)
* `rq:f:1` - Float in [0, 1] encoding expected accuracy
* `sn:B:f,6.70609,12.5494,6.61481,11.4512` - 4 floats for the average signal-to-noise ratio of A, C, G, and T (in that order) over the HQRegion
* `zm:i:4194397` - ZMW hole number

```
samtools view \
D.pealeii/Data/CCS/WoDup/PCLO-CNS-RESUB.C0x1291.ccs.bam \
> D.pealeii/Data/CCS/WoDup/PCLO-CNS-RESUB.C0x1291.ccs.sam

samtools view \
D.pealeii/Data/CCS/WoDup/GRIA-CNS-RESUB.C0x1291.ccs.bam \
> D.pealeii/Data/CCS/WoDup/GRIA-CNS-RESUB.C0x1291.ccs.sam
```

See `Code/Notebooks/test_unmapped_bam_parsing.ipynb` for GENEWIZ-like plots.

Actually, I could probably filter reads by their



### Python script

```
cd /private7/projects/Combinatorics
conda activate combinatorics

mkdir D.pealeii/Data/CCS

nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Raw \
--out_dir D.pealeii/Data/CCS \
>> D.pealeii/Data/CCS/prepare_data.25.1.22.out &

mkdir D.pealeii/Data/CCS.MinPass6

nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Raw \
--out_dir D.pealeii/Data/CCS.MinPass6 \
--min_passes 6 \
>> D.pealeii/Data/CCS.MinPass6/prepare_data.7.2.22.out &



mkdir D.pealeii/Data/CCS.RQ998

nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Raw \
--out_dir D.pealeii/Data/CCS.RQ998 \
--min_rq 0.998 \
>> D.pealeii/Data/CCS.RQ998/prepare_data.8.2.22.out &


mkdir D.pealeii/Data/CCS.MinPass6.RQ998

nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/Raw \
--out_dir D.pealeii/Data/CCS.MinPass6.RQ998 \
--min_passes 6 \
--min_rq 0.998 \
>> D.pealeii/Data/CCS.MinPass6.RQ998/prepare_data.8.2.22.out &


jupyter nbconvert \
--to html \
/private7/projects/Combinatorics/Code/Notebooks/test_mapped_bam_parsing.ipynb \
--TagRemovePreprocessor.remove_input_tags item


jupyter nbconvert \
--to html \
--execute \
/private7/projects/Combinatorics/Code/Notebooks/test_mapped_bam_parsing.ipynb




```

| Sample | Reads in previous CCS files | Reads in new CCS files | Reads in new CCS files, after removing duplicates |
| ------ | --------------------------- | ---------------------- | ------------------------------------------------- |
| PCLO   | 93450                       | 82651                  | 6185                                              |
| GRIA   | 104363                      | 92098                  | 8828                                              |

| Sample | Reads in new CCS files | Reads in new CCS files, `--min_passes 6` |
| ------ | ---------------------- | ---------------------------------------- |
| PCLO   | 82651                  | *                                        |
| GRIA   | 92098                  | *                                        |


## Preparing short reads' data

```bash
nohup \
python Code/prepare_data.py \
--in_dir D.pealeii/Data/RutisReads/Raw/ \
--out_dir D.pealeii/Data/RutisReads \
illumina \
--min_qual_mean 30 \
--trim_to_len 265 \
> D.pealeii/Data/RutisReads/prepare_data.22.5.22.out &
```
* alu 15
* 22.5.22
* 15:13




## O.vulgaris annotations

```bash
mkdir O.vulgaris/Annotations
cd O.vulgaris/Annotations
wget https://springernature.figshare.com/ndownloader/files/13876385  # https://www.nature.com/articles/s41597-019-0017-6#Sec8
tar -xzf 13876385
rm 13876385
rm octopus_abyss41_k41-scaffolds.fa octopus_abyss81_k81-scaffolds.fa scaffolds_81_abyss2_redudance13c.fasta
```

```python
from Code.EditingUtils.gff3 import divide_to_files_by_type
gff3_file = "O.vulgaris/Annotations/gene_models.chromosomer.gff"
out_dir = "O.vulgaris/Annotations"
prefix = "gene_models.chromosomer."
divide_to_files_by_type(gff3_file, out_dir, prefix)
```

Aligning O.vul's gegome against D.pea's transcriptome, using tblastn, to find O.vul's 
editing sites based on the already known D.pea's ones.


```
COMB && cd O.vulgaris/Annotations/tblastn.chromosomer.D.pea

head -n 50107 chromosomer.fa > chromosomer.test.fa
samtools faidx chromosomer.test.fa

makeblastdb \
-in orfs_squ.fa \
-dbtype nucl \
-parse_seqids \
-title orfs_squ.NucDB \
-out orfs_squ.NucDB


tblastx \
-query chromosomer.test.fa \
-db orfs_squ.NucDB \
-out TBlastX.Results.test.tsv \
-outfmt "6 qaccver saccver pident length qframe mismatch gapopen qstart qend sstart send evalue bitscore" \
-num_threads 20 \
-evalue 1e-6
```



<!-- 
```bash
nohup \
python Code/known_sites_alignment.py \
> O.vulgaris/Annotations/known_sites_alignment.out &
```
* alu 17
* 19.01.2023
* 17:31
* 19063 -->


```bash
python Code/orfs_fasta_to_bed.py \
--in_fasta O.vulgaris/Annotations/orfs_oct.fa \
--out_bed O.vulgaris/Annotations/orfs_oct.bed
```
* alu 13
* 7.2.2023


We change direction and use the known vulgaris site in the transcriptome.

```bash
mv O.vulgaris/Annotations O.vulgaris/GenomicAnnotations

mkdir O.vulgaris/Annotations 
cd O.vulgaris/Annotations 
wget https://www.tau.ac.il/~elieis/squid/orfs_oct.fa
samtools faidx orfs_oct.fa

wget https://www.tau.ac.il/~elieis/squid/editing_sites.xlsx

cd -  # return to the main dir of the project

python \
Code/parse_known_sites.py \
--editing_sites_excel_file O.vulgaris/Annotations/editing_sites.xlsx \
--sheet_name "O.vul Editing sites" \
--out_csv_file O.vulgaris/Annotations/O.vul.EditingSites.csv \
--out_bed_file O.vulgaris/Annotations/O.vul.EditingSites.bed
```


## O.vulgaris Iso-Seq data


### Download pre-processed data

https://www.science.org/doi/10.1126/sciadv.add9938#T1

<!-- ```bash
mkdir -p O.vulgaris/Data/PRJNA791920/Raw
accessions=(SRR17321895 SRR17321896 SRR17321897 SRR17321898 SRR17321899 SRR17321900 SRR17321901)
OUTPUT_DIR=O.vulgaris/Data/PRJNA791920/Raw
for acc in ${accessions[@]}; do prefetch.2.10.8 -C yes -p $acc -O $OUTPUT_DIR && fastq-dump.2.10.8 --split-e --skip-technical -O $OUTPUT_DIR $OUTPUT_DIR/$acc.sra && rm $OUTPUT_DIR/$acc.sra; done
gzip O.vulgaris/Data/PRJNA791920/Raw/*
``` -->


```bash
mkdir -p O.vulgaris/Data/PRJNA791920/IsoSeq/Raw

accessions=(SRR17321895 SRR17321896 SRR17321897 SRR17321898 SRR17321899 SRR17321900 SRR17321901)

OUTPUT_DIR=O.vulgaris/Data/PRJNA791920/IsoSeq/Raw

for acc in ${accessions[@]}; do prefetch.2.10.8 -C yes -p $acc -O $OUTPUT_DIR && fastq-dump.2.10.8 --split-e --skip-technical -O $OUTPUT_DIR $OUTPUT_DIR/$acc.sra && rm $OUTPUT_DIR/$acc.sra; done

gzip O.vulgaris/Data/PRJNA791920/IsoSeq/Raw/*
```


```bash
nohup \
python Code/prepare_data.py \
--in_dir O.vulgaris/Data/PRJNA791920/IsoSeq/Raw \
--out_dir O.vulgaris/Data/PRJNA791920/IsoSeq \
--processes 7 \
pacbio_preprocessed_isoseq \
> O.vulgaris/Data/PRJNA791920/IsoSeq/prepare_data.11.12.22.out &
```
* alu 15
* 11.12.1222
* 14:38
* 1399






### Polished CCS data

https://www.science.org/doi/10.1126/sciadv.add9938#T1

This is the polished version of the CCS reads from the paper which I got from the author himself.


```bash
mkdir -p O.vulgaris/Data/PRJNA791920/IsoSeqPolished/CCS

```


```bash
nohup \
python Code/prepare_data.py \
--in_dir O.vulgaris/Data/PRJNA791920/IsoSeqPolished/CCS \
--out_dir O.vulgaris/Data/PRJNA791920/IsoSeqPolished \
--processes 7 \
--threads 6 \
pacbio_polished_ccs_isoseq \
> O.vulgaris/Data/PRJNA791920/IsoSeqPolished/prepare_data.3.4.23.out &
```
* alu 13
* 3.4.23
* 16:45



```
lima \
--num-threads 20 \
--isoseq \
O.vulgaris/Data/PRJNA791920/IsoSeq/Raw/SRR17321896.fastq \
Code/Data/ClontechSMARTer-NEBcDNA.primers.fasta \
O.vulgaris/Data/PRJNA791920/IsoSeq/Raw/SRR17321896.fl.fastq 
```


<!-- ## O.vulgaris FLAM-seq data


```bash
mkdir -p O.vulgaris/Data/PRJNA791920/FLAMSeq/Raw

accessions=(SRR17321871 SRR17321872 SRR17321873 SRR17321874 SRR17321875 SRR17321892 SRR17321893 SRR17321894)

OUTPUT_DIR=O.vulgaris/Data/PRJNA791920/FLAMSeq/Raw

for acc in ${accessions[@]}; do prefetch.2.10.8 -C yes -p $acc -O $OUTPUT_DIR && fastq-dump.2.10.8 --split-e --skip-technical -O $OUTPUT_DIR $OUTPUT_DIR/$acc.sra && rm $OUTPUT_DIR/$acc.sra; done

gzip O.vulgaris/Data/PRJNA791920/FLAMSeq/Raw/*
```

```bash
nohup \
python Code/prepare_data.py \
--in_dir O.vulgaris/Data/PRJNA791920/FLAMSeq/Raw \
--out_dir O.vulgaris/Data/PRJNA791920/FLAMSeq \
--processes 8 \
pacbio_preprocessed_isoseq \
> O.vulgaris/Data/PRJNA791920/FLAMSeq/prepare_data.13.12.22.out &
```
* alu 15
* 13.12.1222
* 13:11
* 47881


```
cd Code
git clone https://github.com/rajewsky-lab/FLAMAnalysis
``` -->


## O.bim annotations

```bash
mkdir O.bim/Annotations
cd O.bim/Annotations
wget https://www.tau.ac.il/~elieis/squid/orfs_bim.fa
samtools faidx orfs_bim.fa

cd - # back to the main directory of the project


```



## Preparing O.bim short reads' data

```bash
nohup \
python Code/prepare_data.py \
--in_dir O.bim/Data/PRJNA285380/Raw/ \
--out_dir O.bim/Data/PRJNA285380 \
illumina \
--min_qual_mean 30 \
> O.bim/Data/PRJNA285380/prepare_data.16.5.23.out &
```
* alu 13
* 15.5.22

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
> O.vulgaris/Annotations/neural_transcripts.18.5.23.out &
```






# Alignment


```

## Python script


```

### --best-n 1


```bash
conda activate combinatorics

cd /private7/projects/Combinatorics

mkdir -p D.pealeii/Alignment/BestN1

nohup \
python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir D.pealeii/Data/CCS/BasicCCS \
--out_dir D.pealeii/Alignment/BestN1 \
> D.pealeii/Alignment/BestN1/align.2.2.22.out &
```



### Aligning short reads

```bash
mkdir -p D.pealeii/Alignment/Illumina/BWAIndex
cd D.pealeii/Alignment/Illumina/BWAIndex

ln -s /private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa /private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/BWAIndex/orfs_squ.fa
samtools faidx /private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/BWAIndex/orfs_squ.fa
bwa-mem2 index /private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/BWAIndex/orfs_squ.fa


gunzip -c \
D.pealeii/Data/RutisReads/TrimmedWoDup/reads_1.fastq.gz \
> D.pealeii/Alignment/Illumina/reads_1.fastq
gunzip -c \
D.pealeii/Data/RutisReads/TrimmedWoDup/reads_2.fastq.gz \
> D.pealeii/Alignment/Illumina/reads_2.fastq

bwa-mem2 mem \
-t 40 \
-M \
D.pealeii/Alignment/Illumina/BWAIndex/orfs_squ.fa \
D.pealeii/Alignment/Illumina/reads_1.fastq \
D.pealeii/Alignment/Illumina/reads_2.fastq \
> D.pealeii/Alignment/Illumina/reads.sam


samtools sort \
-@ 40 \
-o D.pealeii/Alignment/Illumina/reads.sorted.bam \
D.pealeii/Alignment/Illumina/reads.sam

wc -l D.pealeii/Alignment/Illumina/reads_1.fastq
>> 63,574,840 / 4 = 15,893,710

samtools view -c D.pealeii/Alignment/Illumina/reads.sorted.bam
>> 31913225

samtools view -c --require-flags 3 D.pealeii/Alignment/Illumina/reads.sorted.bam
>> 29042109

# 2304 = remove secondary and supplementary (chimeric) alignments
samtools view -c --excl-flags 2304 D.pealeii/Alignment/Illumina/reads.sorted.bam
>> 31,787,420 / 2 = 15,893,710 
# good, that's exactly twice the number of reads in the fastq file

samtools view -c --require-flags 3 --excl-flags 2304 D.pealeii/Alignment/Illumina/reads.sorted.bam
>> 28988674

samtools view \
-@ 40 \
-h \
--require-flags 3 --excl-flags 2304 \
-o D.pealeii/Alignment/Illumina/reads.sorted.aligned.filtered.bam \
D.pealeii/Alignment/Illumina/reads.sorted.bam

samtools index -@ 40 D.pealeii/Alignment/Illumina/reads.sorted.aligned.filtered.bam
```

samtools view --require-flags 3 --excl-flags 2304 D.pealeii/Alignment/Illumina/reads.sorted.bam | cut -f 3 | sort | uniq

comp104747_c0_seq1
comp109402_c0_seq2
comp115544_c0_seq2
comp116833_c0_seq2
comp117823_c0_seq1
comp126362_c0_seq1
comp126798_c4_seq2
comp127721_c0_seq1
comp127748_c1_seq1
comp128294_c1_seq1
comp131231_c1_seq1
comp133111_c1_seq3
comp133135_c0_seq1
comp133607_c0_seq2
comp134045_c0_seq1
comp134322_c0_seq1
comp134452_c1_seq1
comp134572_c0_seq1
comp134586_c0_seq2
comp135071_c0_seq2
comp135521_c1_seq2
comp135778_c6_seq3
comp136010_c0_seq1
comp136058_c0_seq1
comp136085_c0_seq1
comp136709_c1_seq2
comp136810_c2_seq4
comp136967_c2_seq2
comp137317_c2_seq2
comp137622_c1_seq2
comp138035_c0_seq9
comp138474_c10_seq3
comp138788_c10_seq5
comp138881_c1_seq4
comp139266_c0_seq1
comp139467_c7_seq17
comp139758_c0_seq1
comp139867_c0_seq2
comp140029_c0_seq1
comp140190_c1_seq1
comp140363_c1_seq1
comp140439_c0_seq1
comp140666_c0_seq2
comp140712_c0_seq3
comp140897_c1_seq1
comp140910_c2_seq1
comp140987_c3_seq1
comp140989_c0_seq3
comp141044_c0_seq2
comp141053_c1_seq6
comp141158_c1_seq2
comp141378_c0_seq7
comp141455_c0_seq1
comp141495_c4_seq1
comp141517_c0_seq1
comp141532_c3_seq11
comp141565_c6_seq3
comp141574_c0_seq3
comp141634_c1_seq1
comp141640_c0_seq1
comp141684_c0_seq1
comp141689_c2_seq1
comp141840_c0_seq2
comp141880_c1_seq3
comp141881_c0_seq3
comp141882_c0_seq14
comp141914_c0_seq1
comp141919_c0_seq9
comp141964_c0_seq1
comp3559110_c0_seq1
comp3836_c0_seq1
comp49549_c0_seq1
comp49578_c0_seq1
comp49816_c0_seq1
comp88138_c0_seq1
comp90479_c0_seq1


samtools view \
/private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/reads.sorted.aligned.filtered.bam \
comp141882_c0_seq14 | less


```bash
nohup \
python Code/align.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--in_dir D.pealeii/Data/RutisReads/TrimmedWoDup \
--out_dir D.pealeii/Alignment/Illumina \
--threads 40 \
illumina \
--require_flags 3 \
--exclude_flags 2304 \
--separate_by_chrom \
> D.pealeii/Alignment/Illumina/align.23.5.22.out &
```
* alu 15
* 17:05
* 23.5.22


### O. vulgaris

#### Preprocessed isoseq

```bash
mkdir -p O.vulgaris/Alignment/PRJNA791920/IsoSeq

nohup \
python Code/align.py \
--genome O.vulgaris/Annotations/chromosomer.fa \
--in_dir O.vulgaris/Data/PRJNA791920/IsoSeq/Raw \
--out_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq \
--processes 7 \
--threads 10 \
whole_transcriptome_isoseq \
--gff O.vulgaris/Annotations/gene_models.chromosomer.gff \
> O.vulgaris/Alignment/PRJNA791920/IsoSeq/align.30.12.22.out &
```
* alu 13
* 30.12.22
* 14:09
* 44496


Re-aligning using the transcriptomic data.

```bash
mv O.vulgaris/Alignment O.vulgaris/GenomicAlignment

mkdir -p O.vulgaris/Alignment/PRJNA791920/IsoSeq

nohup \
python Code/align.py \
--genome O.vulgaris/Annotations/orfs_oct.fa \
--in_dir O.vulgaris/Data/PRJNA791920/IsoSeq/Raw \
--out_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq \
--processes 7 \
--threads 10 \
pacbio_isoseq_undirected_seq \
--known_sites_bed_file O.vulgaris/Annotations/O.vul.EditingSites.bed \
> O.vulgaris/Alignment/PRJNA791920/IsoSeq/align.7.2.23.out &
```
* alu 13
* 7.2.23
* 13:51
* 44025


#### Polished


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
> O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/align.3.4.23.out &
```
* alu 13
* 3.4.22
* 17:18


# Mpileup & transcripts' creation

Post-processing aligned reads:
* For each file, keep only reads mapped to the relevant region (`comp141693_c0_seq1` for `GRIA`, `comp141882_c0_seq14` for `PCLO`).
* Moreover, keep only reads with sufficient length?
* Chimeric alignments? (Are there any in `pbmm2`?)
* Keep only primary mappings?

BAMs -> edited transcripts:
* [Samtools' mpileup](http://www.htslib.org/doc/samtools-mpileup.html) with `--output-QNAME` such that each each base will be annotated with its read name, which will then allow us to aggregate edited bases per read.





## PacBio (rq 0.998, v. 2)

### Pileup

```bash
conda activate combinatorics
cd /private7/projects/Combinatorics

mkdir -p D.pealeii/MpileupAndTranscripts/RQ998.2

nohup \
python Code/pileup.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/BestN1/data_table_wo_prob_bed.2.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--min_rq 0.998 \
--assurance_factor 1.5 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/pileup.18.5.22.out &
```
* alu 16
* 18.5.22
* 09:12


### Distinct transcripts



```bash
nohup \
julia \
--threads 50 \
Code/UnorderedNaNDepletion/maximal_independent_nan_depletion_4.jl \
--infiles \
D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv \
--samplesnames GRIA.AllRows PCLO.AllRows \
--idcol Transcript \
--editingcol EditedPositions \
--firstcolpos 9 \
--datatype Reads \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_4.Reads.AllRows.out &
```
* alu 13
* 18.5.22
* 18:21



### Distinct proteins

#### Regular

Finding isoforms:

```bash
nohup \
julia \
--threads 50 \
Code/UnorderedNaNDepletion/maximal_independent_nan_depletion_4.jl \
--infiles \
D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
--samplesnames GRIA.AllRows PCLO.AllRows \
--idcol Protein \
--editingcol MinNonSyns \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_4.Proteins.AllRows.out &
```
* alu 16
* 18.5.22
* 19:28

Calculating expression levels:

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueProteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueProteins.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.out &
```
* alu 16
* 18.11.22
* 14:28
* 13596

#### By AA groups (polar/non-polar/positive/negative)

Finding isoforms:

```bash
INFILES=D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv
echo $INFILES

nohup \
julia \
--project=. \
--threads 40 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_add .AAgroups \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--useAAgroups \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.AAgroups.2.out &
```
* alu 13
* 16.11.22
* 15:15
* 1130

Calculating expression levels:

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.AAgroups.16.11.2022-15:40:10.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--useAAgroups \
--postfix_to_add .AAgroups \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.AAgroups.PCLO.out &
```
* alu 13
* 18.11.22
* 14:30
* 45076



#### By BLOSUM62

```bash
INFILES=D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv
echo $INFILES

nohup \
julia \
--project=. \
--threads 40 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_add .Blosum62 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix BLOSUM62 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.Blosum62.out &
```
* alu 13
* 16.11.22
* 15:35
* 8455

Some processes failed.

```
INFILES=D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv
echo $INFILES

nohup \
julia \
--project=. \
--threads 40 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_add .Blosum62 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix BLOSUM62 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.Blosum62.2.out &
```
* alu 13
* 16.11.22
* 17:26
* 27488


Yet again some failed. TODO?


grep "run_fracrepetition failed" D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.Blosum62.2.out | wc -l
> 3



#### AA_groups_Miyata1979

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.2/*.unique_proteins.csv)
echo $INFILES

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .AAgroupsMiyata1979 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--aagroups AA_groups_Miyata1979 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.AAgroupsMiyata1979.out &
```
* alu 13
* 6.12.22
* 20:03
* 12148


Calculating expression levels:

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.AAgroupsMiyata1979.06.12.2022-20:29:59.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.AAgroupsMiyata1979.06.12.2022-20:48:08.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--aagroups AA_groups_Miyata1979 \
--postfix_to_add .AAgroupsMiyata1979 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.AAgroupsMiyata1979.out &
```
* alu 16
* 7.12.22
* 09:30
* 9820



#### By GRANTHAM1974

50-75-100-125-150

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.2/*.unique_proteins.csv)
echo $INFILES
```

> 50

```
nohup \
julia \
--project=. \
--threads 50 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .GRANTHAM1974-50 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 50 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.GRANTHAM1974-50.out &
```
* alu 15
* 6.12.22
* 20:07
* 13835

> 75

```
nohup \
julia \
--project=. \
--threads 50 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .GRANTHAM1974-75 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 75 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.GRANTHAM1974-75.out &
```
* alu 16
* 6.12.22
* 20:09
* 33065

> 100

```
nohup \
julia \
--project=. \
--threads 50 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .GRANTHAM1974-100 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.GRANTHAM1974-100.out &
```
* alu 15
* 6.12.22
* 21:37
* 23392

> 125

```
nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .GRANTHAM1974-125 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 125 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.GRANTHAM1974-125.out &
```
* alu 13
* 6.12.22
* 21:37
* 5445

> 150

```
nohup \
julia \
--project=. \
--threads 50 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .GRANTHAM1974-150 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 150 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/RQ998.2/maximal_independent_nan_depletion_5.Proteins.GRANTHAM1974-150.out &
```
* alu 16
* 6.12.22
* 22:15
* 50969




Calculating expression levels:

> 50

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-50.06.12.2022-20:24:15.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-50.06.12.2022-20:38:37.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"


nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 50 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-50 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.GRANTHAM1974-50.out &
```
* alu 13
* 7.12.22
* 09:37
* 28181

> 75

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-75.06.12.2022-21:09:37.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-75.06.12.2022-22:01:57.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"


nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 75 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-75 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.GRANTHAM1974-75.out &
```
* alu 13
* 7.12.22
* 09:46
* 48923

> 100

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-100.07.12.2022-08:25:55.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-100.07.12.2022-09:37:48.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"


nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-100 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.GRANTHAM1974-100.out &
```
* alu 13
* 7.12.22
* 13:42
* 1596

> 125

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-125.06.12.2022-22:44:40.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-125.07.12.2022-00:06:52.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"


nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 125 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-125 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.GRANTHAM1974-125.out &
```
* alu 15
* 7.12.22
* 13:44
* 24143

> 150

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-150.07.12.2022-00:03:50.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-150.07.12.2022-02:17:58.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
SAMPLESNAMES="GRIA PCLO"


nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 150 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-150 \
> D.pealeii/MpileupAndTranscripts/RQ998.2/expressionlevels.GRANTHAM1974-150.out &
```
* alu 16
* 7.12.22
* 13:46
* 17980

### Deleting this old experiment

```bash
rm -rf D.pealeii/MpileupAndTranscripts/RQ998.2;
```

## PacBio (rq 0.998, 3 top noisy positions)

### Pileup

```bash
conda activate combinatorics
cd /private7/projects/Combinatorics

mkdir -p D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3

nohup \
python Code/pileup.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/BestN1/data_table_wo_prob_bed.2.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--min_rq 0.998 \
--assurance_factor 1.5 \
--min_percent_of_max_coverage 0.1 \
--snp_noise_level 0.1 \
--top_x_noisy_positions 3 \
--parity SE \
--gz_compression \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/pileup.2.3.23.out &
```
* alu 13
* 2.3.23
* 17:57
* 27336


### Distinct proteins

#### Regular

A new tmux session:

```bash
tmux new -s julia13
```

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/DistinctProteins.regular.3.3.23.log
```
* alu 13
* 3.3.23
* 15:20

Calculating expression levels:

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:53:01.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 30 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/expressionlevels.regular.8.5.23.out &
```
* alu 13
* 8.5.22

Distinct proteins, regular, with availble supporting reads:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 8 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .Regular.AvailableReads \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/DistinctProteins.Regular.AvailableReads.22.5.23.log
```

Calculating expression levels for fraction 1, with aggregated supporting reads and proteins:

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.Regular.AvailableReads.22.05.2023-11:08:27.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.Regular.AvailableReads.22.05.2023-11:26:31.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 50 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--postfix_to_add .Regular.SupportingReadsAndProteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/expressionlevels.Regular.SupportingReadsAndProteins.24.5.23.out &
```
* alu 13
* 24.5.22


Regular, fraction 0.1 only:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .Fraction0_1 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--maxfrac 0.1 \
--fracstep 0.1 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/DistinctProteins.regular.Fraction0_1.11.5.23.log
```
* alu 13
* 11.5.23


Calculating expression levels for fraction 0.1 only:


```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:03:49.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:11:00.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--postfix_to_add .Fraction0_1 \
--fractions 0.1 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/expressionlevels.regular.fraction01.14.5.23.out &
```


Calculating expression levels for fraction 0.1 only, with aggregated supporting reads and proteins:


```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:03:49.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:11:00.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 30 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--postfix_to_add .Fraction0_1.SupportingReadsAndProteins \
--fractions 0.1 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/expressionlevels.regular.fraction01.SupportingReadsAndProteins.24.5.23.out &
```
* alu 16
* 24.5.23


#### AA_groups_Miyata1979


A new tmux session:

```bash
tmux new -s julia16
```

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .AAgroupsMiyata1979 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--aagroups AA_groups_Miyata1979 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/DistinctProteins.AAgroupsMiyata1979.3.3.23.log
```
* alu 16
* 3.3.23
* 15:22


Calculating expression levels:

```bash

DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.03.03.2023-15:49:55.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.03.03.2023-16:09:24.csv
"""
ALLROTSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz
"""
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--aagroups AA_groups_Miyata1979 \
--postfix_to_add .AAgroupsMiyata1979 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/expressionlevels.AAgroupsMiyata1979.8.3.23.out &
```
* alu 16
* 8.3.23
* 13:17
* 44868


#### By GRANTHAM1974

```bash
tmux new -s julia15
```

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/*.unique_proteins.csv.gz)
echo $INFILES
```

> 100

```
julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
--postfix_to_add .GRANTHAM1974-100 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/DistinctProteins.GRANTHAM1974-100.8.3.23.log
```
* alu 15
* 8.3.23
* 11:26



Calculating expression levels:

```bash
DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.08.03.2023-12:37:40.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.08.03.2023-13:53:05.csv
"""
ALLROTSFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz
"""
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-100 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/expressionlevels.GRANTHAM1974-100.8.3.23.out &
```
* alu 15
* 8.3.23
* 14:01
* 24812


## PacBio (rq 0.998, 3 top noisy positions, BQ 30)

### Pileup

```bash
conda activate combinatorics
cd /private7/projects/Combinatorics

mkdir -p D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30

python Code/pileup_with_subparsers.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--processes 2 \
--threads 10 \
--keep_pileup_files \
--gz_compression \
directed_sequencing_data \
--data_table D.pealeii/Alignment/BestN1/data_table_wo_prob_bed.2.csv


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
--data_table D.pealeii/Alignment/BestN1/data_table_wo_prob_bed.2.csv \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/pileup.5.2.24.out &
```
* alu 16
* 18:38
* 959729


### Degenerate reads

```bash
nohup \
python Code/Pileup/degenerate_reads.py \
--transcriptome_file D.pealeii/Annotations/orfs_squ.fa \
--out_dir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--reads_first_col_pos 6 \
--postfix_to_remove .csv.gz \
directly_given_variables \
--reads_files D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz \
--chroms comp141693_c0_seq1 comp141882_c0_seq14 \
--starts 170 0 \
--ends 2999 6294 \
--strands "+" "+" \
--names GRIA2 PCLO \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/degenerate_reads.11.3.2024.out &
```
* alu 13
* 11.3.2024

### Distinct proteins

#### Regular

Finding isoforms:

```bash
tmux new -s julia16

INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.regular.6.2.24.log
```
* alu 16
* 2223

Comparing distinct proteins w/o BQ 30:
```bash
# GRIA before
tail -n +2 D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv \
| cut -f 5 | sort -k 1nr | uniq -c | head -n 5
# GRIA after (BQ 30)
tail -n +2 D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:29:20.csv \
| cut -f 5 | sort -k 1nr | uniq -c | head -n 5

# PCLO before
tail -n +2 D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:53:01.csv \
| cut -f 5 | sort -k 1nr | uniq -c | head -n 5
# PCLO after (BQ 30)
tail -n +2 D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:46:24.csv \
| cut -f 5 | sort -k 1nr | uniq -c | head -n 5
```

Calculating epression levels:

```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:29:20.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:46:24.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 30 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.regular.6.2.24.out &
```
* alu 16
* 2773991

#### Regular, fraction 0.1

Regular, fraction 0.1 only:

```bash
tmux new -s julia13

INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.regular.Fraction0_1.6.2.24.log
```
* alu 13


Calculating expression levels for fraction 0.1 only:


```bash
DISTINCTFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.06.02.2024-10:55:32.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.06.02.2024-11:05:02.csv"
ALLROTSFILES="D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
SAMPLESNAMES="GRIA PCLO"

nohup \
julia \
--project=. \
--threads 40 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--postfix_to_add .Fraction0_1 \
--fractions 0.1 \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.regular.fraction01.6.2.24.out &
```
* alu 13
* 11:14
* 101684
  

#### AA_groups_Miyata1979


A new tmux session:

```bash
tmux new -s julia16
```

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)
echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.AAgroupsMiyata1979.6.2.24.log
```
* alu 16
* 12:57


Calculating expression levels:

```bash

DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.06.02.2024-13:43:37.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.06.02.2024-14:14:55.csv
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
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--aagroups AA_groups_Miyata1979 \
--postfix_to_add .AAgroupsMiyata1979 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.AAgroupsMiyata1979.6.2.24.out &
```
* alu 16
* 15:00
* 2815277


#### By GRANTHAM1974

```bash
tmux new -s julia13
```

Finding isoforms:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/*.unique_proteins.csv.gz)
echo $INFILES
```

```bash
julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/DistinctProteins.GRANTHAM1974-100.6.2.24.log
```
* alu 13
* 12:58



Calculating expression levels:

```bash
DISTINCTFILES="""
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.08.03.2023-12:37:40.csv \
D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.08.03.2023-13:53:05.csv
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
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30 \
--substitutionmatrix GRANTHAM1974 \
--similarityscorecutoff 100 \
--similarityvalidator "<" \
--postfix_to_add .GRANTHAM1974-100 \
> D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expressionlevels.GRANTHAM1974-100.6.2.24.out &
```
* alu 13
* 16:13
* 24812



## Illumina (Ruti's Illumina PE READS, take 2 (shortened ids))



#### Complete run

```bash
nohup python Code/pileup.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/Illumina/reads.ByChrom/data_table.csv \
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
>> D.pealeii/MpileupAndTranscripts/Illumina/pileup.12.7.22.out & 
```
 

TODO match project behavior to cloud

```
INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina/*.unique_proteins.csv)
echo $INFILES

JULIA_PROJECT=.

nohup julia \
--threads 50 --proc 3 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
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
> D.pealeii/MpileupAndTranscripts/Illumina80K/pileup.11.4.23.out & 
```
* alu 16
* 11.4.23
* 16:12


### Distinct proteins


A new tmux session:

```bash
tmux new -s julia16
COMB
```

Finding isoforms:

```bash

INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina80K/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 60 --proc 10 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
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
2>&1 | tee D.pealeii/MpileupAndTranscripts/Illumina80K/DistinctProteins.regular.11.4.23.log
```


```bash

# todo run this after the distinct proteins are found
python \
Code/UnorderedNaNDepletion/prepare_fofns_for_expression.py \
--proteins_dir D.pealeii/MpileupAndTranscripts/Illumina80K \
--proteins_prefix "reads.sorted.aligned.filtered." \
--proteins_postfix ".gz"

# DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration/DistinctProteinsForExpressionLevels.txt)
# ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration/UniqueProteinsForExpressionLevels.txt)
# SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration/ChromsNamesForExpressionLevels.txt)

# nohup \
# julia \
# --project=. \
# --threads 20 \
# Code/UnorderedNaNDepletion/expressionlevels.jl \
# --distinctfiles $DISTINCTFILES \
# --allprotsfiles $ALLROTSFILES \
# --samplenames $SAMPLESNAMES \
# --firstcolpos 16 \
# --fractions 0.2 0.4 0.6 0.8 1.0 \
# --outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration \
# > O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/expressionlevels.regular.KnownSitesOnly.30.3.23.out &
```
<!-- * alu 13
* 30.3.22
* 14:41 -->



## O. vulgaris


Calculating expression levels:

```python
from pathlib import Path

proteins_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles")
distinct_proteins_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFiles")

distinct_proteins_files = list(distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv"))
unique_proteins_files = [
    Path(proteins_dir, f"{distinct_proteins_files.name.split('.')[0]}.unique_proteins.csv")
    for distinct_proteins_files in distinct_proteins_files
]
chroms_names = [
    distinct_proteins_files.name.split('.')[0]
    for distinct_proteins_files in distinct_proteins_files
]

distinct_proteins_list = Path(distinct_proteins_dir, "DistinctProteinsForExpressionLevels.txt")
unique_proteins_list = Path(distinct_proteins_dir, "UniqueProteinsForExpressionLevels.txt")
chroms_names_list = Path(distinct_proteins_dir, "ChromsNamesForExpressionLevels.txt")

lists = [distinct_proteins_files, unique_proteins_files, chroms_names]
lists_file_names = [distinct_proteins_list, unique_proteins_list, chroms_names_list]

for _list, list_file_name in zip(lists, lists_file_names):
    with list_file_name.open("w") as list_file_name:
        list_file_name.write(" ".join(map(str, _list)))
```


```bash
DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFiles/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFiles/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFiles/ChromsNamesForExpressionLevels.txt)

nohup \
julia \
--project=. \
--threads 20 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFiles \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/expressionlevels.regular.20.3.23.out &
```
* alu 16
* 20.3.22
* 22:42
* 36821










```bash
mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFilesKnownSites/*.unique_proteins.csv)

echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.regular.KnownSitesOnly.22.3.23.log


```
* alu 16
* 22.3.23



```python
from pathlib import Path

# proteins_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFilesKnownSites")
# distinct_proteins_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites")

proteins_dir = Path("O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFilesKnownSites")
distinct_proteins_dir = Path("O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites")

distinct_proteins_files = list(distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv"))
unique_proteins_files = [
    Path(proteins_dir, f"{distinct_proteins_files.name.split('.')[0]}.unique_proteins.csv")
    for distinct_proteins_files in distinct_proteins_files
]
chroms_names = [
    distinct_proteins_files.name.split('.')[0]
    for distinct_proteins_files in distinct_proteins_files
]

distinct_proteins_list = Path(distinct_proteins_dir, "DistinctProteinsForExpressionLevels.txt")
unique_proteins_list = Path(distinct_proteins_dir, "UniqueProteinsForExpressionLevels.txt")
chroms_names_list = Path(distinct_proteins_dir, "ChromsNamesForExpressionLevels.txt")

lists = [distinct_proteins_files, unique_proteins_files, chroms_names]
lists_file_names = [distinct_proteins_list, unique_proteins_list, chroms_names_list]

for _list, list_file_name in zip(lists, lists_file_names):
    with list_file_name.open("w") as list_file_name:
        list_file_name.write(" ".join(map(str, _list)))
```


```bash
DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites/ChromsNamesForExpressionLevels.txt)

nohup \
julia \
--project=. \
--threads 20 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteinsFilesKnownSites \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/expressionlevels.regular.KnownSitesOnly.22.3.23.out &
```
* alu 16
* 22.3.22
* 17:18














```bash
mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles.KnownSites.NoNoiseFiltration/*.unique_proteins.csv)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 8 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.regular.KnownSites.NoNoiseFiltration.26.3.23.log


```
* alu 16
* 22.3.23


```bash

python \
Code/UnorderedNaNDepletion/prepare_fofns_for_expression.py \
--proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles.KnownSites.NoNoiseFiltration \
--distinct_proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration 

DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration/ChromsNamesForExpressionLevels.txt)

nohup \
julia \
--project=. \
--threads 20 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/DistinctProteins.KnownSites.NoNoiseFiltration \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/expressionlevels.regular.KnownSitesOnly.30.3.23.out &
```
* alu 13
* 30.3.22
* 14:41



## O.vul polished IsoSeq data

### min_samples 7, min_mapped_reads_per_sample 100, total_mapped_reads 0

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered
```

#### Pileup

```bash
# mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered

cp \
O.vulgaris/Data/PRJNA791920/IsoSeq/Raw/samples.csv \
O.vulgaris/Data/PRJNA791920/IsoSeqPolished

nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--min_samples 7 \
--min_mapped_reads_per_sample 100 \
--min_known_sites 0 \
--pooled_transcript_noise_threshold 0.06 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom/ \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/pileup.3.4.23.out & 
```
* alu 13
* 3.4.23
* 18:10


#### Distinct proteins

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 30 --proc 10 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins.regular.26.3.23.log


```
* alu 16
* 22.3.23


```bash

python \
Code/UnorderedNaNDepletion/prepare_fofns_for_expression.py \
--proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/ProteinsFiles \
--proteins_postfix .gz \
--distinct_proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins

DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins/ChromsNamesForExpressionLevels.txt)

nohup \
julia \
--project=. \
--threads 20 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/expressionlevels.regular.27.4.23.out &
```
* alu 13
* 27.4.23
* 14:06


#### Deleting this old experiment

```bash
rm -rf O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered;
```

### min_samples 0, min_mapped_reads_per_sample 0, total_mapped_reads 100

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100
```

#### Pileup

```bash
nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100 \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--min_samples 0 \
--min_mapped_reads_per_sample 0 \
--total_mapped_reads 100 \
--min_known_sites 0 \
--pooled_transcript_noise_threshold 0.06 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom/ \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100/pileup.13.6.23.out &
```
* alu 11


#### Distinct proteins

```bash
# tmux new -s julia13

# COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 20 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100/DistinctProteins.regular.13.6.23.log
```
* alu 13

#### Deleting this old experiment

```bash
rm -rf O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage100;
```

### min_samples 0, min_mapped_reads_per_sample 0, total_mapped_reads 200

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200
```

#### Pileup

```bash
nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200 \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--min_samples 0 \
--min_mapped_reads_per_sample 0 \
--total_mapped_reads 200 \
--min_known_sites 0 \
--pooled_transcript_noise_threshold 0.06 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom/ \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200/pileup.13.6.23.out &
```
* alu 11

```bash
# tmux new -s julia13

# COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 20 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200/DistinctProteins.regular.14.6.23.log
```
* alu 13


#### Deleting this old experiment

```bash
rm -rf O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage200;
```

### min_samples 0, min_mapped_reads_per_sample 0, total_mapped_reads 50

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50
```

#### Pileup

```bash
nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50 \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--min_samples 0 \
--min_mapped_reads_per_sample 0 \
--total_mapped_reads 50 \
--min_known_sites 0 \
--pooled_transcript_noise_threshold 0.06 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/pileup.13.6.23.out & 
```
* alu 11


#### Distinct proteins

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 20 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins.regular.13.6.23.log
```
* alu 13




### min_samples 0, min_mapped_reads_per_sample 0, total_mapped_reads 1000

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000
```

#### Pileup

```bash
nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000 \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--min_samples 0 \
--min_mapped_reads_per_sample 0 \
--total_mapped_reads 1000 \
--min_known_sites 0 \
--pooled_transcript_noise_threshold 0.06 \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom/ \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000/pileup.23.6.23.out &
```
* alu 13

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000/ProteinsFiles/*.unique_proteins.csv.gz)

julia \
--project=. \
--threads 20 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000/DistinctProteins.regular.26.6.23.log
```
* alu 13


```bash
python \
Code/UnorderedNaNDepletion/prepare_fofns_for_expression.py \
--proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/ProteinsFiles \
--proteins_postfix .gz \
--distinct_proteins_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins

DISTINCTFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins/DistinctProteinsForExpressionLevels.txt)
ALLROTSFILES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins/UniqueProteinsForExpressionLevels.txt)
SAMPLESNAMES=$(cat O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/DistinctProteins/ChromsNamesForExpressionLevels.txt)

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/ExpressionLevels

nohup \
julia \
--project=. \
--threads 20 \
Code/UnorderedNaNDepletion/expressionlevels.jl \
--distinctfiles $DISTINCTFILES \
--allprotsfiles $ALLROTSFILES \
--samplenames $SAMPLESNAMES \
--firstcolpos 16 \
--fractions 0.2 0.4 0.6 0.8 1.0 \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/ExpressionLevels \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/expressionlevels.regular.26.9.23.out &
```
* alu 13


### total_mapped_reads 50, BQ 30, BH only, AHL

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH only" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--keep_bam_files \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/pileup.25.2.24.out &
```
* alu 16
* 22:32
* 3620835




```
/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/PositionsFiles/comp181924_c0_seq4.positions.csv.gz


import pandas as pd
df1 = pd.DataFrame(
    {
        "Position": [1, 2, 3, 4, 5],
    }
)
df2 = pd.DataFrame(
    {
        "Position": [1, 5],
        "PV": [0.1, 0.2],
    }
)
df1 = df1.merge(df2, on="Position", how="left")

covered_coding_non_adenosines_positions_df.loc[covered_coding_non_adenosines_positions_df["TotalCoverage"]!=covered_coding_non_adenosines_positions_df["RefBaseCount"]]
covered_coding_non_adenosines_positions_df.loc[covered_coding_non_adenosines_positions_df["AltBaseCount"] > 0]


cd /private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/PositionsFiles/
zcat *.CorrectedNoise.csv.gz | cat | cut -f 5 | uniq -c

from pathlib import Path
import pandas as pd
positions_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/PositionsFiles")

corrected_noise_files = list(positions_dir.glob("*.CorrectedNoise.csv.gz"))
corrected_noise_dfs = [pd.read_table(corrected_noise_file) for corrected_noise_file in corrected_noise_files]
merged_corrected_noise_df = pd.concat(corrected_noise_dfs)
merged_corrected_noise_df["NoisyCorrected"].value_counts()

positions_files = list(positions_dir.glob("*.positions.csv.gz"))
merged_positions_df = pd.concat([pd.read_table(positions_file) for positions_file in positions_files])
merged_positions_df["NoisyCorrected"].value_counts()

corrected_editing_files = list(positions_dir.glob("*.CorrectedEditing.csv.gz"))
merged_corrected_editing_df = pd.concat([pd.read_table(corrected_editing_file) for corrected_editing_file in corrected_editing_files])
merged_corrected_editing_df
merged_corrected_editing_df["EditedCorrected"].value_counts()

```




#### Distinct proteins

##### Original

```bash
tmux new -s julia16

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 20 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins.regular.26.2.24.log
```
* alu 16


```
F=/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins/comp177903_c0_seq3.DistinctUniqueProteins.26.02.2024-09:38:52.csv
F=/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins/comp178324_c0_seq2.DistinctUniqueProteins.26.02.2024-13:15:15.csv


cut -f 5 $F | sort -k 1nr | uniq -c | head -n 5

tail -n 16 $F | cut -f 1-5 | sort -k 1nr | uniq -c | head -n 5
tail -n +16 $F | cut -f 5 | sort -k 1nr | uniq -c | head -n 5
```

##### Fixed reading "Reads" col from files containing a sinle unique protein

```bash
tmux new -s julia1

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins2

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 28 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH/DistinctProteins.regular2.28.2.24.log
```
* alu 13
* 13:20


### total_mapped_reads 1000, BQ 30, BH only, AHL

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH \
--processes 12 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH only" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/pileup.26.2.24.out &
```
* alu 16
* 18:30
* 169548



#### Distinct proteins

##### Original


```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/DistinctProteins.regular.26.2.24.log
```
* alu 13
* 20:57

##### Fixed reading "Reads" col from files containing a sinle unique protein

```bash
tmux new -s julia16


COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/DistinctProteins2

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 28 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/DistinctProteins2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH/DistinctProteins.regular2.28.2.24.log
```
* alu 16
* 12:15

### total_mapped_reads 50, BQ 30, BH after noise, AHL

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise \
--processes 10 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--keep_bam_files \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/pileup.26.2.24.out & 
```
* alu 13
* 11:21
* 3520267


#### Distinct proteins

##### Original

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/DistinctProteins.regular.26.2.24.log
```
* alu 13

##### Fixed reading "Reads" col from files containing a sinle unique protein

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/DistinctProteins2

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 28 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/DistinctProteins2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/DistinctProteins.regular2.28.2.24.log
```
* alu 13
* 12:03

### total_mapped_reads 1000, BQ 30, BH after noise, AHL

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise \
--processes 16 \
--threads 5 \
--keep_pileup_files \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--keep_bam_files \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/pileup.27.2.24.out &
```
* alu 16
* 22:17
* 419866


#### Distinct proteins


##### Original

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/DistinctProteins.regular.28.2.24.log
```
* alu 13
* 11:15

##### Fixed reading "Reads" col from files containing a sinle unique protein

```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/DistinctProteins2


INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/DistinctProteins2 \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/DistinctProteins.regular2.28.2.24.log
```
* alu 13
* 11:15


### total_mapped_reads 50, BQ 30, BH only, AHL, 2

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH only" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
>> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2/pileup.19.3.24.out &
```
* alu 17
* 20:31
* 946222





#### Distinct proteins


```bash
tmux new -s julia17

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 60 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2/DistinctProteins.regular.19.3.24.log
```
* alu 13
* 14:45

#### Pileup test

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test


nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test \
--processes 12 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH only" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
>> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test/pileup.19.3.24.out &
```
* alu 17
* 20:02

#### Distinct proteins test


```bash
tmux new -s julia17

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 40 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BH.2.Test/DistinctProteins.regular.19.3.24.log
```
* alu 17
* 20:07


### total_mapped_reads 1000, BQ 30, BH only, AHL, 2

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH only" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
>> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2/pileup.19.3.24.out &
```
* alu 13
* 20:33



#### Distinct proteins


```bash
tmux new -s julia13

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 40 --proc 8 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BH.2/DistinctProteins.regular.19.3.24.log
```
* alu 17
* 20:57




### total_mapped_reads 50, BQ 30, BH after noise, AHL, 2

#### Pileup

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2
```

```bash
nohup python Code/pileup_with_subparsers.py \
--transcriptome O.vulgaris/Annotations/orfs_oct.fa \
--known_editing_sites O.vulgaris/Annotations/O.vul.EditingSites.bed \
--exclude_flags 2304 \
--parity SE \
--min_rq 0.998 \
--min_bq 30 \
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2 \
--processes 20 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 50 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2/pileup.20.3.24.out & 
```
* alu 17
* 10:30
* 3680182


#### Distinct proteins

```bash
tmux new -s julia16

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.2/DistinctProteins.regular.13.3.24.log
```
* alu 16
* 12:17




### total_mapped_reads 1000, BQ 30, BH after noise, AHL, 2

```bash
mkdir -p O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2
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
--out_dir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2 \
--processes 14 \
--threads 5 \
--gz_compression \
undirected_sequencing_data \
--alignments_stats_table O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv \
--total_mapped_reads 1000 \
--alternative_hypothesis larger \
--final_editing_scheme "BH after noise thresholding" \
--main_by_chrom_dir O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom \
--cds_regions O.vulgaris/Annotations/orfs_oct.bed \
--samples_table O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv \
--min_mapped_reads_per_position 0 \
> O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2/pileup.20.3.24.out &
```
* alu 13
* 10:32
* 4134058


#### Distinct proteins


```bash
tmux new -s julia16

COMB

mkdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2/DistinctProteins

INFILES=$(echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2/ProteinsFiles/*.unique_proteins.csv.gz)

echo $INFILES

julia \
--project=. \
--threads 24 --proc 6 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .unique_proteins.csv.gz \
--idcol Protein \
--firstcolpos 16 \
--datatype Proteins \
--outdir O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2/DistinctProteins \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
2>&1 | tee O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.2/DistinctProteins.regular.13.3.24.log
```
* alu 16
* 11:49