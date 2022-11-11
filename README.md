TODO add a brief explanation about the project.



# Setup

All code was run on CentOS 7.
<br>     

Assume the project resides in PROJECT_DIR.

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




# Mpileup & transcripts' creation

Post-processing aligned reads:
* For each file, keep only reads mapped to the relevant region (`comp141693_c0_seq1` for `GRIA`, `comp141882_c0_seq14` for `PCLO`).
* Moreover, keep only reads with sufficient length?
* Chimeric alignments? (Are there any in `pbmm2`?)
* Keep only primary mappings?

BAMs -> edited transcripts:
* [Samtools' mpileup](http://www.htslib.org/doc/samtools-mpileup.html) with `--output-QNAME` such that each each base will be annotated with its read name, which will then allow us to aggregate edited bases per read.


## Manual tests

### Naive

```
conda activate combinatorics
cd /private7/projects/Combinatorics

samtools view D.pealeii/Alignment/Native/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam comp141882_c0_seq14 | less
```

```
touch D.pealeii/Alignment.Test/two_reads.txt
echo -e "m54278_210226_220806/60358899/ccs" >> D.pealeii/Alignment.Test/two_reads.txt
echo -e "m54278_210226_220806/5374853/ccs" >> D.pealeii/Alignment.Test/two_reads.txt

samtools view \
-N D.pealeii/Alignment/Native/two_reads.txt \
D.pealeii/Alignment/Native/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam comp141882_c0_seq14 \
| less







samtools mpileup \
--positions D.pealeii/Annotations/D.pea.EditingSites.bed \
--region comp141882_c0_seq14 \
--fasta-ref D.pealeii/Annotations/orfs_squ.fa \
--no-BAQ \
--no-output-ins --no-output-ins \
--no-output-del --no-output-del \
--no-output-ends \
--output-QNAME \
D.pealeii/Alignment/Native/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
| less




samtools mpileup \
--positions D.pealeii/Annotations/D.pea.EditingSites.bed \
--region comp141882_c0_seq14 \
--fasta-ref D.pealeii/Annotations/orfs_squ.fa \
--no-BAQ \
--no-output-ins --no-output-ins \
--no-output-del --no-output-del \
--no-output-ends \
--output-QNAME \
D.pealeii/Alignment/Native/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
--output D.pealeii/Alignment.Test/PCLO.mpileup




samtools mpileup \
--positions D.pealeii/Annotations/D.pea.EditingSites.bed \
--region comp141882_c0_seq14 \
--fasta-ref D.pealeii/Annotations/orfs_squ.fa \
--no-BAQ \
--no-output-ins --no-output-ins \
--no-output-del --no-output-del \
--no-output-ends \
--output-QNAME \
--excl-flags 2304 \
D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
--output D.pealeii/Alignment.Test/PCLO.primary.linear.mpileup



samtools view -F 256 D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
| grep "m54278_210226_220806/19268008/ccs" \
| less




samtools view D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
| wc -l
>>> 84980

samtools view -F 2304 D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
| wc -l
>>> 82595

samtools view \
-F 2304 \
-c \
D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
comp141882_c0_seq14
>>> 82575
# it means most of the aligned reads fall within the required region
```


comp141882_c0_seq14


### --best-n 1

```bash
conda activate combinatorics
cd /private7/projects/Combinatorics

samtools view \
-F 2304 \
-c \
D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
comp141882_c0_seq14
>>> 82548

mkdir D.pealeii/Alignment.Test.N1

samtools mpileup \
--positions D.pealeii/Annotations/D.pea.EditingSites.bed \
--region comp141882_c0_seq14 \
--fasta-ref D.pealeii/Annotations/orfs_squ.fa \
--no-BAQ \
--no-output-ins --no-output-ins \
--no-output-del --no-output-del \
--no-output-ends \
--output-QNAME \
--excl-flags 2304 \
D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
--output D.pealeii/Alignment.Test.N1/PCLO.primary.linear.mpileup

cut -f 7 D.pealeii/Alignment.Test.N1/PCLO.primary.linear.mpileup | \
tr ',' '\n' | \
sort | \
uniq | \
wc -l
>>> 8758
```

So the change from 82548 aligned reads in the bam file to just 8758 happens in the creation of the pileup file, so it is not a matter of wrong parsing of the pileup file itself.


```bash
samtools mpileup \
--positions D.pealeii/Annotations/D.pea.EditingSites.bed \
--region comp141882_c0_seq14 \
--fasta-ref D.pealeii/Annotations/orfs_squ.fa \
--no-BAQ \
--no-output-ins --no-output-ins \
--no-output-del --no-output-del \
--no-output-ends \
--output-QNAME \
--excl-flags 2304 \
--max-depth 90000 \
D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
--output D.pealeii/Alignment.Test.N1/PCLO.primary.linear.max_depth_90k.mpileup

cut -f 7 D.pealeii/Alignment.Test.N1/PCLO.primary.linear.max_depth_90k.mpileup | \
tr ',' '\n' | \
sort | \
uniq | \
wc -l
>>> 82548
```

So that problem is fixed. The default `--max-depth 8000` didn't suite my deep sequencing within specific regions.


## Python script


### W/O max depth per position

```bash
conda activate combinatorics

cd /private7/projects/Combinatorics

mkdir -p D.pealeii/MpileupAndTranscripts/UnmappedNaN.BestN1.WoMaxDepth

nohup \
python Code/create_transcripts.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/BestN1/data_table_wo_prob_bed.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/UnmappedNaN.BestN1.WoMaxDepth \
> D.pealeii/MpileupAndTranscripts/UnmappedNaN.BestN1.WoMaxDepth/create_transcripts.3.3.22.out &
```

### W/O max depth per position, use only reads with rq 0.998

```bash
conda activate combinatorics
cd /private7/projects/Combinatorics

mkdir -p D.pealeii/MpileupAndTranscripts/UnmappedNaN.BestN1.WoMaxDepth.RQ998

nohup \
python Code/create_transcripts.py \
--genome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/BestN1/data_table_wo_prob_bed.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/UnmappedNaN.BestN1.WoMaxDepth.RQ998 \
--min_rq 0.998 \
> D.pealeii/MpileupAndTranscripts/UnmappedNaN.BestN1.WoMaxDepth.RQ998/create_transcripts.3.3.22.out &
```


### rq 0.998, v. 2

#### Pileup

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


#### Reads



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



#### Proteins


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



### Ruti's Illumina PE READS

#### Pileup - example

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/IlluminaExample

nohup \
python Code/pileup.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/Illumina/reads.ByChrom/data_table.example.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--top_x_noisy_positions 5 \
--assurance_factor 1.5 \
--include_flags 3 \
--exclude_flags 2304 \
--min_bq 30 \
--parity PE \
--processes 4 \
--threads 20 \
--keep_pileup_files \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/pileup.example.1.6.22.out & 
```

#### Reads & proteins - example


```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/IlluminaExample/*.unique_reads.csv)
echo $INFILES

nohup \
julia \
--threads 20 \
Code/UnorderedNaNDepletion/maximal_independent_nan_depletion_4.jl \
--infiles $INFILES \
--postfix .unique_reads.csv \
--idcol Transcript \
--editingcol EditedPositions \
--firstcolpos 9 \
--datatype Reads \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
# --fracstep 1.0 \
--fracrepetitions 2 \
--algrepetitions 2 \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_nan_depletion_4.Reads.Test.out &

nohup \
julia \
--threads 20 \
Code/UnorderedNaNDepletion/maximal_independent_nan_depletion_4.jl \
--infiles $INFILES \
--postfix .unique_reads.csv \
--idcol Transcript \
--editingcol EditedPositions \
--firstcolpos 9 \
--datatype Reads \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 0.5 \
--fracrepetitions 1 \
--algrepetitions 1 \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_nan_depletion_4.Reads.Test.out &
```
* alu 13
* 2.6.22
* 17:05


```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/IlluminaExample/*.unique_proteins.csv)

JULIA_PROJECT=.

nohup \
julia \
--threads 20 \
Code/UnorderedNaNDepletion/maximal_independent_nan_depletion_4.jl \
--infiles $INFILES \
--postfix .unique_proteins.csv \
--idcol Protein \
--editingcol MinNonSyns \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 1.0 \
--fracrepetitions 1 \
--algrepetitions 1 \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_nan_depletion_4.Proteins.Test.out &
```
* alu 13
* 2.6.22
* 18:51


More proteins tests for fraction 1.0 (the complete data) W/O parallelization:

```bash
INFILES=$(echo D.pealeii/MpileupAndTranscripts/IlluminaExample/*.unique_proteins.csv)
echo $INFILES

JULIA_PROJECT=.

# small test

nohup julia \
--threads 10 --proc 2 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix reads.sorted.aligned.filtered. \
--postfix .unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 1.0 \
--fracrepetitions 2 \
--algrepetitions 1 \
--testfraction 0.0001 \
--run_fraction_parallelism sequential \
--run_fracrepetition_parallelism distributed \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_set_5.Proteins.SmallTest.out &


nohup julia \
--threads 20 --proc 5 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .SmallTestB \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 0.25 \
--fracrepetitions 2 \
--testfraction 0.01 \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_set_5.Proteins.SmallTestB.out &


# medium test

nohup julia \
--threads 20 --proc 2 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix reads.sorted.aligned.filtered. \
--postfix .unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 0.5 \
--fracrepetitions 4 \
--algrepetitions 2 \
--testfraction 0.2 \
--run_fraction_parallelism sequential \
--run_fracrepetition_parallelism distributed \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_set_5.Proteins.MediumTest.out &


nohup julia \
--threads 20 --proc 5 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .MediumTestB \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 0.5 \
--fracrepetitions 8 \
--testfraction 0.2 \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_set_5.Proteins.MediumTestB.out &


# medium test 2

nohup julia \
--threads 40 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix reads.sorted.aligned.filtered. \
--postfix .unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 0.5 \
--fracrepetitions 4 \
--algrepetitions 2 \
--testfraction 0.4 \
--run_fraction_parallelism sequential \
--run_fracrepetition_parallelism distributed \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_set_5.Proteins.MediumTest2.out &


# medium test 3

nohup julia \
--threads 40 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix reads.sorted.aligned.filtered. \
--postfix .unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample \
--fracstep 0.25 \
--fracrepetitions 4 \
--algrepetitions 2 \
--testfraction 0.6 \
--run_fraction_parallelism sequential \
--run_fracrepetition_parallelism distributed \
> D.pealeii/MpileupAndTranscripts/IlluminaExample/maximal_independent_set_5.Proteins.MediumTest3.out &
```





### Ruti's Illumina PE READS, take 2 (shortened ids)

#### Example

##### Pileup

```bash
mkdir -p D.pealeii/MpileupAndTranscripts/IlluminaExample2

nohup \
python Code/pileup.py \
--transcriptome D.pealeii/Annotations/orfs_squ.fa \
--data_table D.pealeii/Alignment/Illumina/reads.ByChrom/data_table.example.csv \
--known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
--out_dir D.pealeii/MpileupAndTranscripts/IlluminaExample2 \
--top_x_noisy_positions 3 \
--assurance_factor 1.5 \
--include_flags 3 \
--exclude_flags 2304 \
--min_bq 30 \
--parity PE \
--processes 4 \
--threads 20 \
--keep_pileup_files \
> D.pealeii/MpileupAndTranscripts/IlluminaExample2/pileup.example.3.7.22.out & 
```
* alu 16
* 3.7.22
* 13:00
 
##### Julia

```
INFILES=$(echo D.pealeii/MpileupAndTranscripts/IlluminaExample2/*.unique_proteins.csv)
echo $INFILES

JULIA_PROJECT=.

nohup julia \
--threads 30 --proc 5 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .0_25.2.3.0_01 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample2 \
--fracstep 0.25 \
--fracrepetitions 2 \
--algrepetitions 3 \
--testfraction 0.01 \
--algs Ascending Descending Unordered \
> D.pealeii/MpileupAndTranscripts/IlluminaExample2/maximal_independent_set_5.Proteins.0_25.2.3.0_01.out &


nohup julia \
--threads 50 --proc 2 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .Full \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample2 \
--fracstep 0.25 \
--fracrepetitions 2 \
--algrepetitions 2 \
--algs Ascending Descending Unordered \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/IllumiqnaExample2/maximal_independent_set_5.Proteins.Full.out &

nohup julia \
--threads 60 --proc 4 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--postfix_to_add .Full2 \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/IlluminaExample2 \
--fracstep 0.1 \
--fracrepetitions 5 \
--algrepetitions 3 \
--algs Ascending Descending \
--run_solve_threaded \
> D.pealeii/MpileupAndTranscripts/IllumiqnaExample2/maximal_independent_set_5.Proteins.Full2.out &
```


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
* alu 13
* 12.7.22
* 10:07

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
* alu 13
* 12.7.22
* 17:54

Running again, using a single core only, on files that weren't completed, or that weren't reached before the main process was killed.

```
<!-- INFILES=$(echo D.pealeii/MpileupAndTranscripts/Illumina/*.unique_proteins.csv)
echo $INFILES -->


INFILES="\
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp140712_c0_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141574_c0_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141532_c3_seq11.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141882_c0_seq14.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141880_c1_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141565_c6_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp140987_c3_seq1.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp136058_c0_seq1.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141378_c0_seq7.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141684_c0_seq1.unique_proteins.csv
"

echo $INFILES

JULIA_PROJECT=.

nohup julia \
--threads 50 \
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
> D.pealeii/MpileupAndTranscripts/Illumina/maximal_independent_set_5.Proteins.Addendum.out &
```
* alu 13
* 17.7.22
* 19:35



Continuing rerunning, using a single core only, on files that weren't completed, or that weren't reached before the main process was killed.
Now, without `run_solve_threaded` and without comp140987_c3_seq1 that was completed in the previous try.

```
INFILES="\
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp140712_c0_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141574_c0_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141532_c3_seq11.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141882_c0_seq14.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141880_c1_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141565_c6_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp136058_c0_seq1.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141378_c0_seq7.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141684_c0_seq1.unique_proteins.csv"

echo $INFILES

JULIA_PROJECT=.

nohup julia \
--threads 50 \
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
> D.pealeii/MpileupAndTranscripts/Illumina/maximal_independent_set_5.Proteins.Addendum2.out &
```
* alu 13
* 18.7.22
* 13:55


```
INFILES="\
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp140712_c0_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141574_c0_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141532_c3_seq11.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141882_c0_seq14.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141880_c1_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141565_c6_seq3.unique_proteins.csv \
D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141684_c0_seq1.unique_proteins.csv"

echo $INFILES

JULIA_PROJECT=.

nohup julia \
--threads 50 \
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
> D.pealeii/MpileupAndTranscripts/Illumina/maximal_independent_set_5.Proteins.Addendum3.out &
```
* alu 13
* 18.7.22
* 13:55



















   








