import subprocess

from icecream import ic

conditions_to_ignore = ["RUSC2", "RIMS2"]

big_chroms = [
    "comp140712_c0_seq3",
    "comp141532_c3_seq11",
    "comp141565_c6_seq3",
    "comp141574_c0_seq3",
    "comp141684_c0_seq1",
    "comp141880_c1_seq3",
    "comp141882_c0_seq14",
]

conditions = [
    "RUSC2",
    "TRIM2",
    "CA2D3",
    "ABL",
    "DGLA",
    "K0513",
    "KCNAS",
    "ACHA4",
    "ANR17",
    "TWK7",
    "SCN1",
    "CACB2",
    "RIMS2",
    "PCLO",
    "DOP1",
    "IQEC1",
    "CSKI1",
    "MTUS2",
    "ROBO2",
]

chroms = [
    "comp141881_c0_seq3",
    "comp141044_c0_seq2",
    "comp140439_c0_seq1",
    "comp126362_c0_seq1",
    "comp141517_c0_seq1",
    "comp141840_c0_seq2",
    "comp141640_c0_seq1",
    "comp140987_c3_seq1",
    "comp140910_c2_seq1",
    "comp136058_c0_seq1",
    "comp141378_c0_seq7",
    "comp141158_c1_seq2",
    "comp140712_c0_seq3",
    "comp141882_c0_seq14",
    "comp141880_c1_seq3",
    "comp141565_c6_seq3",
    "comp141684_c0_seq1",
    "comp141532_c3_seq11",
    "comp141574_c0_seq3",
]

reads_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.reads.csv"
    for chrom in chroms
]

unique_proteins_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_proteins.csv"
    for chrom in chroms
]

distinct_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.DistinctUniqueProteins.12.07.2022-20:54:38.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.DistinctUniqueProteins.13.07.2022-06:33:23.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.DistinctUniqueProteins.12.07.2022-22:51:22.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.DistinctUniqueProteins.15.07.2022-06:11:18.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.DistinctUniqueProteins.14.07.2022-07:43:15.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.DistinctUniqueProteins.13.07.2022-20:30:25.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.DistinctUniqueProteins.12.07.2022-19:44:02.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.DistinctUniqueProteins.18.07.2022-07:50:43.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.DistinctUniqueProteins.13.07.2022-16:15:35.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.DistinctUniqueProteins.21.07.2022-07:57:53.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.DistinctUniqueProteins.19.07.2022-08:12:24.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.DistinctUniqueProteins.13.07.2022-01:54:59.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140712_c0_seq3.DistinctUniqueProteins.05.11.2022-02_11_33.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141882_c0_seq14.DistinctUniqueProteins.11.11.2022-22_23_54.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141880_c1_seq3.DistinctUniqueProteins.14.11.2022-20_58_20.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141565_c6_seq3.DistinctUniqueProteins.16.11.2022-12_42_59.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141684_c0_seq1.DistinctUniqueProteins.19.11.2022-01_51_21.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141532_c3_seq11.DistinctUniqueProteins.05.12.2022-20_38_20.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141574_c0_seq3.DistinctUniqueProteins.12.12.2022-19_31_01.csv",
]


for (
    chrom,
    condition,
    reads_file,
    unique_proteins_file,
    distinct_unique_proteins_file,
) in zip(
    chroms,
    conditions,
    reads_files,
    unique_proteins_files,
    distinct_unique_proteins_files,
):
    if condition in conditions_to_ignore:
        ic(condition)
        continue

    print(f"Uploading files for condition: {condition}")

    # gs_bucket_path = f"gs://kobi-rna-comb-bucket/main/Illumina/{condition}/"
    sample_size_interfix = "BigSamples" if chrom in big_chroms else "SmallSamples"
    gs_bucket_path = (
        f"gs://kobi-rna-comb-bucket/main/Illumina.{sample_size_interfix}/{condition}/"
    )

    for file in [distinct_unique_proteins_file, unique_proteins_file, reads_file]:

        subprocess.run(
            [
                "gsutil",
                "cp",
                file,
                gs_bucket_path,
            ]
        )

    print(f"Finished uploading files for condition: {condition}\n")
