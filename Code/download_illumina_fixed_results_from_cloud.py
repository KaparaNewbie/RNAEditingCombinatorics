import subprocess
from pathlib import Path

from icecream import ic

out_dir = Path(
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina",
    "FinalFixedExpressionFromCloud",
)
out_dir.mkdir(exist_ok=True)


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


for condition in conditions:
    print(f"Downloading files for condition: {condition}")
    gs_bucket_path = f"gs://kobi-rna-comb-bucket/main/results/{condition}"

    subprocess.run(["gsutil", "cp", "-r", gs_bucket_path, out_dir])

    print(f"Finished downloading files for condition: {condition}\n")
