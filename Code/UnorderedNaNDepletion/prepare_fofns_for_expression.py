"""
Prepare files of file-names needed to run `expressionlevels.jl`.
Useful for undirected seq data where there are many transcripts with sufficient mapping.
"""


from pathlib import Path
import sys
import argparse

sys.path.append(str(Path("Code").absolute()))
from General.argparse_utils import abs_path_from_str, expanded_path_from_str


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Prepare files of file-names needed to run `expressionlevels.jl`.",
)
parser.add_argument(
    "--proteins_dir",
    required=True,
    type=abs_path_from_str,
    help="A folder with files of unique proteins.",
)
parser.add_argument(
    "--distinct_proteins_dir",
    # required=True,
    type=abs_path_from_str,
    default=None,
    help="A folder with files with distinct isoforms. Defaults to `protein_dir`.",
)
parser.add_argument(
    "--distinct_proteins_postfix", default="DistinctUniqueProteins.*.csv"
)
parser.add_argument(
    "--out_dir",
    type=abs_path_from_str,
    default=None,
    help="Write the fofns to `out_dir`. Defaults to `distinct_proteins_dir`.",
)

args = parser.parse_args()

proteins_dir = args.proteins_dir
distinct_proteins_dir = args.distinct_proteins_dir
distinct_proteins_postfix = args.distinct_proteins_postfix
out_dir = args.out_dir

distinct_proteins_dir = (
    distinct_proteins_dir if distinct_proteins_dir is not None else proteins_dir
)
out_dir = out_dir if out_dir is not None else distinct_proteins_dir

distinct_proteins_files = list(
    distinct_proteins_dir.glob(f"*{distinct_proteins_postfix}")
)
unique_proteins_files = [
    Path(
        proteins_dir,
        f"{distinct_proteins_files.name.split('.')[0]}.unique_proteins.csv",
    )
    for distinct_proteins_files in distinct_proteins_files
]
chroms_names = [
    distinct_proteins_files.name.split(".")[0]
    for distinct_proteins_files in distinct_proteins_files
]

distinct_proteins_list = Path(out_dir, "DistinctProteinsForExpressionLevels.txt")
unique_proteins_list = Path(out_dir, "UniqueProteinsForExpressionLevels.txt")
chroms_names_list = Path(out_dir, "ChromsNamesForExpressionLevels.txt")

lists = [distinct_proteins_files, unique_proteins_files, chroms_names]
lists_file_names = [distinct_proteins_list, unique_proteins_list, chroms_names_list]

for _list, list_file_name in zip(lists, lists_file_names):
    with list_file_name.open("w") as list_file_name:
        list_file_name.write(" ".join(map(str, _list)))
