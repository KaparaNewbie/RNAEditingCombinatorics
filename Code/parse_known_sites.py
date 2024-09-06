import argparse
import pandas as pd

from General.argparse_utils import abs_path_from_str

parser = argparse.ArgumentParser()
parser.add_argument("--editing_sites_excel_file", type=abs_path_from_str, help="Known sites file downloaded from https://www.tau.ac.il/~elieis/squid/")
parser.add_argument("--sheet_name", help="For example, `D.pea Editing sites` or `O.vul Editing sites`")
parser.add_argument("--out_csv_file", type=abs_path_from_str, help="Save specifc sheet_name in a separate csv file")
parser.add_argument("--out_bed_file", type=abs_path_from_str, help="Parse editing sites to a bed file in this path")
parser.add_argument("--ignore_score", action="store_true", help="Set score to `.` instead of given editing frequency")
args = parser.parse_args()

df = pd.read_excel(args.editing_sites_excel_file, sheet_name=args.sheet_name)
df.to_csv(args.out_csv_file, index=False)
if args.ignore_score:

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
else:
    bed_cols = ["Trinity id", "Editing location (base1)", "SwissProt name", "Editing Level"]
    df2 = df.filter(bed_cols).rename(
        columns={
            "Trinity id": "#TrinityID",
            "Editing location (base1)": "End",
            "SwissProt name": "SwissProtName",
            "Editing Level": "EditingLevel",
        }
    )
    df2.insert(1, "Start", df2["End"] - 1)
    df2["Strand"] = "+"

df2.to_csv(args.out_bed_file, index=False, sep="\t")