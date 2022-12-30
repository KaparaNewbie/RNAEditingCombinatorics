from pathlib import Path

import pandas as pd


GFF_COLS = [
    "Chrom",
    "Source",
    "Type",
    "Start",
    "End",
    "Score",
    "Strand",
    "Phase",
    "Attributes",
]


def read_gff(gff_file: Path, names: list[str] = GFF_COLS):
    gff_df = pd.read_csv(gff_file, sep="\t", names=names, comment="#")
    # gff_df["Start"] = gff_df["Start"].astype(int)
    # gff_df["End"] = gff_df["End"].astype(int)
    return gff_df


# def attributes_to_dict(attributes):
#     """
#     Parse a string of a gff file into a dictionary.
#     @param attributes: attributes of a single feature (one line) in a gff file, as they appear in the 9th column
#     @type attributes: str
#     @return: attributes
#     @rtype: dict
#     """
#     attributes = [attribute.split("=") for attribute in attributes.split(";")]
#     attributes = {
#         attribute[0]: attribute[1] for attribute in attributes
#     }  # {k: v for (k, v) in attributes}
#     return attributes


def types_terms_in_gff3_file(gff3_file):
    """Return a set of the terms found in the 'type' column (the 3rd column) of a GFF3 file.
    A GFF3 is a tab-delimited 9-cols file for storing genomic features.
    Documentation:
    (1)  https://www.ensembl.org/info/website/upload/gff3.html
    (2)  http://gmod.org/wiki/GFF3"""
    gff_df = read_gff(gff3_file)
    return set(gff_df["Type"].unique())



def divide_to_files_by_type(gff_file, out_dir, prefix=""):
    """Divide a gff3_file to sub files by their type (3rd col), and write them to out_dir."""
    gff_types = types_terms_in_gff3_file(gff_file)
    
    gff_df = read_gff(gff_file)
    
    type_gff_dfs = {
        gff_type: gff_df.loc[gff_df["Type"] == gff_type]
        for gff_type in gff_types
    }
    
    for gff_type in gff_types:
        out_file = Path(out_dir, f"{prefix}{gff_type}.gff")
        gff_type_df = type_gff_dfs[gff_type]
        gff_type_df.to_csv(out_file, sep="\t", index=False, header=False)
    

def sub_gff_dfs_by_genes(gff_df: pd.DataFrame) -> list[pd.DataFrame]:
    """Divide gff_df to multiple dataframes, where each sub df contains one gene alongside its exons, introns, etc."""
    genes_indices = gff_df.loc[gff_df["Type"] == "gene"].index
    sub_gff_dfs = []
    for i in range(len(genes_indices)):
        start = genes_indices[i]
        if i < len(genes_indices) - 1:
            end = genes_indices[i + 1]
        else:
            end = len(gff_df)
        sub_gff_df = gff_df.iloc[start:end]
        sub_gff_dfs.append(sub_gff_df)
    return sub_gff_dfs
