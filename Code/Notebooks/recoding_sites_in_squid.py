# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Recoding potential at annotated squid RNA-editing sites
#
# One row per triplet in each BED-defined transcript CDS. Only **annotated editing sites** whose coding-strand reference base is A may undergo A→G substitution (inosine is read as G). Enumerate every nonempty subset of those sites, without assuming that the edits actually co-occur or have equal probability.
#
# Use the [NCBI standard genetic code (table 1)](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1). This measures codon translation changes, not initiation efficiency, RNA structure, or editing probability. BED coordinates are zero-based, end-exclusive; codon indices and positions within a codon are one-based. Reverse-complement negative-strand CDSs before translation. The CDS BED must contain one contiguous, in-frame interval per transcript; it is not a genomic exon BED.
#
# BED names are homology labels, not established unique squid gene IDs. Keep transcript IDs separate, including transcripts sharing a label. Reference stop codons (`*`) remain in the codon table; stop loss is reported separately from the number of recodable amino acids.

# %%
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
from IPython.display import display

# Locate the project from either the notebook directory or the project root.
project_root = next(
    p for p in (Path.cwd(), *Path.cwd().parents)
    if (p / "D.pealeii/Annotations/Jan2025/orfs_squ.fa").is_file()
)
annotation_dir = project_root / "D.pealeii/Annotations/Jan2025"
orfs_fasta_file = annotation_dir / "orfs_squ.fa"
orfs_bed_file = annotation_dir / "orfs_squ.bed"
editing_sites_bed_file = annotation_dir / "D.pea.EditingSites.bed"
output_dir = annotation_dir / "RecodingPotential"

# %% [markdown]
# ## Classification
#
# - `not_editable`: no usable annotated A in the codon.
# - `editable_not_recoded`: every permitted nonempty edit combination preserves translation.
# - `always_recoded_if_edited`: every permitted nonempty edit combination changes translation.
# - `combination_dependent_recoding`: some combinations preserve translation and others change it.
#
# `min_edits_for_recoding` is the minimum over all translation-changing combinations; `requires_multiple_edits_for_recoding` means that minimum is at least two. Separately, `min_edits_by_recoded_amino_acid` records the minimum for **each target**: a particular amino acid may require two edits even when another is reachable with one. `editing_outcomes` records every combination as JSON, including synonymous outcomes. For stop codons these fields describe stop loss; the transcript amino-acid count excludes them.

# %%
# NCBI table 1, in T/C/A/G order at each codon position.
CODONS = ["".join(bases) for bases in product("TCAG", repeat=3)]
GENETIC_CODE = dict(zip(
    CODONS,
    "FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG",
))
COMPLEMENT = str.maketrans("ACGT", "TGCA")

def classify_codon(codon, editable_positions):
    """Positions are one-based and refer only to usable annotated coding-strand As."""
    positions = tuple(sorted(set(editable_positions)))
    if codon not in GENETIC_CODE:
        raise ValueError(f"Invalid codon: {codon!r}")
    if any(p not in (1, 2, 3) or codon[p - 1] != "A" for p in positions):
        raise ValueError(f"Editable positions must be As: {codon}, {positions}")
    ref_aa = GENETIC_CODE[codon]
    outcomes = []
    for count in range(1, len(positions) + 1):
        for edited_positions in combinations(positions, count):
            edited = list(codon)
            for position in edited_positions:
                edited[position - 1] = "G"
            edited = "".join(edited)
            aa = GENETIC_CODE[edited]
            outcomes.append({
                "positions": list(edited_positions), "codon": edited,
                "amino_acid": aa, "recoded": aa != ref_aa,
            })
    recoded = [o for o in outcomes if o["recoded"]]
    synonymous = [o for o in outcomes if not o["recoded"]]
    minimum_by_aa = {}
    for o in recoded:
        aa = o["amino_acid"]
        minimum_by_aa[aa] = min(minimum_by_aa.get(aa, 4), len(o["positions"]))
    minimum = min(minimum_by_aa.values(), default=None)
    if not outcomes:
        category = "not_editable"
    elif not recoded:
        category = "editable_not_recoded"
    elif not synonymous:
        category = "always_recoded_if_edited"
    else:
        category = "combination_dependent_recoding"
    compact_json = lambda value: json.dumps(value, separators=(",", ":"), sort_keys=True)
    return {
        "reference_codon": codon, "reference_amino_acid": ref_aa,
        "annotated_editable_positions_1based": compact_json(positions),
        "n_editable_sites": len(positions), "category": category,
        "can_be_edited": bool(outcomes), "can_be_recoded": bool(recoded),
        "can_be_edited_without_recoding": bool(synonymous),
        "always_recoded_if_edited": bool(outcomes) and not synonymous,
        "n_editing_combinations": len(outcomes),
        "n_recoding_combinations": len(recoded),
        "n_synonymous_combinations": len(synonymous),
        "min_edits_for_recoding": minimum,
        "requires_multiple_edits_for_recoding": minimum is not None and minimum >= 2,
        "has_target_requiring_multiple_edits": any(v >= 2 for v in minimum_by_aa.values()),
        "recoded_amino_acids": compact_json(sorted(minimum_by_aa)),
        "min_edits_by_recoded_amino_acid": compact_json(minimum_by_aa),
        "editing_outcomes": compact_json(outcomes),
        "is_reference_stop": ref_aa == "*",
        "can_recode_amino_acid": ref_aa != "*" and bool(recoded),
        "can_lose_stop": ref_aa == "*" and bool(recoded),
    }

# Precompute all valid codon / annotated-site-mask combinations.
# A bit set at position p-1 means that base p is an annotated editable A.
pattern_rows = []
for codon in CODONS:
    for mask in range(8):
        positions = [p for p in (1, 2, 3) if mask & (1 << (p - 1))]
        if all(codon[p - 1] == "A" for p in positions):
            pattern_rows.append({"pattern_key": codon + str(mask), **classify_codon(codon, positions)})
codon_patterns_df = pd.DataFrame(pattern_rows)
codon_patterns_df["min_edits_for_recoding"] = codon_patterns_df["min_edits_for_recoding"].astype("Int8")
pattern_ids = {key: i for i, key in enumerate(codon_patterns_df["pattern_key"])}
display(codon_patterns_df.loc[
    codon_patterns_df["pattern_key"].isin(["CCC0", "GCA4", "AAG1", "AAA7", "TAA6"]),
    ["reference_codon", "annotated_editable_positions_1based", "category", "min_edits_for_recoding", "min_edits_by_recoded_amino_acid"],
])


# %%
def read_fasta(path):
    sequences = {}
    name, pieces = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(pieces).upper()
                name = line[1:].split()[0]
                if name in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                pieces = []
            else:
                if name is None:
                    raise ValueError("Sequence before first FASTA header")
                pieces.append(line)
    if name is not None:
        sequences[name] = "".join(pieces).upper()
    if not sequences:
        raise ValueError("Empty FASTA")
    return sequences

BED_COLUMNS = ["transcript_id", "cds_start_0based", "cds_end_0based", "gene_label", "score", "strand"]

def read_bed6(path):
    # Both headerless CDS BED and the #TrinityID-prefixed editing BED are supported.
    frame = pd.read_csv(path, sep="\t", comment="#", header=None)
    if frame.shape[1] != 6:
        raise ValueError(f"Expected BED6: {path}")
    frame.columns = BED_COLUMNS
    for column in ["cds_start_0based", "cds_end_0based"]:
        values = pd.to_numeric(frame[column], errors="raise")
        if values.isna().any() or (values % 1 != 0).any():
            raise ValueError(f"Noninteger BED coordinates: {path}")
        frame[column] = values.astype("int64")
    if not frame["strand"].isin(["+", "-"]).all():
        raise ValueError(f"BED strand must be + or -: {path}")
    return frame

def extract_cds(sequence, start, end, strand):
    if strand not in ("+", "-") or not 0 <= start < end <= len(sequence):
        raise ValueError(f"Invalid CDS bounds/strand: {start}, {end}, {strand}")
    cds = sequence[start:end]
    if len(cds) % 3:
        raise ValueError("CDS length must be divisible by three; no silent trimming")
    if set(cds) - set("ACGT"):
        raise ValueError("Ambiguous CDS bases; cannot assign exact recoding potential")
    return cds if strand == "+" else cds.translate(COMPLEMENT)[::-1]

def coding_offset(position, start, end, strand):
    return position - start if strand == "+" else end - 1 - position

def transcript_codon_start(offsets, start, end, strand):
    return start + offsets if strand == "+" else end - offsets - 3

transcript_sequences = read_fasta(orfs_fasta_file)
cds_bed_df = read_bed6(orfs_bed_file).drop(columns="score")
if cds_bed_df["transcript_id"].duplicated().any():
    raise ValueError("Expected exactly one contiguous CDS interval per transcript")
missing = set(cds_bed_df["transcript_id"]) - transcript_sequences.keys()
if missing:
    raise ValueError(f"CDS transcripts absent from FASTA: {sorted(missing)[:10]}")
cds_sequences = {}
for row in cds_bed_df.itertuples(index=False):
    try:
        cds_sequences[row.transcript_id] = extract_cds(
            transcript_sequences[row.transcript_id], row.cds_start_0based,
            row.cds_end_0based, row.strand,
        )
    except ValueError as exc:
        raise ValueError(f"{row.transcript_id}: {exc}") from exc
print(f"Loaded {len(cds_sequences):,} validated transcript CDSs.")

# %% [markdown]
# ## Audit annotated sites
#
# The editing BED is in FASTA coordinates. A valid negative-strand site has T in the stored FASTA and A after reverse complementation. Sites outside the CDS, absent from this FASTA/CDS set, on an inconsistent strand, or with a non-A coding reference are excluded with explicit reasons. Duplicate usable coordinates contribute only once. No unannotated A is introduced as an editable site.

# %%
editing_sites_df = read_bed6(editing_sites_bed_file).rename(columns={
    "cds_start_0based": "site_start_0based", "cds_end_0based": "site_end_0based",
    "gene_label": "site_gene_label", "strand": "site_strand",
}).drop(columns="score")
site_audit_df = editing_sites_df.merge(cds_bed_df, on="transcript_id", how="left", validate="many_to_one")
site_status, reference_bases = [], []
seen_sites = set()
editable_offsets = {transcript: set() for transcript in cds_sequences}
for row in site_audit_df.itertuples(index=False):
    seq = transcript_sequences.get(row.transcript_id)
    position = row.site_start_0based
    base = seq[position] if seq is not None and 0 <= position < len(seq) else None
    reference_bases.append(base)
    if seq is None:
        status = "transcript_missing_from_fasta"
    elif pd.isna(row.cds_start_0based):
        status = "transcript_missing_from_cds_bed"
    elif row.site_end_0based != position + 1:
        status = "not_single_base_interval"
    elif not 0 <= position < len(seq):
        status = "outside_transcript"
    elif row.site_strand != row.strand:
        status = "strand_mismatch"
    elif not row.cds_start_0based <= position < row.cds_end_0based:
        status = "outside_cds"
    elif base != ("A" if row.strand == "+" else "T"):
        status = "reference_not_editable_A"
    elif (row.transcript_id, position) in seen_sites:
        status = "duplicate_site"
    else:
        status = "used"
        seen_sites.add((row.transcript_id, position))
        editable_offsets[row.transcript_id].add(int(coding_offset(
            position, row.cds_start_0based, row.cds_end_0based, row.strand,
        )))
    site_status.append(status)
site_audit_df["fasta_reference_base"] = reference_bases
site_audit_df["status"] = site_status
display(site_audit_df["status"].value_counts().rename_axis("status").to_frame("sites"))
excluded_sites_df = site_audit_df.loc[site_audit_df["status"] != "used"].copy()
if not excluded_sites_df.empty:
    warnings.warn(f"Excluded {len(excluded_sites_df):,} annotated rows; see excluded_sites_df and the exported audit.")
    display(excluded_sites_df.head(10))

# %%
# Small independent examples cover annotation restriction, subset effects, and orientation.
assert classify_codon("AAA", [])['category'] == "not_editable"
assert classify_codon("AAA", [3])['category'] == "editable_not_recoded"
assert classify_codon("AAG", [1])['category'] == "always_recoded_if_edited"
assert classify_codon("AAA", [1, 2, 3])['category'] == "combination_dependent_recoding"
assert json.loads(classify_codon("AAA", [1, 2, 3])["min_edits_by_recoded_amino_acid"]) == {"E": 1, "G": 2, "R": 1}
assert classify_codon("TAA", [2, 3])["requires_multiple_edits_for_recoding"]
assert classify_codon("TAA", [2])["category"] == "editable_not_recoded"
assert not classify_codon("TAA", [2, 3])["can_recode_amino_acid"]
assert classify_codon("TAA", [2, 3])["can_lose_stop"]
# The internal CDS is TTTCAT on the FASTA strand and ATGAAA on the coding strand.
assert extract_cds("CCTTTCATGG", 2, 8, "-") == "ATGAAA"
assert coding_offset(4, 2, 8, "-") == 3
assert transcript_codon_start(np.array([0, 3]), 2, 8, "-").tolist() == [5, 2]
for invalid_seq in ["ATGA", "ATN"]:
    try:
        extract_cds(invalid_seq, 0, len(invalid_seq), "+")
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid CDS was silently accepted")
for row in codon_patterns_df.itertuples(index=False):
    assert row.n_editing_combinations == 2 ** row.n_editable_sites - 1
    assert row.n_editing_combinations == row.n_recoding_combinations + row.n_synonymous_combinations
    for outcome in json.loads(row.editing_outcomes):
        differences = [i + 1 for i, (a, b) in enumerate(zip(row.reference_codon, outcome["codon"])) if a != b]
        assert differences == outcome["positions"]
        assert all(row.reference_codon[p - 1] == "A" and outcome["codon"][p - 1] == "G" for p in differences)
print("Codon classification and coordinate checks passed.")

# %% [markdown]
# ## One row per CDS triplet
#
# `cds_offset_0based` counts from the first base of the oriented CDS. `transcript_start_0based` and `transcript_end_0based` locate the triplet in the original FASTA. On the negative strand, codon position 1 is `transcript_end_0based - 1`; on the positive strand it is `transcript_start_0based`. Add/subtract the one-based positions in `annotated_editable_positions_1based` to recover the exact editing coordinates.
#
# Repeated sequence/outcome fields use pandas categorical storage to keep the full table manageable.

# %%
transcript_index_parts, offset_parts, pattern_parts = [], [], []
for transcript_index, row in enumerate(cds_bed_df.itertuples(index=False)):
    cds = cds_sequences[row.transcript_id]
    n_codons = len(cds) // 3
    masks = np.zeros(n_codons, dtype=np.uint8)
    for offset in editable_offsets[row.transcript_id]:
        masks[offset // 3] |= 1 << (offset % 3)
    patterns = np.fromiter(
        (pattern_ids[cds[i * 3:i * 3 + 3] + str(int(masks[i]))] for i in range(n_codons)),
        dtype=np.uint8, count=n_codons,
    )
    transcript_index_parts.append(np.full(n_codons, transcript_index, dtype=np.int32))
    offset_parts.append(np.arange(n_codons, dtype=np.int32) * 3)
    pattern_parts.append(patterns)
transcript_indices = np.concatenate(transcript_index_parts)
offsets = np.concatenate(offset_parts)
pattern_indices = np.concatenate(pattern_parts)
del transcript_index_parts, offset_parts, pattern_parts

codons_df = pd.DataFrame({
    "transcript_id": pd.Categorical.from_codes(transcript_indices, cds_bed_df["transcript_id"]),
    "codon_index_1based": offsets // 3 + 1,
    "cds_offset_0based": offsets,
})
for column in ["gene_label", "strand"]:
    values = pd.Categorical(cds_bed_df[column])
    codons_df[column] = pd.Categorical.from_codes(values.codes[transcript_indices], values.categories)
for column in ["cds_start_0based", "cds_end_0based"]:
    codons_df[column] = cds_bed_df[column].to_numpy()[transcript_indices]
starts = np.where(
    codons_df["strand"] == "+",
    codons_df["cds_start_0based"].to_numpy() + offsets,
    codons_df["cds_end_0based"].to_numpy() - offsets - 3,
)
codons_df["transcript_start_0based"] = starts
codons_df["transcript_end_0based"] = starts + 3
for column in codon_patterns_df.columns.drop("pattern_key"):
    source = codon_patterns_df[column]
    if source.dtype == object:
        values = pd.Categorical(source)
        codons_df[column] = pd.Categorical.from_codes(values.codes[pattern_indices], values.categories)
    else:
        codons_df[column] = source.array.take(pattern_indices)
del starts, transcript_indices, offsets, pattern_indices
print(f"Built {len(codons_df):,} codon rows; memory {codons_df.memory_usage(deep=True).sum() / 1e6:,.0f} MB.")
display(codons_df.head())
display(codons_df["category"].value_counts().rename_axis("category").to_frame("codons"))

# %%
transcript_summary_df = codons_df.groupby("transcript_id", observed=True, sort=False).agg(
    n_cds_codons=("reference_codon", "size"),
    n_reference_stop_codons=("is_reference_stop", "sum"),
    n_annotated_editable_sites=("n_editable_sites", "sum"),
    n_editable_codons=("can_be_edited", "sum"),
    n_recodable_codons_including_stops=("can_be_recoded", "sum"),
    n_recodable_amino_acids=("can_recode_amino_acid", "sum"),
    n_stop_loss_codons=("can_lose_stop", "sum"),
    n_codons_requiring_multiple_edits=("requires_multiple_edits_for_recoding", "sum"),
    n_codons_with_target_requiring_multiple_edits=("has_target_requiring_multiple_edits", "sum"),
).reset_index()
category_counts = pd.crosstab(codons_df["transcript_id"], codons_df["category"]).add_prefix("n_").reset_index()
transcript_summary_df = cds_bed_df.merge(transcript_summary_df, on="transcript_id", validate="one_to_one").merge(
    category_counts, on="transcript_id", validate="one_to_one",
)
transcript_summary_df["n_amino_acids"] = transcript_summary_df["n_cds_codons"] - transcript_summary_df["n_reference_stop_codons"]
transcript_summary_df["fraction_amino_acids_recodable"] = (
    transcript_summary_df["n_recodable_amino_acids"] / transcript_summary_df["n_amino_acids"].replace(0, np.nan)
)
# Validate completeness and that every retained editing site appears exactly once.
assert len(codons_df) == ((cds_bed_df["cds_end_0based"] - cds_bed_df["cds_start_0based"]) // 3).sum()
assert len(transcript_summary_df) == len(cds_bed_df)
assert codons_df["n_editable_sites"].sum() == (site_audit_df["status"] == "used").sum()
assert (transcript_summary_df["n_cds_codons"] == (transcript_summary_df["cds_end_0based"] - transcript_summary_df["cds_start_0based"]) // 3).all()
assert (transcript_summary_df["n_recodable_amino_acids"] <= transcript_summary_df["n_amino_acids"]).all()
display(transcript_summary_df.sort_values("n_recodable_amino_acids", ascending=False).head(20))
print(f"Recodable amino-acid positions across transcripts: {transcript_summary_df['n_recodable_amino_acids'].sum():,}")
print(f"Transcripts with at least one recodable amino acid: {(transcript_summary_df['n_recodable_amino_acids'] > 0).sum():,}")

# %%
transcript_summary_df

# %%
transcript_summary_df["n_recodable_codons_including_stops"].describe().round(2)

# %%
transcript_summary_df.loc[
    transcript_summary_df["n_recodable_codons_including_stops"].eq(
        transcript_summary_df["n_recodable_codons_including_stops"].max()
    )
]

# %%
transcript_summary_df.loc[
    transcript_summary_df["transcript_id"].eq(
        "comp141882_c0_seq14"
    )
]

# %% [markdown]
# ## Export
#
# The compressed TSV contains every CDS triplet, including uneditable codons and reference stops. JSON-valued fields can be decoded with `json.loads`. The transcript summary gives the requested amino-acid counts while preserving the BED gene label; counts across transcripts are not deduplicated counts of genomic gene positions. The audit includes every input editing-site row and its inclusion/exclusion status. The small pattern table provides a convenient way to inspect combination-level effects without loading all codons.

# %%
output_dir.mkdir(parents=True, exist_ok=True)
codon_table_file = output_dir / "squid_cds_codon_recoding.tsv.gz"
transcript_summary_file = output_dir / "squid_transcript_recoding_summary.tsv"
site_audit_file = output_dir / "squid_editing_site_audit.tsv"
pattern_table_file = output_dir / "squid_codon_editing_patterns.tsv"
codons_df.to_csv(codon_table_file, sep="\t", index=False, chunksize=100_000,
                 compression={"method": "gzip", "compresslevel": 1})
transcript_summary_df.to_csv(transcript_summary_file, sep="\t", index=False)
site_audit_df.to_csv(site_audit_file, sep="\t", index=False)
codon_patterns_df.to_csv(pattern_table_file, sep="\t", index=False)
for path in [codon_table_file, transcript_summary_file, site_audit_file, pattern_table_file]:
    print(f"{path.name}: {path.stat().st_size / 1e6:,.2f} MB")
