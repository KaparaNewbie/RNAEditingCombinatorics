#!/usr/bin/env python
"""Find the largest extension of each FIXED, saved protein independent set.

Scientific question
-------------------
Each vertex is one original unique protein pattern. At each amino-acid position
the pattern contains a set of possibilities, such as {K, R}. Two vertices have
an edge when their possibility sets overlap at EVERY position: the patterns are
indistinguishable. One disjoint position makes them distinguishable (no edge).
An independent set therefore consists of pairwise distinguishable patterns.

Let A be the vertices available to a saved solution, S its selected vertices,
and N(v) the neighbors of vertex v. First verify that S is independent. Then
define C = A minus (S union all N(s), for s in S). Every member of C can be added
individually to S, but members of C may conflict with each other. Find a largest
independent subset T of C. The answer is |T|, and the extended size is |S| + |T|.

S is maximal exactly when C is empty: no additional vertex can be included.
Maximum instead means largest cardinality among ALL independent sets of a graph.
Here only the graph induced by C is optimized. Every original member of S stays
fixed; none is removed or replaced. Any feasible extension must lie in C, so an
exact maximum independent set of C gives the largest possible extension of S.
It does not establish the global maximum independent set of the original graph.

Entropy filtering belongs to the separate expression/read-reassignment analysis.
It would change this diagnostic's original graph and is never applied here.

Workflow and provenance
-----------------------
1. Read the explicit 24-cohort squid configuration and saved Figure 6 octopus
   populations. Match canonical rows by protein IDs AND source-run evidence.
2. Verify octopus eligibility using original saved position calls and coverage.
   Record imputed one-isoform genes in the manifest, without inventing solutions.
3. Load one original protein matrix per worker; build compact AA compatibility
   lookups and query neighbors lazily. Reuse these lookups within the cohort.
4. Use the full matrix population at Fraction=1; otherwise map saved AvailableReads
   to their original proteins. Equal fractions need not mean equal populations.
5. Check S, construct C, solve exactly on C, and verify the extension witness.
   Prioritize canonical rows globally, then other full-fraction and smaller-fraction
   rows within each cohort. Commit each result and witness together to SQLite.
6. Export tables. COMPLETE means exact search and witness checks succeeded;
   INCOMPLETE means not attempted or interrupted, with a feasible lower bound when
   available; NOT_VERIFIABLE means exact sampled membership is missing;
   INVALID_INPUT means an input invariant failed, including a non-independent S.

NumIndividuallyAddable is |C|; MaxAdditionalTogether is |T| only when proven exact;
ExtendedIndependentSetSize is |S| + |T|. BestSavedSolutionSizeForSameGraph is an
observed heuristic size, not a global optimum. GapToBestSavedSolution subtracts
the original selected size, not the extended size. Unknown numbers remain blank.
SourceRow is the one-based CSV data-row number; OriginalSolutionIndex preserves
the source's Index value if present and is otherwise blank. GraphID hashes the
available-vertex mask in original matrix-row order. It is meaningful WITHIN the
same dataset/transcript and unchanged matrix; it is not a global graph identity
or a content hash of the amino-acid matrix.

Reading integer vertex sets: a four-protein example
---------------------------------------------------
Assign P0, P1, P2, P3 to matrix rows 0, 1, 2, 3. Binary displays put bit 0 on the
RIGHT: P0=0001, P1=0010, P2=0100, P3=1000. Thus {P0,P2}=0101 (integer 5).
Suppose the only edges are P0--P1 and P2--P3. N(P0)=0010 and N(P2)=1000.
For S={P0}=0001 and available=1111, blocked = S | N(P0) = 0011 and
candidates = available & ~blocked = 1100, representing {P2,P3}. Either candidate
can join S, but they conflict: |C|=2, maximum additions=1, extended size=2.

The operators are set operations: & is intersection, | is union, and x & ~y
removes y from x. Python's ~ has infinitely many leading one bits, so always
intersect it with a finite vertex set; do not count or iterate ~y on its own.
The shift 1 << i builds the singleton containing vertex i. No dense pairwise
distance matrix is needed. The bit helpers and graph/search functions below
explain how these same operations implement the biological definitions.
"""
import os

# Set before importing NumPy/pandas, including in spawned cohort workers.
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import argparse
import csv
import functools
import hashlib
import json
import multiprocessing
import resource
import signal
import sqlite3
import sys
import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS = ROOT / "Code/Notebooks"
OUTPUT = NOTEBOOKS / "JointMISAnalysis"
OCTOPUS = ROOT / "O.vulgaris/MpileupAndTranscripts/PRJNA791920"
PARAMETERS = ("Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition")
IDENTITY = ("Species", "Dataset", "Platform", "Gene", "TranscriptID",
            "CoverageThreshold", "SolutionFile", "ProteinMatrixFile")
RESULT_COLUMNS = ("SolutionID", *IDENTITY, "SourceRow", "OriginalSolutionIndex",
    *PARAMETERS, "CanonicalStatus", "CanonicalSource", "NumVerticesInAvailableGraph",
    "OriginalSelectedSize", "IsIndependent", "IsMaximal", "NumIndividuallyAddable",
    "MaxAdditionalTogether", "ExtendedIndependentSetSize", "BestSavedSolutionSizeForSameGraph",
    "GapToBestSavedSolution", "ResultStatus", "Notes", "GraphID",
    "AdditionalLowerBound", "ElapsedSeconds")
MANIFEST_COLUMNS = (*IDENTITY[:6], "EligibleForManuscript", *IDENTITY[6:],
    "HasOriginalSolverResult", "CanonicalIdentificationStatus", "ExclusionOrMissingReason",
    "AdditionalSolutionFile", "CanonicalSourceRow", "CanonicalSource", "EligibilityCheck")
WITNESS_COLUMNS = ("SolutionID", "Dataset", "Gene", "SourceRow", "ProteinID", "ExtensionRank")
TERMINAL = {"COMPLETE", "NOT_VERIFIABLE", "INVALID_INPUT"}
csv.field_size_limit(256 * 1024**2)


def read_tsv(path, **kwargs):
    """Return a pandas table from a path, preserving literal IDs such as 'NA'.

    Extra pandas reader arguments select columns, dtypes, or chunks. Missing-value
    inference is disabled so protein/read identifiers are not silently converted.
    """
    return pd.read_csv(path, sep="\t", keep_default_na=False, **kwargs)


def write_tsv(path, rows, columns):
    """Write iterable row dictionaries in the requested column order to a Path.

    Return None. A sibling temporary file replaces the destination atomically;
    callers must supply a generated-output path, never an original scientific input.
    """
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def stamp(path):
    """Return a path's byte size and nanosecond mtime for resume checks.

    These inexpensive metadata checks detect ordinary file changes; they are not
    content-level SHA-256 verification.
    """
    stat = Path(path).stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def check(condition, message):
    """Return None for a satisfied invariant; otherwise raise ValueError(message).

    Checks remain active even when Python is run with assertion optimization.
    """
    if not condition:
        raise ValueError(message)


def protein_ids(text):
    """Split a comma-delimited original ID string into a list without coercion.

    Also used for read IDs. An empty field becomes an empty list; order and any
    duplicates are preserved for the caller to validate.
    """
    return text.split(",") if text else []


def solution_key(cohort, path, row_number):
    """Return a stable solution ID from dataset, source path, and one-based row.

    cohort is a manifest/config dictionary; path is below ROOT. The row number
    excludes the header and is independent of any stored Index column.
    """
    relative = Path(path).relative_to(ROOT)
    return f"{cohort['Dataset']}:{relative}:row={row_number}"


def source_rows(cohort):
    """Read one cohort's original solution files into a list of row dictionaries.

    cohort supplies SolutionFile and optionally AdditionalSolutionFile. Each row
    retains source strings and gains its actual SolutionFile, SourceRow and
    SolutionID. Separate PacBio1 0.1 rows are omitted only if the main file already
    has identical run parameters, selected IDs AND AvailableReads. Different saved
    runs within a file remain separate, even when their selected sets coincide.
    No random subsampling or reinterpretation of the original Index is performed.
    """
    rows, main_signatures = [], set()
    for file_number, path in enumerate(filter(None, [cohort["SolutionFile"],
                                                    cohort.get("AdditionalSolutionFile")])):
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = set(PARAMETERS) | {"NumUniqueSamples", "UniqueSamples"}
            check(required <= set(reader.fieldnames), f"Wrong solution schema: {path}")
            for number, row in enumerate(reader, 1):
                if file_number:
                    check(float(row["Fraction"]) == 0.1, f"Non-0.1 row in {path}")
                # Only compare supplemental rows with the main file. Repeated
                # original heuristic runs remain distinct observations.
                signature = (tuple(row[k] for k in PARAMETERS),
                    tuple(sorted(protein_ids(row["UniqueSamples"]))),
                    tuple(sorted(protein_ids(row.get("AvailableReads", "")))))
                if file_number and signature in main_signatures:
                    continue
                if not file_number:
                    main_signatures.add(signature)
                row.update(SolutionFile=path, SourceRow=number,
                           SolutionID=solution_key(cohort, path, number))
                rows.append(row)
    return rows


def same_parameters(row, saved):
    """Return whether two mapping-like rows describe the same four run parameters.

    Algorithm is compared as text; numeric fields accept equivalent CSV forms
    such as '4' and '4.0'. Matching parameters alone does not identify a protein set.
    """
    return all(row[k] == str(saved[k]) if k == "Algorithm"
               else float(row[k]) == float(saved[k]) for k in PARAMETERS)


def identify_canonical(cohort, rows, bindings, pacbio3, saved_octopus):
    """Return (SolutionID or None, EXACT/UNKNOWN, evidence text) for a cohort.

    rows contains original source dictionaries. bindings maps (platform, gene)
    to ordinary squid figure bindings; pacbio3 is the saved WithUMIs export table;
    saved_octopus maps dataset/transcript to Figure 6 rows. Each branch requires
    exact protein-ID set identity together with its available run/source-row
    evidence. Equal cardinality is insufficient: tied sets can have different
    extensions. Exactly one matching source row is required; ambiguity stays UNKNOWN.
    Expression exports establish identity only; entropy flags do not filter IDs.
    """
    dataset, gene = cohort["Dataset"], cohort["Gene"]
    if dataset.startswith("Octopus"):
        saved = saved_octopus[dataset][cohort["TranscriptID"]]
        ids = set(protein_ids(saved["Proteins"]))
        matches = [r for r in rows if same_parameters(r, saved)
                   and set(protein_ids(r["UniqueSamples"])) == ids]
        source = (f"MaxDistinctProtsForFig6.TMR{cohort['CoverageThreshold']}.Octopus.Pooled.csv; "
                  "exact saved Proteins and all four run parameters; "
                  "notebook sort(Fraction descending)/groupby(Chrom)/nlargest(1) tie resolved by saved row")
    elif dataset == "PacBio3":
        saved = pacbio3.loc[pacbio3.Chrom.eq(cohort["TranscriptID"])]
        ids = set(saved.Protein)
        indices = set(saved["#Solution"])
        matches = [r for r in rows if len(indices) == 1
                   and str(r["SourceRow"]) in indices and float(r["Fraction"]) == 1
                   and set(protein_ids(r["UniqueSamples"])) == ids]
        source = ("AssignedExpression.PacBio3.WithUMIs.tsv: exact Protein set and #Solution; "
                  "MergedSamples notebook cells d3f36c3a/8fb523c5/2df283d5; "
                  "Figure4b source cells 53b4f680/eb4fe8e2/8429762c; no entropy filtering")
    else:
        saved = bindings[(cohort["Platform"], gene)]
        ids = set(saved["ordered_ids"])
        matches = [r for r in rows if same_parameters(r, saved["parameters"])
                   and str(r["SourceRow"]) == str(saved["solution_id"])
                   and set(protein_ids(r["UniqueSamples"])) == ids]
        source = ("figure_rebuild_work/canonical_inputs/pattern_bindings.json: "
                  "ordinary_full_fraction, exact ordered_ids set, solution_id and run parameters")
    if len(matches) != 1:
        return None, "UNKNOWN", f"{source}; matching original rows={len(matches)}"
    return matches[0]["SolutionID"], "EXACT", source


def connect():
    """Open and return the writable SQLite result checkpoint with a busy timeout.

    Cohort workers share this database; each solution and witness are committed
    together. Calling this helper is not a read-only inspection operation.
    """
    connection = sqlite3.connect(str(OUTPUT / "checkpoint.sqlite"), timeout=120)
    connection.execute("PRAGMA busy_timeout=120000")
    return connection


def build_inventory():
    """Build/refresh the configured inventory and return (cohorts, provenance notes).

    No arguments: paths come from ROOT/OUTPUT. input_config.json explicitly fixes
    2 PacBio1, 3 March-2026 PacBio3 and 19 Illumina gene/dataset pairs from their
    notebooks. Figure 6 Chrom columns fix 1,711 TMR50 and 86 TMR1000 octopus genes;
    each branch uses its own original solution and protein directory. Imputed
    one-isoform entries have no solver rows. Updated/entropy/distance variants are
    excluded, and duplicate or ambiguous source selections fail explicitly.

    Canonical identity is established before graph work. OriginalSolutionIndex is
    copied only from an actual Index field. The function writes inventory/source
    metadata and inserts missing checkpoint rows without replacing existing results.
    It checks recorded source signatures/stamps on resume. It is intentionally NOT
    called by the documentation-only targeted verification procedure.
    """
    config = json.loads((OUTPUT / "input_config.json").read_text())
    for path, expected in config["source_code_sha256"].items():
        check(hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected,
              f"Provenance source changed since input configuration: {path}")
    bindings_path = NOTEBOOKS / "figure_rebuild_work/canonical_inputs/pattern_bindings.json"
    bindings = {(b["platform"], b["gene"]): b for b in json.loads(bindings_path.read_text())["cohorts"]
                if b["cohort"] == "ordinary_full_fraction"}
    pb3_path = NOTEBOOKS / "AssignedExpression.PacBio3.WithUMIs.tsv"
    pacbio3 = read_tsv(pb3_path, usecols=["Chrom", "Protein", "#Solution"], dtype=str)
    cohorts, saved_octopus = config["cohorts"], {}
    sources = {str(p): stamp(p) for p in [bindings_path, pb3_path, *config["source_code_sha256"]]}
    for threshold, expected_size in [(50, 1711), (1000, 86)]:
        dataset = f"Octopus_TMR{threshold}"
        table_path = NOTEBOOKS / f"MaxDistinctProtsForFig6.TMR{threshold}.Octopus.Pooled.csv"
        saved = read_tsv(table_path, dtype=str)
        check(saved.Chrom.is_unique and len(saved) == expected_size,
              f"Unexpected saved manuscript population: {dataset}, {len(saved)}")
        saved_octopus[dataset] = saved.set_index("Chrom").to_dict("index")
        sources[str(table_path)] = stamp(table_path)
        base = OCTOPUS / f"IsoSeq.Polished.Unclustered.TotalCoverage{threshold}.PooledSamples"
        for record in saved.to_dict("records"):
            chrom = record["Chrom"]
            paths = [p for p in (base / "DistinctProteins").glob(f"{chrom}.DistinctUniqueProteins.*.csv")
                     if not any(t in p.name.lower() for t in
                                ["updated", "entropy", "expression", "miyata", "grantham", "aagroups"])]
            check(len(paths) <= 1, f"Ambiguous original files: {dataset}/{chrom}: {paths}")
            if not paths:
                check(float(record["NumOfProteins"]) == 1 and not record["Proteins"],
                      f"Missing solver file for a non-imputed Figure 6 row: {dataset}/{chrom}")
            cohorts.append(dict(Species="Octopus", Dataset=dataset, Platform="PacBio",
                Gene=record["Transcript"], TranscriptID=chrom, CoverageThreshold=threshold,
                SolutionFile=str(paths[0]) if paths else "",
                ProteinMatrixFile=str(base / "ProteinsFiles" / f"{chrom}.unique_proteins.csv.gz")))
    database = connect()
    database.execute("PRAGMA journal_mode=WAL")
    database.execute("CREATE TABLE IF NOT EXISTS solutions (id TEXT PRIMARY KEY, cohort TEXT, "
                     "canonical INTEGER, metadata TEXT, result TEXT, witness TEXT)")
    for cohort in cohorts:
        cohort.setdefault("CoverageThreshold", "")
        cohort.update(EligibleForManuscript=True, EligibilityCheck="PENDING" if cohort["Species"] == "Octopus" else "EXPLICIT_INVENTORY")
        cohort["HasOriginalSolverResult"] = bool(cohort["SolutionFile"])
        cohort["CohortID"] = cohort["Dataset"] + ":" + cohort["TranscriptID"]
        if not cohort["SolutionFile"]:
            cohort.update(CanonicalIdentificationStatus="NOT_APPLICABLE",
                          ExclusionOrMissingReason="IMPUTED_ONE_ISOFORM / NO_SOLVER_RESULT")
            continue
        for path in filter(None, [cohort["SolutionFile"], cohort["ProteinMatrixFile"], cohort.get("AdditionalSolutionFile")]):
            check(Path(path).is_file(), f"Missing original input: {path}")
            sources[path] = stamp(path)
        rows = source_rows(cohort)
        canonical, status, evidence = identify_canonical(cohort, rows, bindings, pacbio3, saved_octopus)
        cohort.update(CanonicalID=canonical, CanonicalIdentificationStatus=status,
                      CanonicalSource=evidence, ExclusionOrMissingReason="" if canonical else evidence,
                      CanonicalSourceRow=next((r["SourceRow"] for r in rows if r["SolutionID"] == canonical), ""))
        full_best = max(int(r["NumUniqueSamples"]) for r in rows if float(r["Fraction"]) == 1)
        for row in rows:
            metadata = {k: cohort[k] for k in IDENTITY}
            metadata.update({k: row[k] for k in (*PARAMETERS, "SolutionFile", "SourceRow", "SolutionID")})
            metadata.update(OriginalSolutionIndex=row.get("Index", ""),
                OriginalSelectedSize=int(row["NumUniqueSamples"]),
                CanonicalStatus="EXACT" if row["SolutionID"] == canonical else ("UNKNOWN" if not canonical else "NOT_APPLICABLE"),
                CanonicalSource=evidence if row["SolutionID"] == canonical or not canonical else "")
            result = dict(ResultStatus="INCOMPLETE", Notes="Not yet attempted")
            if float(row["Fraction"]) == 1:
                result.update(BestSavedSolutionSizeForSameGraph=full_best,
                              GapToBestSavedSolution=full_best - metadata["OriginalSelectedSize"])
            elif not row.get("AvailableReads"):
                result.update(ResultStatus="NOT_VERIFIABLE",
                    Notes="Original solution has no AvailableReads; no exact saved read/vertex subset is identified. Full graph is not substituted.")
            database.execute("INSERT OR IGNORE INTO solutions VALUES (?,?,?,?,?,?)",
                (row["SolutionID"], cohort["CohortID"], int(row["SolutionID"] == canonical),
                 json.dumps(metadata), json.dumps(result), "[]"))
    database.commit()
    database.close()
    actual_paths = [c["SolutionFile"] for c in cohorts if c["SolutionFile"]]
    check(len(actual_paths) == len(set(actual_paths)), "Duplicate selected source file")
    sources_file = OUTPUT / "input_sources.json"
    if sources_file.exists():
        old_sources = json.loads(sources_file.read_text())
        for path, old in old_sources.items():
            check(stamp(path) == old, f"Original input changed; cannot reuse checkpoint: {path}")
        sources.update(old_sources)
    if (OUTPUT / "eligibility_check.tsv").exists():
        eligibility = read_tsv(OUTPUT / "eligibility_check.tsv")
        membership = {(r.Dataset, r.Chrom): r.ReconstructedEligible
                      for r in eligibility.itertuples(index=False)}
        for cohort in cohorts:
            if cohort["Species"] == "Octopus":
                cohort["EligibilityCheck"] = "MATCH" if membership.get(
                    (cohort["Dataset"], cohort["TranscriptID"]), False) else "DISCREPANCY"
    sources_file.write_text(json.dumps(sources, indent=2) + "\n")
    (OUTPUT / "cohorts.json").write_text(json.dumps(cohorts, indent=2) + "\n")
    write_tsv(OUTPUT / "input_manifest.tsv", cohorts, MANIFEST_COLUMNS)
    return cohorts, config["discrepancies"]


def verify_octopus_eligibility(cohorts, deadline):
    """Write the reconstructed eligibility comparison, updating cohort dictionaries.

    cohorts is the inventory list; deadline is an absolute time.time() timestamp.
    Return None, reusing an existing eligibility_check.tsv when present. Otherwise
    stream original TMR50 position files, retain NoisyFinal calls, and count positions
    with MismatchFrequency >=0.05. Use the original alternative-base count rule,
    including EditedFinal, without rerunning editing/noise detection. Alternative
    bases tied for the maximum have the same count, so the frequency needs no draw.

    Both notebook branches join this SAME SNP table to alignment coverage, fill
    uncounted SNP totals with zero, and require SNPs<=3 plus MappedReads>=50 or
    >=1000. Compare with saved manuscript sets; do not replace them on disagreement.
    This writes checkpoints and reads substantial position data when uncached;
    targeted verification instead inspects the existing comparison table.
    """
    output = OUTPUT / "eligibility_check.tsv"
    if output.exists():
        return  # Original input stamps are verified by build_inventory on resume.
    alignment_path = ROOT / "O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv"
    alignment = read_tsv(alignment_path, usecols=["Chrom", "MappedReads"])
    check(alignment.Chrom.is_unique, "Alignment Chrom key is not unique")
    sources = json.loads((OUTPUT / "input_sources.json").read_text())
    sources[str(alignment_path)] = stamp(alignment_path)
    positions_dir = OCTOPUS / "IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples/PositionsFiles"
    counts_path = OUTPUT / "snp_counts_checkpoint.json"
    counts = json.loads(counts_path.read_text()) if counts_path.exists() else {}
    columns = ["Chrom", "RefBase", "A", "T", "C", "G", "NoisyFinal", "EditedFinal"]
    # The notebook right-joins alignment genes with coverage >=50 before reading positions.
    for chrom in alignment.loc[alignment.MappedReads.ge(50), "Chrom"]:
        path = positions_dir / f"{chrom}.positions.csv.gz"
        if chrom in counts:
            continue
        if time.time() >= deadline:
            raise TimeoutError("Total run budget during SNP eligibility reconstruction")
        if not path.exists():
            counts[chrom] = 0  # Same missing-position/outer-join/fillna(0) behavior.
            continue
        sources[str(path)] = stamp(path)
        count = 0
        for chunk in read_tsv(path, usecols=columns, chunksize=20000):
            check(chunk.Chrom.eq(chrom).all(), f"Positions Chrom mismatch: {path}")
            check(chunk.NoisyFinal.isin([True, False]).all(), f"Invalid NoisyFinal: {path}")
            noisy = chunk.loc[chunk.NoisyFinal.eq(True)]
            for ref, group in noisy.groupby("RefBase"):
                check(ref in "ATCG", f"Unexpected reference base {ref}")
                alternative = group[[b for b in "ATCG" if b != ref]].max(axis=1)
                edited = group.EditedFinal.eq(True)
                if edited.any():
                    check(ref == "A", "EditedFinal row must have reference A, as in the notebook")
                    alternative = alternative.where(~edited, group.G)
                frequency = alternative / (group[ref] + alternative)
                count += int((alternative.gt(0) & frequency.ge(0.05)).sum())
        counts[chrom] = count
        if len(counts) % 100 == 0:
            counts_path.write_text(json.dumps(counts))
            (OUTPUT / "input_sources.json").write_text(json.dumps(sources, indent=2))
            print(f"SNP eligibility: {len(counts)} alignment genes", flush=True)
    counts_path.write_text(json.dumps(counts))
    (OUTPUT / "input_sources.json").write_text(json.dumps(sources, indent=2))
    alignment["SNPs"] = alignment.Chrom.map(counts).fillna(0).astype(int)
    records = []
    for threshold in (50, 1000):
        manuscript = {c["TranscriptID"] for c in cohorts if c["Dataset"] == f"Octopus_TMR{threshold}"}
        reconstructed = set(alignment.loc[alignment.MappedReads.ge(threshold) & alignment.SNPs.le(3), "Chrom"])
        lookup = alignment.set_index("Chrom").to_dict("index")
        for chrom in sorted(manuscript | reconstructed):
            records.append(dict(Dataset=f"Octopus_TMR{threshold}", Chrom=chrom,
                InManuscript=chrom in manuscript, ReconstructedEligible=chrom in reconstructed,
                **lookup[chrom]))
        for cohort in cohorts:
            if cohort["Dataset"] == f"Octopus_TMR{threshold}":
                cohort["EligibilityCheck"] = "MATCH" if cohort["TranscriptID"] in reconstructed else "DISCREPANCY"
        print(f"TMR{threshold} eligibility: saved={len(manuscript)}, reconstructed={len(reconstructed)}, "
              f"symmetric difference={len(manuscript ^ reconstructed)}", flush=True)
    write_tsv(output, records, ["Dataset", "Chrom", "InManuscript", "ReconstructedEligible", "MappedReads", "SNPs"])
    write_tsv(OUTPUT / "input_manifest.tsv", cohorts, MANIFEST_COLUMNS)
    (OUTPUT / "cohorts.json").write_text(json.dumps(cohorts, indent=2))


def bit_indices(mask):
    """Yield the zero-based protein-row indices present in a nonnegative integer.

    For mask=0101 (binary), yield 0 then 2. The input must be a finite vertex set,
    never the negative integer produced by an unbounded bitwise complement.
    """
    while mask:
        # Isolate the lowest occupied bit: 0101 & -0101 leaves 0001.
        # In binary arithmetic, negation flips the bits and adds one, so all
        # higher bits cancel under AND while the lowest set bit remains.
        bit = mask & -mask
        yield bit.bit_length() - 1
        # XOR toggles this known-present bit off, removing the yielded vertex.
        mask ^= bit


def bit_count(mask):
    """Return the cardinality of a finite, nonnegative vertex mask.

    Each occupied bit is one protein, so the number of binary '1' digits is the
    set size. Use this Python-3.8-compatible form instead of int.bit_count().
    """
    return bin(mask).count("1")  # Python 3.8 compatibility.


def vertex_mask(indices, size):
    """Return the integer set of row indices in a matrix with size vertices.

    indices is an iterable of valid integers in [0,size); repeated indices collapse
    to one vertex. Little-endian bit packing makes row i correspond to 1 << i.
    The temporary vector is length N, not an N-by-N pairwise matrix.
    """
    values = np.zeros(size, dtype=np.uint8)
    values[list(indices)] = 1
    return int.from_bytes(np.packbits(values, bitorder="little").tobytes(), "little")


def build_graph(path, need_reads):
    """Build lazy indistinguishability queries from an original protein matrix.

    path names the original TSV, possibly gzipped; need_reads requests an exact
    read-ID-to-protein mapping for saved subsamples. Return a six-item tuple:
    (ids, index, read_index, full, neighbors, direct_masks). ids maps row numbers
    back to original Protein IDs; index is its inverse; read_index is empty when
    not requested. full contains every row. neighbors(i) returns an integer mask
    of vertices indistinguishable from row i, excluding i itself. direct_masks
    holds per-column uint32 AA-set encodings for a separate witness check.

    Two different kinds of bits are used: vertex masks have one bit per protein;
    direct_masks have one bit per amino-acid symbol, including '*'. Protein IDs
    must be unique, AA cells must contain valid nonempty possibility sets, and
    each saved read must map to one protein. No entropy filter or distance graph
    is introduced. Only a bounded cache of queried neighbors is retained.
    """
    header = read_tsv(path, nrows=0).columns.tolist()
    # Original codon-position headers end in a parenthesized reference AA.
    # Cells such as 'K,R' are possibility SETS, not an ordered amino-acid sequence.
    aa_columns = [c for c in header if "(" in c and c.endswith(")")]
    check(aa_columns and header.index(aa_columns[0]) in (14, 15),
          f"Unexpected original AA schema: {path}")
    columns = ["Protein", *aa_columns] + (["Reads", "NumOfReads"] if need_reads else [])
    proteins = read_tsv(path, usecols=columns, dtype=str)
    ids = proteins.Protein.tolist()
    check(len(ids) == len(set(ids)), f"Duplicate protein IDs: {path}")
    size = len(ids)
    # N low bits set to one: for four rows, (1 << 4) - 1 == 1111.
    full = (1 << size) - 1
    index = {name: i for i, name in enumerate(ids)}
    read_index = {}
    if need_reads:
        for i, (reads, expected) in enumerate(zip(proteins.Reads, proteins.NumOfReads)):
            names = protein_ids(reads)
            check(len(names) == int(expected), f"Read count mismatch for protein {ids[i]}")
            for name in names:
                check(name not in read_index, f"Read belongs to multiple original proteins: {name}")
                read_index[name] = i
    lookups, direct_masks = [], []
    alphabet = "ACDEFGHIKLMNPQRSTVWY*"
    for col in aa_columns:
        # At one position, categorical code k labels one observed possibility
        # string. It is a lookup index, not an amino acid or a numerical distance.
        categories = pd.Categorical(proteins[col])
        codes = categories.codes.copy()
        possibilities = [set(value.split(",")) for value in categories.categories]
        check(all(value and value <= set(alphabet) for value in possibilities),
              f"Invalid amino-acid possibility set in {col}")
        # masks[k] is the set of protein rows having category k at this position.
        # These category masks partition the rows: each protein has exactly one.
        masks = [int.from_bytes(np.packbits(codes == k, bitorder="little").tobytes(), "little")
                 for k in range(len(possibilities))]
        # compatible[k] unions all categories sharing at least one AA with k.
        # Integer sum equals bitwise OR here because category masks are disjoint:
        # no occupied bit is added twice, so there are no binary carries.
        compatible = [sum(mask for mask, other in zip(masks, possibilities) if value & other)
                      for value in possibilities]
        # Smaller compatible populations tend to eliminate rows sooner. Sorting
        # positions by this score changes only evaluation order: set intersection
        # is commutative, so requiring overlap at ALL positions is unchanged.
        score = sum(bit_count(mask) for mask in compatible) / len(compatible)
        lookups.append((score, codes, compatible))
        # Separate AA-symbol masks support direct comparisons of witness patterns.
        aa_bits = np.array([sum(1 << alphabet.index(a) for a in value) for value in possibilities], dtype=np.uint32)
        direct_masks.append(aa_bits[codes])
    lookups.sort(key=lambda item: item[0])
    # Bound the neighbor cache to both 1024 vertices and roughly 64 MiB.
    cache_size = max(1, min(1024, (64 * 1024**2) // max(1, (size + 7) // 8)))

    @functools.lru_cache(maxsize=cache_size)
    def neighbors(vertex):
        """Return the finite neighbor mask of one zero-based protein-row index.

        An edge requires AA overlap at every position; a single disjoint position
        removes that other protein permanently. Self-overlap is not a graph edge.
        """
        remaining = full
        # Left shift builds the singleton {vertex}; e.g. 1 << 2 == 0100.
        self_bit = 1 << vertex
        for _, codes, compatible in lookups:
            # Intersect possible neighbors with proteins compatible at THIS
            # position. Survivors after ALL positions are indistinguishable.
            remaining &= compatible[codes[vertex]]
            if remaining == self_bit:
                # Only self remains; later intersections cannot add neighbors.
                return 0
        # Set difference removes self, keeping only other indistinguishable rows.
        return remaining & ~self_bit

    return ids, index, read_index, full, neighbors, direct_masks


def exact_extension(candidates, neighbors, deadline):
    """Return (added row indices, proven_exact, interruption_reason) for C.

    candidates is a finite integer mask whose vertices are already compatible
    with the fixed original S. neighbors(i) supplies conflict edges; deadline is
    an absolute time.time() timestamp. This is an UNWEIGHTED cardinality problem.

    Each stack entry holds remaining (candidate vertices still available) and
    chosen (an independent list of additions already taken on that branch).
    Remaining vertices have no edge to chosen. best is the largest feasible list
    found anywhere so far. Even taking every remaining vertex cannot exceed
    len(chosen)+|remaining|, giving a safe upper bound for pruning.

    Isolated vertices can be included without losing anything. For a degree-one
    vertex v with sole neighbor u, take a maximum independent set: if it contains
    u, replace u with v, which conflicts with no other remaining vertex. If it
    contains neither, v could be added, contradicting maximal cardinality. Thus
    some optimum contains v, so including v and deleting u is exact. This argument
    does not apply to arbitrary weighted objectives.

    Otherwise branch on a vertex: exclude it, or include it and remove all its
    neighbors. Every independent set follows one branch. Safe reductions and
    upper-bound pruning preserve an optimum; exhausting the stack proves the
    returned size is maximum within C. On time/memory interruption, return the
    best feasible set with proven_exact=False; its size is only a lower bound.
    """
    best = []
    stack = [(candidates, [])]
    try:
        while stack:
            if time.time() >= deadline:
                raise TimeoutError("Exact candidate search time limit")
            remaining, chosen = stack.pop()
            # Upper bound: even an edge-free remainder cannot improve best.
            if len(chosen) + bit_count(remaining) <= len(best):
                continue
            while remaining:
                pivot, highest_degree, reduced = None, -1, False
                for vertex in bit_indices(remaining):
                    # Degree is measured only in the current candidate subgraph.
                    adjacent = neighbors(vertex) & remaining
                    degree = bit_count(adjacent)
                    if degree <= 1:
                        chosen = chosen + [vertex]
                        # Include vertex; subtract both it and its current
                        # neighbors. Everything left stays compatible with chosen.
                        remaining &= ~(adjacent | (1 << vertex))
                        reduced = True
                        break
                    if degree > highest_degree:
                        pivot, highest_degree = vertex, degree
                if len(chosen) > len(best):
                    best = chosen
                if not reduced:
                    break
            if not remaining or len(chosen) + bit_count(remaining) <= len(best):
                continue
            without = remaining & ~(1 << pivot)
            # Exclude pivot. The other branch below includes pivot and subtracts
            # its neighbors: without & ~neighbors(pivot) is a set difference.
            # The stack is LIFO, so the include branch is explored first.
            stack.append((without, chosen))
            stack.append((without & ~neighbors(pivot), chosen + [pivot]))
        return best, True, ""
    except (TimeoutError, MemoryError) as error:
        return best, False, f"{type(error).__name__}: {error}"


def verify_witness(added, selected, direct_masks):
    """Check an extension directly from AA-symbol masks; return None or raise.

    added and selected are lists of original matrix-row indices. direct_masks is
    the list of per-position uint32 arrays returned by build_graph. For each added
    row, compare with all of S and earlier additions, retaining a Boolean flag
    only while the AA sets overlap at every position. Any final True flag is a
    forbidden indistinguishable pair. This does not query cached vertex neighbors;
    it verifies feasibility, not the optimality established by exact_extension.
    """
    previous = list(selected)
    for vertex in added:
        overlap = np.ones(len(previous), dtype=bool)
        for masks in direct_masks:
            # Here & intersects AA-symbol sets, not sets of protein vertices.
            overlap &= (masks[previous] & masks[vertex]) != 0
        check(not overlap.any(), "Extension witness contains indistinguishable proteins")
        previous.append(vertex)


def analyze_solution(row, graph, available, best_saved, seconds, global_deadline):
    """Check one saved S and return (result dictionary, added original protein IDs).

    row is an original solution dictionary; graph is build_graph's tuple;
    available is that row's exact vertex mask; best_saved is the largest original
    saved size for THIS population in THIS cohort. seconds is a relative timeout
    and global_deadline an absolute timestamp. Selected IDs must be unique, match
    the saved count, exist in the matrix and lie in available. An original conflict
    returns INVALID_INPUT. Otherwise construct C, search it exactly and check the
    witness. COMPLETE requires both exactness and successful verification.

    GraphID hashes little-endian membership bytes in this fixed matrix order;
    use it with dataset/transcript/matrix provenance, never across unrelated graphs.
    Interrupted results keep unknown exact quantities absent and a feasible lower
    bound when available. The caller, not this function, persists the result.
    """
    ids, index, _, _, neighbors, direct_masks = graph
    started = time.time()
    result = dict(ResultStatus="INCOMPLETE", Notes="", AdditionalLowerBound=None)
    added = []
    with time_limit(min(seconds, max(0.01, global_deadline - time.time()))):
        try:
            names = protein_ids(row["UniqueSamples"])
            check(len(names) == len(set(names)) == int(row["NumUniqueSamples"]), "Selected protein count/uniqueness mismatch")
            selected = [index[name] for name in names]
            selected_mask = vertex_mask(selected, len(ids))
            check(not selected_mask & ~available, "Selected protein outside the exact available graph")
            result.update(NumVerticesInAvailableGraph=bit_count(available),
                          GraphID=hashlib.sha256(available.to_bytes((len(ids) + 7) // 8, "little")).hexdigest(),
                          BestSavedSolutionSizeForSameGraph=best_saved,
                          GapToBestSavedSolution=best_saved - len(selected))
            # blocked = S union N(S): selected vertices are unavailable for addition,
            # and any neighbor would conflict with at least one fixed member of S.
            blocked = selected_mask
            for vertex in selected:
                adjacent = neighbors(vertex)
                if adjacent & selected_mask:
                    other = next(bit_indices(adjacent & selected_mask))
                    result.update(IsIndependent=False, ResultStatus="INVALID_INPUT",
                                  Notes=f"Original selected conflict: {ids[vertex]},{ids[other]}")
                    return result, []
                # blocked |= neighbors(vertex): union this member's forbidden
                # neighbors into the running set; OR never removes prior entries.
                blocked |= adjacent
            # C = available minus blocked. The finite available mask also removes
            # Python's infinite leading ones from the complement ~blocked.
            candidates = available & ~blocked
            result.update(IsIndependent=True, IsMaximal=not candidates,
                          NumIndividuallyAddable=bit_count(candidates), AdditionalLowerBound=0)
            # Reserve a little time to verify and persist a feasible witness after timeout.
            added, exact, reason = exact_extension(candidates, neighbors,
                min(started + max(0.01, seconds - 2), global_deadline - 2))
            verify_witness(added, selected, direct_masks)
            result["AdditionalLowerBound"] = len(added)
            if exact:
                result.update(MaxAdditionalTogether=len(added),
                    ExtendedIndependentSetSize=len(selected) + len(added), ResultStatus="COMPLETE")
            else:
                result["Notes"] = reason
        except (TimeoutError, MemoryError) as error:
            result["Notes"] = f"{type(error).__name__}: {error}"
            # The exact search constructs a feasible set; an interrupted
            # independent witness check still cannot be called COMPLETE.
            result["AdditionalLowerBound"] = len(added) if result.get("IsIndependent") else None
        except (ValueError, KeyError) as error:
            result.update(ResultStatus="INVALID_INPUT", Notes=f"{type(error).__name__}: {error}")
    result["ElapsedSeconds"] = round(time.time() - started, 4)
    return result, [ids[i] for i in added]


@contextmanager
def time_limit(seconds):
    """Yield inside a seconds-long wall-clock alarm, raising TimeoutError on expiry.

    Linux/main-thread context manager used for a single solution or preparation
    block. Cancel the alarm and restore the prior signal handler on every exit.
    """
    def expired(signum, frame):
        """Translate the SIGALRM callback into the exception handled by the search."""
        raise TimeoutError("Per-solution or preparation time limit")
    previous = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def worker(cohort, phase, limits, deadline):
    """Process one cohort/phase and return (CohortID, completed count, message).

    cohort is an inventory dictionary; phase is 'canonical' or 'remaining'; limits
    gives memory_gib, canonical_seconds and other_seconds; deadline is absolute
    wall-clock time. Load only this cohort's matrix and reuse it across its target
    rows. A separate global canonical pass may require reloading the matrix later.
    COMPLETE, NOT_VERIFIABLE and INVALID_INPUT rows are skipped on resume; unfinished
    rows can be retried. Remaining targets run in descending Fraction order.

    Full-fraction rows use all original vertices. Lower fractions use the actual
    saved read IDs, mapped to proteins with at least one sampled supporting read.
    best_by_population groups by exact membership masks, not fraction labels or
    repetition numbers. Each result/witness pair is committed immediately to SQLite.
    Time/memory limits produce INCOMPLETE records, never approximate COMPLETE ones.
    """
    cap = int(limits["memory_gib"] * 1024**3)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (min(cap, hard) if hard >= 0 else cap, hard))
    database = connect()
    saved = {sid: json.loads(result) for sid, result in database.execute(
        "SELECT id,result FROM solutions WHERE cohort=?", (cohort["CohortID"],))}
    rows = source_rows(cohort)
    targets = [r for r in rows if (r["SolutionID"] == cohort.get("CanonicalID")) == (phase == "canonical")
               and saved[r["SolutionID"]]["ResultStatus"] not in TERMINAL]
    targets.sort(key=lambda r: (-float(r["Fraction"]), r["SolutionFile"], r["SourceRow"]))
    if not targets:
        database.close()
        return cohort["CohortID"], 0, "Already completed or not verifiable"
    seconds = limits["canonical_seconds"] if phase == "canonical" else limits["other_seconds"]
    try:
        with time_limit(min(seconds, max(0.01, deadline - time.time()))):
            needs_reads = any(float(r["Fraction"]) < 1 and r.get("AvailableReads") for r in rows)
            graph = build_graph(cohort["ProteinMatrixFile"], needs_reads)
            ids, index, reads, full, neighbors, direct_masks = graph
            populations, best_by_population, population_errors = {}, {}, {}
            # Exact masks, not nominal fractions, define which saved sizes are comparable.
            for row in rows:
                sid = row["SolutionID"]
                if float(row["Fraction"]) == 1:
                    available = full
                elif row.get("AvailableReads"):
                    try:
                        available = vertex_mask({reads[name] for name in protein_ids(row["AvailableReads"])}, len(ids))
                    except KeyError as error:
                        population_errors[sid] = f"Saved AvailableReads contains an unmapped original read: {error}"
                        continue
                else:
                    continue
                populations[sid] = available
                best_by_population[available] = max(best_by_population.get(available, 0), int(row["NumUniqueSamples"]))
        completed = 0
        for row in targets:
            if time.time() >= deadline:
                break
            sid = row["SolutionID"]
            if sid in population_errors:
                result, witness = dict(ResultStatus="INVALID_INPUT", Notes=population_errors[sid]), []
            else:
                available = populations[sid]
                result, witness = analyze_solution(row, graph, available,
                    best_by_population[available], seconds, deadline)
            database.execute("UPDATE solutions SET result=?,witness=? WHERE id=?",
                             (json.dumps(result), json.dumps(witness), sid))
            database.commit()  # Each solution and its witness are one transaction.
            completed += result["ResultStatus"] == "COMPLETE"
        return cohort["CohortID"], completed, "Processed"
    except (MemoryError, TimeoutError, ValueError, KeyError) as error:
        status = "INCOMPLETE" if isinstance(error, (MemoryError, TimeoutError)) else "INVALID_INPUT"
        for row in targets:
            current = json.loads(database.execute("SELECT result FROM solutions WHERE id=?",
                                 (row["SolutionID"],)).fetchone()[0])
            if current["ResultStatus"] not in TERMINAL:
                result = dict(ResultStatus=status, Notes=f"Cohort preparation: {type(error).__name__}: {error}")
                database.execute("UPDATE solutions SET result=? WHERE id=?", (json.dumps(result), row["SolutionID"]))
        database.commit()
        return cohort["CohortID"], 0, f"{status}: {error}"
    finally:
        database.close()


def export_canonical():
    """Refresh canonical_solutions.tsv from committed exact-identity rows; return None.

    'canonical=1' denotes verified source-row identity, not computation completion.
    Consequently the table also preserves identified canonicals still INCOMPLETE.
    This lightweight export keeps canonical results visible during other work.
    """
    database = connect()
    rows = [{**json.loads(metadata), **json.loads(result)} for metadata, result in
            database.execute("SELECT metadata,result FROM solutions WHERE canonical=1 ORDER BY rowid")]
    database.close()
    write_tsv(OUTPUT / "canonical_solutions.tsv", rows, RESULT_COLUMNS)


def export_results(cohorts, notes, started, command, limits, final=False):
    """Export committed results and return their combined metadata/result dictionaries.

    cohorts and notes describe inventory/provenance; started is invocation start
    time; command and limits document execution. final additionally checks input
    size/mtime stamps. Write all_solutions, canonical_solutions, extension_witnesses,
    incomplete_extension_witnesses and run_summary under OUTPUT. Unknown fields
    remain blank. Exact witnesses and incomplete lower-bound witnesses are separated.
    Existing generated tables are replaced; original input files are never written.

    Canonical identity and computational status are distinct. GraphID is scoped
    to its original cohort/matrix. The summary reports observed saved sizes and
    fixed-set extensions without interpreting either as a global MIS optimum.
    """
    database = connect()
    records, witnesses, lower_witnesses = [], [], []
    for sid, metadata, result, witness in database.execute("SELECT id,metadata,result,witness FROM solutions ORDER BY rowid"):
        record = {**json.loads(metadata), **json.loads(result)}
        records.append(record)
        for rank, protein in enumerate(json.loads(witness), 1):
            target = witnesses if record["ResultStatus"] == "COMPLETE" else lower_witnesses
            target.append({**record, "ProteinID": protein, "ExtensionRank": rank})
    database.close()
    canonical = [r for r in records if r["CanonicalStatus"] == "EXACT"]
    write_tsv(OUTPUT / "all_solutions.tsv", records, RESULT_COLUMNS)
    write_tsv(OUTPUT / "canonical_solutions.tsv", canonical, RESULT_COLUMNS)
    write_tsv(OUTPUT / "extension_witnesses.tsv", witnesses, WITNESS_COLUMNS)
    write_tsv(OUTPUT / "incomplete_extension_witnesses.tsv", lower_witnesses, WITNESS_COLUMNS)
    summary = ["JointMISAnalysis: maximum extension of each fixed original set; no global MIS claim.",
        f"Script: {Path(__file__).resolve()}", f"Command: {command}",
        f"Python: {sys.executable}; Python {sys.version.split()[0]}; NumPy {np.__version__}; pandas {pd.__version__}",
        f"Elapsed this invocation: {time.time() - started:.1f} seconds; limits: {json.dumps(limits)}",
        "Checkpoint: checkpoint.sqlite; each solution commits independently; COMPLETE rows are skipped on resume.",
        "GapToBestSavedSolution = best saved size for exactly the same vertex population minus original selected size.",
        "UNKNOWN canonical identity is not resolved by choosing an arbitrary maximum/tie.",
        "Input paths: input_manifest.tsv and input_sources.json (including eligibility inputs).",
        "Discrepancies and provenance notes:", *notes, "", "Dataset coverage:"]
    for dataset in dict.fromkeys(c["Dataset"] for c in cohorts):
        subset = [r for r in records if r["Dataset"] == dataset]
        cs = [r for r in subset if r["CanonicalStatus"] == "EXACT"]
        counts = {s: sum(r["ResultStatus"] == s for r in subset) for s in ["COMPLETE", "INCOMPLETE", "NOT_VERIFIABLE", "INVALID_INPUT"]}
        additions = [r["MaxAdditionalTogether"] for r in subset if r.get("MaxAdditionalTogether") is not None]
        summary.append(f"{dataset}: cohorts={sum(c['Dataset'] == dataset for c in cohorts)}, "
            f"solver cohorts={sum(c['Dataset'] == dataset and c['HasOriginalSolverResult'] for c in cohorts)}, "
            f"solutions={len(subset)}, statuses={counts}, verified canonicals={len(cs)}, "
            f"completed canonicals={sum(r['ResultStatus'] == 'COMPLETE' for r in cs)}, "
            f"nonmaximal canonicals={sum(r.get('IsMaximal') is False for r in cs)}, "
            f"nonmaximal known noncanonicals={sum(r.get('IsMaximal') is False and r['CanonicalStatus'] == 'NOT_APPLICABLE' for r in subset)}, "
            f"maximum exact additions={max(additions) if additions else 'unknown'}")
    summary.extend(["", "Canonical sizes (all identified rows; blank values mean unfinished):",
                    "Dataset\tGene\tTranscriptID\tOriginal\tIndividuallyAddable\tMaxAdditional\tExtended\tIsMaximal\tStatus"])
    for r in canonical:
        summary.append("\t".join(str(r.get(k, "")) for k in ["Dataset", "Gene", "TranscriptID", "OriginalSelectedSize",
            "NumIndividuallyAddable", "MaxAdditionalTogether", "ExtendedIndependentSetSize", "IsMaximal", "ResultStatus"]))
    summary.extend(["", "Completed solutions per cohort:"])
    counts = {}
    for r in records:
        key = (r["Dataset"], r["TranscriptID"])
        counts[key] = counts.get(key, 0) + (r["ResultStatus"] == "COMPLETE")
    summary.extend(f"{dataset}\t{chrom}\t{count}" for (dataset, chrom), count in counts.items())
    reasons = {}
    for r in records:
        if r["ResultStatus"] != "COMPLETE":
            key = (r["ResultStatus"], r.get("Notes", ""))
            reasons[key] = reasons.get(key, 0) + 1
    summary.extend(["", "Unfinished/unverifiable/invalid reasons:"])
    summary.extend(f"{count}\t{status}\t{reason}" for (status, reason), count in reasons.items())
    summary.append("Research inputs are opened read-only. No original solution, notebook, matrix, figure, or manuscript is written.")
    if final:
        changed = [p for p, before in json.loads((OUTPUT / "input_sources.json").read_text()).items() if stamp(p) != before]
        summary.append(f"Final original-input size/mtime verification: {'PASS' if not changed else changed}")
    (OUTPUT / "run_summary.txt").write_text("\n".join(summary) + "\n")
    return records


def self_check():
    """Run the small synthetic checks; print successes and return None or raise.

    Compare exact_extension with exhaustive enumeration on 108 seeded small graphs,
    including restricted candidate masks. A four-pattern fixture then checks AA
    compatibility, read mapping, fixed-S extension and witness IDs. This is not a
    rerun of the real datasets. The function WRITES synthetic_proteins.tsv to OUTPUT;
    redirect OUTPUT to a temporary directory when preserving existing artifacts.
    """
    rng = np.random.RandomState(1892)
    for size in range(1, 10):
        for repeat in range(12):
            adjacency = [0] * size
            for i in range(size):
                for j in range(i):
                    if rng.rand() < 0.35:
                        adjacency[i] |= 1 << j
                        adjacency[j] |= 1 << i
            candidate = (1 << size) - 1 if repeat % 2 else int(rng.randint(1 << size))
            optimum = max(bit_count(mask) for mask in range(1 << size) if not mask & ~candidate
                          if all(not adjacency[i] & mask for i in bit_indices(mask)))
            witness, exact, _ = exact_extension(candidate, adjacency.__getitem__, time.time() + 5)
            check(exact and len(witness) == optimum, "Synthetic exact-MIS validation failed")
    print("PASS: 108 labelled synthetic graphs, exact MIS agrees with exhaustive enumeration", flush=True)
    columns = ["Protein", "Reads", "NumOfReads", *[f"metadata_{i}" for i in range(11)], "0:3(K)", "3:6(A)"]
    patterns = [("K", "A"), ("K,R", "A,G"), ("R", "G"), ("N", "A")]
    path = OUTPUT / "synthetic_proteins.tsv"
    rows = [{**dict.fromkeys(columns, "0"), "Protein": str(i), "Reads": f"r{i}",
             "NumOfReads": "1", "0:3(K)": a, "3:6(A)": b} for i, (a, b) in enumerate(patterns)]
    write_tsv(path, rows, columns)
    graph = build_graph(path, True)
    for i in range(len(patterns)):
        for j in range(len(patterns)):
            direct = i != j and all(set(a.split(",")) & set(b.split(","))
                                   for a, b in zip(patterns[i], patterns[j]))
            check(bool(graph[4](i) & (1 << j)) == direct, "Synthetic AA compatibility mismatch")
    result, witness = analyze_solution(dict(UniqueSamples="0", NumUniqueSamples="1"),
        graph, graph[3], 2, 5, time.time() + 5)
    check(result["ResultStatus"] == "COMPLETE" and result["MaxAdditionalTogether"] == 2,
          "Synthetic fixed-set extension failed")
    check(set(witness) == {"2", "3"}, "Synthetic witness IDs mismatch")
    print("PASS: AA-set compatibility, read membership, fixed-set extension and direct witness checks", flush=True)


def main():
    """Parse command-line arguments and run the requested workflow; return None.

    The default command prepares/resumes the inventory and processes canonical
    solutions before remaining rows. --prepare-only still writes inventory/exports;
    it is not read-only. --self-check invokes only the synthetic checks and writes
    their fixture. A normal resume skips completed graph calculations but refreshes
    generated metadata/tables. The foreground process has no promised disconnect
    survival; committed SQLite rows provide the recovery point after interruption.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--memory-gib", type=float, default=32)
    parser.add_argument("--canonical-seconds", type=float, default=3600)
    parser.add_argument("--other-seconds", type=float, default=1200)
    parser.add_argument("--hours", type=float, default=24)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    started = time.time()
    deadline = started + args.hours * 3600
    limits = vars(args)
    command = " ".join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]])
    check(1 <= args.workers <= 8, "Use between 1 and 8 cohort workers")
    meminfo = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    available_gib = int(meminfo["MemAvailable"].split()[0]) / 1024**2
    workers = min(args.workers, len(os.sched_getaffinity(0)), int(0.70 * available_gib / args.memory_gib))
    check(workers >= 1, "Available RAM is insufficient for one worker within the aggregate limit")
    limits["effective_workers"] = workers
    (OUTPUT / "execution_contract.json").write_text(json.dumps(dict(command=command, host=os.uname().nodename,
        available_gib=available_gib, load=os.getloadavg(), limits=limits, started=started,
        lifecycle="Foreground process; disconnect survival is not assumed"), indent=2))
    cohorts, notes = build_inventory()
    export_results(cohorts, notes, started, command, limits)
    verify_octopus_eligibility(cohorts, deadline)
    if args.prepare_only:
        export_results(cohorts, notes, started, command, limits, final=True)
        return
    # First reproduce the two requested PacBio1 canonical cases, then prioritize
    # all other canonicals before opening the remaining-solution phase.
    solver_cohorts = [c for c in cohorts if c["HasOriginalSolverResult"]]
    try:
        for cohort in solver_cohorts[:2]:
            print(worker(cohort, "canonical", limits, deadline), flush=True)
            records = export_results(cohorts, notes, started, command, limits)
            result = next(r for r in records if r["SolutionID"] == cohort["CanonicalID"])
            print("PacBio1 comparison", cohort["Gene"], {k: result.get(k) for k in
                ["OriginalSelectedSize", "NumIndividuallyAddable", "MaxAdditionalTogether", "ResultStatus"]}, flush=True)
        context = multiprocessing.get_context("spawn")
        for phase in ["canonical", "remaining"]:
            jobs = [c for c in solver_cohorts if phase != "canonical" or c.get("CanonicalID")]
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
                pending = {pool.submit(worker, cohort, phase, limits, deadline): cohort for cohort in jobs}
                last_export = time.time()
                while pending:
                    done, _ = wait(pending, timeout=10, return_when=FIRST_COMPLETED)
                    for future in done:
                        cohort = pending.pop(future)
                        print(phase, future.result(), flush=True)
                    if done and phase == "canonical":
                        export_canonical()
                    if time.time() - last_export > 60:
                        export_results(cohorts, notes, started, command, limits)
                        last_export = time.time()
                    if time.time() >= deadline:
                        for future in pending:
                            future.cancel()
                        break
            if time.time() >= deadline:
                notes.append("24-hour/invocation execution budget reached; unattempted rows remain INCOMPLETE")
                break
    finally:
        export_results(cohorts, notes, started, command, limits, final=True)


if __name__ == "__main__":
    main()
