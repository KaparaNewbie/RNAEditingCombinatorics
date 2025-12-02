# Squid vs Octopus ORF homolog analysis using reciprocal BLAST
import subprocess
import tempfile
import os
import shutil
import re
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd


# Global dictionaries to store UniProt mappings
squid_uniprot_map = {}
octopus_uniprot_map = {}


def parse_uniprot_from_header(header):
    """Parse UniProt accession from FASTA header"""
    # Look for pattern like 'sp|Q99MU3|DSRAD_MOUSE' or 'tr|A0A1Y3VGE5|A0A1Y3VGE5_OCTBI'
    uniprot_match = re.search(r"(?:sp|tr)\|[A-Z0-9]+\|([A-Z0-9_]+)", header)
    if uniprot_match:
        return uniprot_match.group(1)
    return None


def load_uniprot_mappings(squid_fasta_file, octopus_fasta_file):
    """Load UniProt mappings from FASTA files"""
    global squid_uniprot_map, octopus_uniprot_map

    print("Loading UniProt mappings from FASTA headers...")

    # Load squid mappings
    squid_count = 0
    for record in SeqIO.parse(squid_fasta_file, "fasta"):
        uniprot_id = parse_uniprot_from_header(record.description)
        if uniprot_id:
            squid_uniprot_map[record.id] = uniprot_id
            squid_count += 1

    # Load octopus mappings
    octopus_count = 0
    for record in SeqIO.parse(octopus_fasta_file, "fasta"):
        uniprot_id = parse_uniprot_from_header(record.description)
        if uniprot_id:
            octopus_uniprot_map[record.id] = uniprot_id
            octopus_count += 1

    print(
        f"Loaded {squid_count} squid UniProt mappings and {octopus_count} octopus UniProt mappings"
    )


def extract_uniprot_names(transcript_id):
    """Extract UniProt names from loaded mappings"""
    # Check squid mappings first, then octopus mappings
    if transcript_id in squid_uniprot_map:
        return squid_uniprot_map[transcript_id]
    elif transcript_id in octopus_uniprot_map:
        return octopus_uniprot_map[transcript_id]
    return None


def get_target_gene_status(transcript_id, target_squid_genes):
    """Determine which target gene this transcript corresponds to"""
    transcript_to_target = {
        "comp134400_c0_seq1_extended": "ADAR1",
        "comp141693_c0_seq1": "GRIA2",
        "comp141881_c0_seq3": "RUSC2",
        "comp141044_c0_seq2": "TRIM2",
        "comp140439_c0_seq1": "CA2D3",
        "comp126362_c0_seq1": "ABL",
        "comp141517_c0_seq1": "DGLA",
        "comp141840_c0_seq2": "K0513",
        "comp141640_c0_seq1": "KCNAS",
        "comp140987_c3_seq1": "ACHA4",
        "comp140910_c2_seq1": "ANR17",
        "comp136058_c0_seq1": "TWK7",
        "comp141378_c0_seq7": "SCN1",
        "comp141158_c1_seq2": "CACB2",
        "comp140712_c0_seq3": "RIMS2",
        "comp141882_c0_seq14": "PCLO",
        "comp141880_c1_seq3": "DOP1",
        "comp141565_c6_seq3": "IQEC1",
        "comp141684_c0_seq1": "CSKI1",
        "comp141532_c3_seq11": "MTUS2",
        "comp141574_c0_seq3": "ROBO2",
    }

    return transcript_to_target.get(transcript_id, None)


def parse_blast_results(blast_output_file):
    """Parse BLAST output into a pandas DataFrame"""
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
    ]

    try:
        df = pd.read_csv(blast_output_file, sep="\t", names=columns, header=None)
        return df
    except pd.errors.EmptyDataError:
        print(f"No BLAST results found in {blast_output_file}")
        return pd.DataFrame(columns=columns)
    except Exception as e:
        print(f"Error reading BLAST results from {blast_output_file}: {e}")
        return pd.DataFrame(columns=columns)


def create_comprehensive_report(
    squid_to_octopus_df,
    octopus_to_squid_df,
    reciprocal_df,
    target_squid_genes,
    target_squid_chroms,
    use_bitscore=False,
):
    """Create comprehensive report with best matches for all target genes, including reverse matches"""

    # Get reciprocal gene pairs
    reciprocal_genes = set()
    if not reciprocal_df.empty:
        reciprocal_genes.update(reciprocal_df["SquidGene"].tolist())
        reciprocal_genes.update(reciprocal_df["OctopusGene"].tolist())

    # Create mapping from gene to transcript ID
    gene_to_transcript = dict(zip(target_squid_genes, target_squid_chroms))
    transcript_to_gene = dict(zip(target_squid_chroms, target_squid_genes))

    # Process all target genes using their specific transcript IDs
    comprehensive_hits = []

    # Track which combinations we've already added to avoid duplicates
    added_pairs = set()

    # Part 1: Squid -> Octopus matches (original logic)
    for target_gene in target_squid_genes:
        squid_transcript_id = gene_to_transcript[target_gene]

        # Find ALL hits for this transcript (not just filtered ones)
        squid_hits = squid_to_octopus_df[
            squid_to_octopus_df["qseqid"] == squid_transcript_id
        ]

        # Sort by bitscore or coverage based on flag
        if use_bitscore:
            squid_hits = squid_hits.sort_values("bitscore", ascending=False)
        else:
            squid_hits = squid_hits.sort_values("qcov", ascending=False)

        # Extract UniProt name for squid gene
        squid_uniprot_full = extract_uniprot_names(squid_transcript_id)
        target_status = get_target_gene_status(squid_transcript_id, target_squid_genes)
        is_target_gene = target_status is not None

        if not squid_hits.empty:
            # Get the best hit (highest bitscore)
            best_hit = squid_hits.iloc[0]
            octopus_gene = best_hit["sseqid"]

            # Extract UniProt name for octopus gene
            octopus_uniprot_full = extract_uniprot_names(octopus_gene)

            # Check reciprocal status
            is_reciprocal = squid_transcript_id in reciprocal_genes

            # Get reverse hit info if it exists
            reverse_hits = octopus_to_squid_df[
                (octopus_to_squid_df["qseqid"] == octopus_gene)
                & (octopus_to_squid_df["sseqid"] == squid_transcript_id)
            ]

            if not reverse_hits.empty:
                # Choose best reverse hit based on flag (note: for reverse direction, we want best octopus coverage)
                if use_bitscore:
                    best_reverse = reverse_hits.loc[reverse_hits["bitscore"].idxmax()]
                else:
                    best_reverse = reverse_hits.loc[reverse_hits["qcov"].idxmax()]

                reverse_identity = best_reverse["pident"]
                reverse_coverage = best_reverse["qcov"]
            else:
                reverse_identity = None
                reverse_coverage = None
            # reverse_bitscore = (
            #     reverse_hits.iloc[0]["bitscore"] if not reverse_hits.empty else None
            # )

            pair_key = (squid_transcript_id, octopus_gene)
            added_pairs.add(pair_key)

            comprehensive_hits.append(
                {
                    "TargetGene": target_gene,
                    "IsTargetGene": is_target_gene,
                    "MatchedTarget": target_status,
                    "SquidGene": squid_transcript_id,
                    "SquidUniProt": squid_uniprot_full,
                    "OctopusGene": octopus_gene,
                    "OctopusUniProt": octopus_uniprot_full,
                    "IsReciprocal": is_reciprocal,
                    "MatchDirection": "Squid->Octopus",
                    "SquidToOctopusIdentity": best_hit["pident"],
                    "SquidToOctopusCoverage": best_hit["qcov"],
                    # "SquidToOctopusBitscore": best_hit["bitscore"],
                    "OctopusToSquidIdentity": reverse_identity,
                    "OctopusToSquidCoverage": reverse_coverage,
                    # "OctopusToSquidBitscore": reverse_bitscore,
                }
            )
        else:
            # No hits found - still include in report
            comprehensive_hits.append(
                {
                    "TargetGene": target_gene,
                    "IsTargetGene": is_target_gene,
                    "MatchedTarget": target_status,
                    "SquidGene": squid_transcript_id,
                    "SquidUniProt": squid_uniprot_full,
                    "OctopusGene": None,
                    "OctopusUniProt": None,
                    "IsReciprocal": False,
                    "MatchDirection": "No match",
                    "SquidToOctopusIdentity": None,
                    "SquidToOctopusCoverage": None,
                    # "SquidToOctopusBitscore": None,
                    "OctopusToSquidIdentity": None,
                    "OctopusToSquidCoverage": None,
                    # "OctopusToSquidBitscore": None,
                }
            )

    # Part 2: Octopus -> Squid matches (new logic)
    # Find octopus genes that have good matches to target squid genes
    octopus_to_target_count = 0
    added_from_octopus = 0

    for _, octopus_hit in octopus_to_squid_df.iterrows():
        octopus_gene = octopus_hit["qseqid"]
        squid_target = octopus_hit["sseqid"]

        # Check if this squid gene is one of our targets
        if squid_target in transcript_to_gene:
            octopus_to_target_count += 1
            target_gene = transcript_to_gene[squid_target]
            pair_key = (squid_target, octopus_gene)

            # Skip if we already added this pair from the squid->octopus direction
            if pair_key in added_pairs:
                continue

            added_from_octopus += 1

            # Extract UniProt names
            squid_uniprot_full = extract_uniprot_names(squid_target)
            octopus_uniprot_full = extract_uniprot_names(octopus_gene)

            # Check reciprocal status
            is_reciprocal = squid_target in reciprocal_genes

            # Get forward hit info if it exists
            forward_hits = squid_to_octopus_df[
                (squid_to_octopus_df["qseqid"] == squid_target)
                & (squid_to_octopus_df["sseqid"] == octopus_gene)
            ]

            if not forward_hits.empty:
                # Choose best forward hit based on flag (for forward direction, we want best squid coverage)
                if use_bitscore:
                    best_forward = forward_hits.loc[forward_hits["bitscore"].idxmax()]
                else:
                    best_forward = forward_hits.loc[forward_hits["qcov"].idxmax()]

                forward_identity = best_forward["pident"]
                forward_coverage = best_forward["qcov"]
            else:
                forward_identity = None
                forward_coverage = None
            # forward_bitscore = (
            #     forward_hits.iloc[0]["bitscore"] if not forward_hits.empty else None
            # )

            comprehensive_hits.append(
                {
                    "TargetGene": target_gene,
                    "IsTargetGene": True,
                    "MatchedTarget": target_gene,
                    "SquidGene": squid_target,
                    "SquidUniProt": squid_uniprot_full,
                    "OctopusGene": octopus_gene,
                    "OctopusUniProt": octopus_uniprot_full,
                    "IsReciprocal": is_reciprocal,
                    "MatchDirection": "Octopus->Squid",
                    "SquidToOctopusIdentity": forward_identity,
                    "SquidToOctopusCoverage": forward_coverage,
                    # "SquidToOctopusBitscore": forward_bitscore,
                    "OctopusToSquidIdentity": octopus_hit["pident"],
                    "OctopusToSquidCoverage": octopus_hit["qcov"],
                    # "OctopusToSquidBitscore": octopus_hit["bitscore"],
                }
            )

    print(
        f"Debug: Found {octopus_to_target_count} octopus->target hits, added {added_from_octopus} new ones"
    )

    df = pd.DataFrame(comprehensive_hits)

    # # Sort by target gene status, then by reciprocal status, then by best available bitscore
    # if not df.empty:
    #     # Use the maximum bitscore from either direction for sorting
    #     df["BestBitscore"] = df[
    #         ["SquidToOctopusBitscore", "OctopusToSquidBitscore"]
    #     ].max(axis=1)
    #     df = df.sort_values(
    #         ["IsTargetGene", "IsReciprocal",
    #          "BestBitscore"
    #          ],
    #         ascending=[False, False, False],
    #     )
    #     # Remove the temporary sorting column
    #     df = df.drop("BestBitscore", axis=1)

    # Don't drop IsTargetGene and MatchedTarget columns as they're needed for analysis
    # df = df.drop(columns=["IsTargetGene", "MatchedTarget"])

    return df


def run_blast(query_file, subject_file, output_file, blast_type="blastn"):
    """Run BLAST alignment"""
    cmd = [
        blast_type,
        "-query",
        query_file,
        "-subject",
        subject_file,
        "-out",
        output_file,
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-evalue",
        "1e-5",
        "-max_target_seqs",
        "50",
    ]
    subprocess.run(cmd, check=True)


def find_reciprocal_homologs_clean(
    squid_transcriptome_file,
    octopus_transcriptome_file,
    target_squid_chroms,
    target_squid_genes,
    min_identity=70,
    min_coverage=0.5,
    blast_output_dir=None,
    use_bitscore=False,
):
    """Find reciprocal best hits with cleaned output format"""

    # Create temporary directory for BLAST output
    if blast_output_dir is None:
        temp_dir = tempfile.mkdtemp()
    else:
        temp_dir = blast_output_dir
        os.makedirs(temp_dir, exist_ok=True)

    try:
        # Parse transcriptome lengths
        print("Loading transcriptome sequences...")
        squid_lengths = {}
        octopus_lengths = {}

        for record in SeqIO.parse(squid_transcriptome_file, "fasta"):
            squid_lengths[record.id] = len(record.seq)

        for record in SeqIO.parse(octopus_transcriptome_file, "fasta"):
            octopus_lengths[record.id] = len(record.seq)

        print(
            f"Loaded {len(squid_lengths)} squid transcripts and {len(octopus_lengths)} octopus transcripts"
        )

        # Load UniProt mappings from FASTA headers
        load_uniprot_mappings(squid_transcriptome_file, octopus_transcriptome_file)

        # Run BLAST in both directions
        print("Running reciprocal BLAST alignment...")
        squid_to_octopus_file = os.path.join(temp_dir, "squid_to_octopus.blast")
        print(f"Running BLAST: squid transcriptome -> octopus transcriptome...")
        run_blast(
            squid_transcriptome_file,
            octopus_transcriptome_file,
            squid_to_octopus_file,
            blast_type="blastn",
        )

        squid_to_octopus_df = parse_blast_results(squid_to_octopus_file)

        if squid_to_octopus_df.empty:
            print("No BLAST hits found from squid to octopus")
            return pd.DataFrame(), pd.DataFrame()

        octopus_to_squid_file = os.path.join(temp_dir, "octopus_to_squid.blast")
        print(f"Running BLAST: octopus transcriptome -> squid transcriptome...")
        run_blast(
            octopus_transcriptome_file,
            squid_transcriptome_file,
            octopus_to_squid_file,
            blast_type="blastn",
        )

        octopus_to_squid_df = parse_blast_results(octopus_to_squid_file)

        if octopus_to_squid_df.empty:
            print("No BLAST hits found from octopus to squid")
            return pd.DataFrame(), pd.DataFrame()

        # Add sequence lengths and calculate coverage
        print("Calculating sequence coverage...")
        squid_to_octopus_df["qlen"] = squid_to_octopus_df["qseqid"].map(squid_lengths)
        squid_to_octopus_df["slen"] = squid_to_octopus_df["sseqid"].map(octopus_lengths)
        squid_to_octopus_df["qcov"] = (
            squid_to_octopus_df["length"] / squid_to_octopus_df["qlen"]
        )
        squid_to_octopus_df["scov"] = (
            squid_to_octopus_df["length"] / squid_to_octopus_df["slen"]
        )

        octopus_to_squid_df["qlen"] = octopus_to_squid_df["qseqid"].map(octopus_lengths)
        octopus_to_squid_df["slen"] = octopus_to_squid_df["sseqid"].map(squid_lengths)
        octopus_to_squid_df["qcov"] = (
            octopus_to_squid_df["length"] / octopus_to_squid_df["qlen"]
        )
        octopus_to_squid_df["scov"] = (
            octopus_to_squid_df["length"] / octopus_to_squid_df["slen"]
        )

        # Apply filters for reciprocal analysis
        print(
            f"Applying filters: min_identity={min_identity}%, min_coverage={min_coverage}"
        )

        squid_to_octopus_filtered = squid_to_octopus_df[
            (squid_to_octopus_df["pident"] >= min_identity)
            & (squid_to_octopus_df["qcov"] >= min_coverage)
        ]

        octopus_to_squid_filtered = octopus_to_squid_df[
            (octopus_to_squid_df["pident"] >= min_identity)
            & (octopus_to_squid_df["qcov"] >= min_coverage)
        ]

        print(
            f"After filtering: {len(squid_to_octopus_filtered)} squid->octopus hits, {len(octopus_to_squid_filtered)} octopus->squid hits"
        )

        # Find reciprocal best hits
        reciprocal_hits = []

        for _, s2o_hit in squid_to_octopus_filtered.iterrows():
            squid_gene = s2o_hit["qseqid"]
            octopus_gene = s2o_hit["sseqid"]

            # Check if octopus gene has squid gene as best hit
            o2s_matches = octopus_to_squid_filtered[
                (octopus_to_squid_filtered["qseqid"] == octopus_gene)
                & (octopus_to_squid_filtered["sseqid"] == squid_gene)
            ]

            if not o2s_matches.empty:
                # Choose best match based on flag
                if use_bitscore:
                    best_o2s = o2s_matches.loc[o2s_matches["bitscore"].idxmax()]
                else:
                    best_o2s = o2s_matches.loc[o2s_matches["qcov"].idxmax()]

                # Extract UniProt names
                squid_uniprot_full = extract_uniprot_names(squid_gene)
                octopus_uniprot_full = extract_uniprot_names(octopus_gene)

                # Check if this is a targeted gene
                target_status = get_target_gene_status(squid_gene, target_squid_genes)
                is_target_gene = target_status is not None

                reciprocal_hits.append(
                    {
                        "IsTargetGene": is_target_gene,
                        "MatchedTarget": target_status,
                        "SquidGene": squid_gene,
                        "SquidUniProt": squid_uniprot_full,
                        "OctopusGene": octopus_gene,
                        "OctopusUniProt": octopus_uniprot_full,
                        "SquidToOctopusIdentity": s2o_hit["pident"],
                        "SquidToOctopusCoverage": s2o_hit["qcov"],
                        # "SquidToOctopusBitscore": s2o_hit["bitscore"],
                        "OctopusToSquidIdentity": best_o2s["pident"],
                        "OctopusToSquidCoverage": best_o2s["qcov"],
                        # "OctopusToSquidBitscore": best_o2s["bitscore"],
                    }
                )

        reciprocal_df = (
            pd.DataFrame(reciprocal_hits) if reciprocal_hits else pd.DataFrame()
        )

        # Sort reciprocal results
        if not reciprocal_df.empty:
            reciprocal_df = reciprocal_df.sort_values(
                [
                    "IsTargetGene",
                    #  "SquidToOctopusBitscore"
                ],
                ascending=[
                    False,
                    #    False
                ],
            )

        # Create comprehensive report (includes ALL target genes with best matches)
        comprehensive_df = create_comprehensive_report(
            squid_to_octopus_df,
            octopus_to_squid_df,
            reciprocal_df,
            target_squid_genes,
            target_squid_chroms,
            use_bitscore,
        )

        # Save results
        if blast_output_dir:
            reciprocal_output_file = os.path.join(
                blast_output_dir, "reciprocal_homologs_clean.csv"
            )
            comprehensive_output_file = os.path.join(
                blast_output_dir, "comprehensive_homologs_clean.csv"
            )

            if not reciprocal_df.empty:
                reciprocal_df.to_csv(reciprocal_output_file, index=False)
                print(f"Clean reciprocal homologs saved to {reciprocal_output_file}")

            comprehensive_df.to_csv(comprehensive_output_file, index=False)
            print(f"Clean comprehensive report saved to {comprehensive_output_file}")

            # Also save raw BLAST results
            squid_to_octopus_df.to_csv(
                os.path.join(blast_output_dir, "squid_to_octopus_raw.csv"), index=False
            )
            octopus_to_squid_df.to_csv(
                os.path.join(blast_output_dir, "octopus_to_squid_raw.csv"), index=False
            )

        return reciprocal_df, comprehensive_df

    except Exception as e:
        print(f"Error in reciprocal homolog analysis: {e}")
        return pd.DataFrame(), pd.DataFrame()

    finally:
        # Clean up temp directory if we created it
        if blast_output_dir is None and "temp_dir" in locals():
            shutil.rmtree(temp_dir, ignore_errors=True)


def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description="Squid vs Octopus ORF homolog analysis using reciprocal BLAST"
    )

    parser.add_argument(
        "--squid-fasta",
        default="/private6/projects/Combinatorics/D.pealeii/Annotations/Jan2025/orfs_squ.fa",
        help="Path to squid transcriptome FASTA file (default: %(default)s)",
    )

    parser.add_argument(
        "--octopus-fasta",
        default="/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.fa",
        help="Path to octopus transcriptome FASTA file (default: %(default)s)",
    )

    parser.add_argument(
        "--output-dir",
        default="/private6/projects/Combinatorics/SquidVsOctopusHomologs",
        help="Output directory for results (default: %(default)s)",
    )

    parser.add_argument(
        "--min-identity",
        type=float,
        default=70.0,
        help="Minimum identity threshold for BLAST hits (default: %(default)s%%)",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.5,
        help="Minimum coverage threshold for BLAST hits (default: %(default)s)",
    )

    parser.add_argument(
        "--use-bitscore",
        action="store_true",
        help="Use bitscore instead of coverage for selecting best matches (default: False, uses coverage)",
    )

    args = parser.parse_args()

    # Define target genes and transcript IDs
    target_squid_chroms = [
        "comp141693_c0_seq1",
        "comp134400_c0_seq1_extended",
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

    target_squid_genes = [
        "GRIA2",
        "ADAR1",
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

    # Validate input files
    if not os.path.exists(args.squid_fasta):
        print(f"Error: Squid FASTA file not found: {args.squid_fasta}")
        return 1

    if not os.path.exists(args.octopus_fasta):
        print(f"Error: Octopus FASTA file not found: {args.octopus_fasta}")
        return 1

    # Create output directory
    blast_dir = Path(args.output_dir)
    blast_dir.mkdir(parents=True, exist_ok=True)

    print("Running clean reciprocal analysis with UniProt parsing...")
    print(f"Squid FASTA: {args.squid_fasta}")
    print(f"Octopus FASTA: {args.octopus_fasta}")
    print(f"Output directory: {args.output_dir}")
    print(f"Min identity: {args.min_identity}%")
    print(f"Min coverage: {args.min_coverage}")
    print(f"Selection criterion: {'Bitscore' if args.use_bitscore else 'Coverage'}")

    # Run the clean analysis with updated functions
    clean_reciprocal_df, clean_comprehensive_df = find_reciprocal_homologs_clean(
        args.squid_fasta,
        args.octopus_fasta,
        target_squid_chroms,
        target_squid_genes,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        blast_output_dir=blast_dir,
        use_bitscore=args.use_bitscore,
    )

    print("\n" + "=" * 80)
    print("CLEAN ANALYSIS RESULTS")
    print("=" * 80)

    print(f"Target squid genes analyzed: {len(target_squid_genes)}")

    # Analyze comprehensive results (includes ALL target genes with best matches)
    if not clean_comprehensive_df.empty:
        total_genes = len(clean_comprehensive_df)
        genes_with_hits = len(
            clean_comprehensive_df[clean_comprehensive_df["OctopusGene"].notna()]
        )
        target_genes_with_hits = len(
            clean_comprehensive_df[
                (clean_comprehensive_df["IsTargetGene"] == True)
                & (clean_comprehensive_df["OctopusGene"].notna())
            ]
        )

        print(f"\nCOMPREHENSIVE RESULTS (includes best matches even if not great):")
        print(f"- Total target genes: {total_genes}")
        print(f"- Target genes with octopus hits: {target_genes_with_hits}")
        print(
            f"- Target gene success rate: {target_genes_with_hits/len(target_squid_genes)*100:.1f}%"
        )

    # Analyze reciprocal results (strict criteria)
    if not clean_reciprocal_df.empty:
        reciprocal_total = len(clean_reciprocal_df)
        reciprocal_targeted = len(
            clean_reciprocal_df[clean_reciprocal_df["IsTargetGene"] == True]
        )
        print(f"\nRECIPROCAL BEST HITS (strict criteria):")
        print(f"- Total reciprocal best hits: {reciprocal_total}")
        print(f"- Reciprocal hits for target genes: {reciprocal_targeted}")
        print(
            f"- Target gene reciprocal rate: {reciprocal_targeted/len(target_squid_genes)*100:.1f}%"
        )

        # Show targeted gene reciprocal hits
        print(f"\nRECIPROCAL HITS FOR TARGET GENES:")
        target_reciprocals = clean_reciprocal_df[
            clean_reciprocal_df["IsTargetGene"] == True
        ]
        if not target_reciprocals.empty:
            display_cols = [
                "MatchedTarget",
                "SquidUniProt",
                "SquidToOctopusIdentity",
                "OctopusToSquidIdentity",
            ]
            print(target_reciprocals[display_cols].to_string(index=False))
        else:
            print("No reciprocal hits found for target genes")
    else:
        print("\nNo reciprocal best hits found")

    # Show comprehensive results for target genes
    print(f"\nCOMPREHENSIVE TARGET GENE RESULTS:")
    target_comprehensive = clean_comprehensive_df[
        clean_comprehensive_df["IsTargetGene"] == True
    ]
    if not target_comprehensive.empty:
        display_cols = [
            "TargetGene",
            "SquidUniProt",
            "IsReciprocal",
            "MatchDirection",
            "SquidToOctopusIdentity",
            "OctopusToSquidIdentity",
            "OctopusGene",
        ]
        print(target_comprehensive[display_cols].to_string(index=False))

    print(f"\nFiles saved:")
    print(f"- Clean reciprocal homologs: {blast_dir}/reciprocal_homologs_clean.csv")
    print(f"- Clean comprehensive report: {blast_dir}/comprehensive_homologs_clean.csv")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    exit(main())
