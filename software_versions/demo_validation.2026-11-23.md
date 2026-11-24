# Bounded demo validation, 23 November 2026

The tracked 10,918,280-byte `PCLO-CNS-RESUB.C0x1291.ccs.bam` is a real small demonstration dataset. At source commit `2955c17eee69b15e9c121b92de190ea1bfac983b`, one isolated run completed alignment and pileup/translation, but did **not** complete Julia isoform selection within a five-minute total cap. The historical “about 5 minutes” claim and six final selection counts in [DEMO.md](../Demo/DEMO.md) are not verified for this current setup.

The test used the existing environments on alu15, Rocky Linux 9.8 x86_64: Python 3.13.14 (`combinatorics2`), samtools 1.18 explicitly selected from `combinatorics`, pbmm2 26.1.99 from `pacbiocomb`, MultiQC 1.35, and Julia 1.7.2. See [complete version inventories](README.md). Four CPUs were allowed, numerical-library threads were restricted to one, and the 24 GiB memory ceiling was not approached: measured peak process-tree PSS was 2.27 GiB. Alignment took 19.12 seconds and pileup 109.58 seconds. Julia was stopped after its remaining 172.74 seconds; total elapsed time including termination was 301.93 seconds.

To preserve the installed Julia depot, the test used `--compiled-modules=no` and a new empty writable depot overlay ahead of the existing read-only depot, with offline package mode. No package installation or environment update occurred. This cold-load condition differs from a warm compiled Julia session; the result does not establish a minimum runtime or disprove a five-minute runtime on every machine. No desktop install/runtime test was performed.

| Completed output | Rows | Analysis unit / check |
|---|---:|---|
| reads | 6,582 | Original reads retained for the demo analysis |
| unique_reads | 5,487 | Unique nucleotide patterns; `NumOfReads` sums to 6,582 |
| proteins | 5,487 | Translated unique nucleotide patterns; original read multiplicities sum to 6,582 |
| unique_proteins | 5,339 | Unique protein patterns; original read multiplicities sum to 6,582 |

These are intermediate tables, **not** selected distinct isoforms. No final selection CSV was produced before the timeout. The historical final counts 3102, 3102, 3093, 3094, 3103 and 3102 are documented examples, not a newly validated exact-output assertion.

## Current command corrections and prerequisites

Obtain the original squid reference `D.pealeii/Annotations/orfs_squ.fa` and paired `D.pea.EditingSites.bed` as described in [README.md](../README.md#reference-data-and-annotations). They are external inputs, not bundled by the small BAM. Use the repository root as the working directory, install/resolve the documented environments before the demo, and verify executable versions. Do not silently replace the reference with the extended Jan2025 version.

The old alignment command omits the now-required `pacbio --known_sites_bed_file` argument. The old Julia source path `Code/UnorderedNaNDepletion/maximal_independent_set_5.jl` no longer exists; the tracked path is `Code/Simulations/maximal_independent_set_5.jl`. Use a new output directory and update a **copy** of the demo data-table `Path` column to its new aligned BAM. Run stages sequentially in the foreground. The old `nohup ... &` pileup recipe does not itself wait for completion before the selection command.

The validated alignment/pileup invocation structure was:

```bash
# Set PYTHON, SAMTOOLS, MULTIQC and DEMO_WORK to existing executables/new output path.
"$PYTHON" Code/align.py --genome D.pealeii/Annotations/orfs_squ.fa \
  --in_dir Demo --out_dir "$DEMO_WORK/alignment" --processes 1 --threads 4 \
  --samtools_path "$SAMTOOLS" --multiqc_path "$MULTIQC" \
  pacbio --known_sites_bed_file D.pealeii/Annotations/D.pea.EditingSites.bed

# First copy Demo/DataTable.Squid.LongReads.Demo.csv and change its Path to
# $DEMO_WORK/alignment/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam (an absolute path).
"$PYTHON" Code/pileup_with_subparsers.py \
  --transcriptome D.pealeii/Annotations/orfs_squ.fa \
  --known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
  --exclude_flags 2304 --parity SE --min_rq 0.998 --min_bq 30 \
  --out_dir "$DEMO_WORK/pileup" --processes 1 --threads 4 \
  --samtools_path "$SAMTOOLS" --gz_compression directed_sequencing_data \
  --data_table "$DEMO_WORK/DataTable.Squid.LongReads.Demo.csv"
```

For the subsequent Julia stage, retain the original demo thresholds, fraction step 0.2, four fraction repetitions, two algorithm repetitions, Ascending/Descending choices, Protein ID and first feature column 15. Correct only the source/input/output paths and resource settings. Its current full completion and expected counts remain to be measured before claiming an end-to-end demo or a desktop runtime. The review archive records the exact attempted command and timeout; the original demo and approved outputs were preserved.

The demo uses CPU tools and required no GPU or cluster scheduler. This is separate from reproducing the full manuscript, which requires much larger native inputs, retained scientific tables/solutions and the documented external figure assets. Measured demo memory is not a full-pipeline resource requirement or a certified minimum.

## Follow-up: complete validation, 24 November 2026

**VERIFIED.** A fresh isolated run at the same source revision completed all three stages with normal Julia compiled-module caching and a private writable depot overlay before the existing installed depot. Alignment took 20.13 seconds, pileup/translation 108.57 seconds and Julia selection 120.91 seconds; the full measured run was 249.75 seconds (4.16 minutes). A second full selection from the identical newly generated table passed in 115.74 seconds. Both processes could read existing compiled caches; the private overlay stayed empty, so neither is a cold-depot benchmark.

All 80 rows in each run contain valid IDs/cardinalities, the complete fraction/repetition grid and valid independent sets under an independently constructed native AA-compatibility graph. Native multiplicities remain 6,582 original reads, 5,487 RNA/protein-pattern records and 5,339 unique proteins. The first full-fraction counts span 3,093–3,103; the second spans 3,093–3,104. The six corresponding historical-example rows in the first run are 3100, 3102, 3094, 3095, 3099 and 3103. Global RNG sampling and heuristic tie breaking are stochastic; no deterministic identity with the historical example is asserted.

The existing software versions above were verified again. On alu15 (Intel Xeon E5-2650 v4, 2.20 GHz), both runs used CPU affinity 0–3, Julia four threads and one BLAS thread. Peak sampled process-tree PSS was 2.27 GiB for the pipeline and 1.34 GiB for the repeat. The supervised bounds were four cumulative demo hours, 32 GiB process-tree PSS and 20 GiB new audit files; none was approached. This is a server-based desktop-resource proxy, not an ordinary-desktop measurement. No installation or source/environment change was made. The source BAM/reference/BED hashes remained unchanged.

The original November 23 five-minute attempt above remains historical evidence of a timer stop, not a Julia algorithm crash. Its `--compiled-modules=no` condition is different from the new runs. Current foreground commands, source paths, copied sample-table handling, prerequisites and final output interpretation are now in [Demo/DEMO.md](../Demo/DEMO.md). Full reports, both output tables, hashes, resource samples and exact commands are local audit artifacts under `Code/Notebooks/scientific_audit_2026-10-14/final_statistics_and_demo_2026-11-24/`; this ignored audit directory is not promised as part of a clean clone.
