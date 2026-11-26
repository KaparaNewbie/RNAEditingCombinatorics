# Running the demo

The complete tracked small-BAM demo was **verified on 24 November 2026**, using source commit `2955c17eee69b15e9c121b92de190ea1bfac983b`. Alignment, pileup/translation and Julia selection completed in **249.75 seconds (4.16 minutes)** on alu15 with four allowed CPU cores. This is a server-based desktop-resource proxy, not a measured ordinary-desktop runtime or a universal five-minute promise. An earlier attempt on 23 November disabled Julia compiled-module caching (`--compiled-modules=no`) and was stopped by its five-minute audit timer. That timeout was not an algorithm failure or a comparable runtime benchmark; the superseded report remains in Git history.

## Inputs and existing environments

Run from the repository root. The 10,918,280-byte `Demo/PCLO-CNS-RESUB.C0x1291.ccs.bam` and `Demo/DataTable.Squid.LongReads.Demo.csv` are tracked demonstration inputs. A normal repository clone includes them; do not rerun the historical subsampling recipe below. Obtain the **original** squid FASTA and corresponding editing-site BED using [Reference data and annotations](../README.md#reference-data-and-annotations), including the original editing-site workbook-to-BED conversion. The external reference/BED are not bundled with the BAM. Do not substitute the extended Jan2025 files; verify downloaded/generated copies against these validated input identities:

| Input | SHA-256 |
|---|---|
| `Demo/PCLO-CNS-RESUB.C0x1291.ccs.bam` | `8eb77c82b565c529931d00f5d0059dc4a41be49dfc76c76bb1f63513da39d6ce` |
| `D.pealeii/Annotations/orfs_squ.fa` | `49f094d8f83aabb5961760353aca9d2d349c185a19a14451e78dd0c45f1c6e09` |
| `D.pealeii/Annotations/D.pea.EditingSites.bed` | `7d5787218c2664806e0024e344c65379970f2887b75d464005f16f5a9d265471` |

Use already installed environments described in [README.md](../README.md) and [software_versions](../software_versions/README.md). The validation used Python 3.13.14 (`combinatorics2`), samtools 1.18, pbmm2 26.1.99 (`pacbiocomb`), MultiQC 1.35 and Julia 1.7.2 with the repository Project/Manifest. The alignment script activates the existing `pacbiocomb` environment under `~/anaconda3` by default; on another host explicitly set its `pacbio --base_conda_env_dir`, `--pb_conda_env_name` and `--pbmm2_path` options to the installed locations. No package installation is part of this recipe.

## Sequential commands

Choose a new output directory. Set the executable/depot paths and allowed CPU IDs for your machine; the values below are the measured alu15 configuration. Normal Julia compiled-module caching is enabled. A private first depot keeps any new cache files separate; offline mode assumes the dependencies already exist in the second depot.

```bash
set -euo pipefail
PYTHON=/home/alu/kobish/anaconda3/envs/combinatorics2/bin/python
SAMTOOLS=/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools
MULTIQC=/home/alu/kobish/anaconda3/envs/combinatorics2/bin/multiqc
JULIA=/home/alu/kobish/Julia/julia-1.7.2/bin/julia
JULIA_SHARED_DEPOT=/home/alu/kobish/.julia
DEMO_CPUS=0-3
DEMO_WORK="$PWD/Demo/new_validation_run"  # must not already exist
mkdir "$DEMO_WORK"
mkdir "$DEMO_WORK/alignment" "$DEMO_WORK/pileup" "$DEMO_WORK/distinct" \
      "$DEMO_WORK/tmp" "$DEMO_WORK/julia_overlay"
export PATH="$(dirname "$PYTHON"):$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 MULTIQC_DISABLE_VERSION_CHECK=1
export TMPDIR="$DEMO_WORK/tmp" MPLCONFIGDIR="$DEMO_WORK/tmp/matplotlib"
export JULIA_DEPOT_PATH="$DEMO_WORK/julia_overlay:$JULIA_SHARED_DEPOT" JULIA_PKG_OFFLINE=true

# Alignment: retain the original reference and its paired editing annotations.
taskset -c "$DEMO_CPUS" "$PYTHON" Code/align.py \
  --genome D.pealeii/Annotations/orfs_squ.fa --in_dir Demo \
  --out_dir "$DEMO_WORK/alignment" --processes 1 --threads 4 \
  --samtools_path "$SAMTOOLS" --multiqc_path "$MULTIQC" \
  pacbio --known_sites_bed_file D.pealeii/Annotations/D.pea.EditingSites.bed \
  > "$DEMO_WORK/alignment.log" 2>&1

# Update only a COPY of the sample table to the newly generated aligned BAM.
"$PYTHON" - "$DEMO_WORK" <<'PYTABLE'
import csv, sys
from pathlib import Path
work = Path(sys.argv[1]).resolve()
rows = list(csv.DictReader(Path("Demo/DataTable.Squid.LongReads.Demo.csv").open()))
assert len(rows) == 1
bam = work / "alignment/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam"
assert bam.is_file()
rows[0]["Path"] = str(bam)
with (work / "DataTable.Squid.LongReads.Demo.csv").open("x") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
PYTABLE

# Pileup/translation: wait for successful completion before selection.
taskset -c "$DEMO_CPUS" "$PYTHON" Code/pileup_with_subparsers.py \
  --transcriptome D.pealeii/Annotations/orfs_squ.fa \
  --known_editing_sites D.pealeii/Annotations/D.pea.EditingSites.bed \
  --exclude_flags 2304 --parity SE --min_rq 0.998 --min_bq 30 \
  --out_dir "$DEMO_WORK/pileup" --processes 1 --threads 4 \
  --samtools_path "$SAMTOOLS" --gz_compression directed_sequencing_data \
  --data_table "$DEMO_WORK/DataTable.Squid.LongReads.Demo.csv" \
  > "$DEMO_WORK/pileup.log" 2>&1

# Full selection: all reads, five fractions, four fraction repetitions,
# two algorithms and two algorithm repetitions. No background jobs/workers.
taskset -c "$DEMO_CPUS" "$JULIA" --startup-file=no --project=. --threads 4 \
  Code/Simulations/maximal_independent_set_5.jl \
  --infiles "$DEMO_WORK/pileup/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz" \
  --postfix_to_remove .C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz \
  --idcol Protein --firstcolpos 15 --datatype Proteins \
  --outdir "$DEMO_WORK/distinct" --fracstep 0.2 --fracrepetitions 4 \
  --algrepetitions 2 --algs Ascending Descending --run_solve_threaded \
  > "$DEMO_WORK/selection.log" 2>&1

# Inspect the final tab-separated CSV (one timestamped file in this fresh directory).
"$PYTHON" - "$DEMO_WORK" <<'PYRESULT'
import csv, sys
from pathlib import Path
files = list((Path(sys.argv[1]) / "distinct").glob("*.DistinctUniqueProteins.*.csv"))
assert len(files) == 1
rows = list(csv.DictReader(files[0].open(), delimiter="\t"))
assert len(rows) == 80
print(files[0])
for row in rows[-6:]:
    print(*(row[k] for k in ["Fraction", "FractionRepetition", "Algorithm",
                             "AlgorithmRepetition", "NumUniqueSamples"]), sep="\t")
PYRESULT
```

These shell commands run in the foreground. The measured audit additionally supervised process-tree memory, disk and elapsed time (32 GiB, 20 GiB, four-hour safety ceiling), retaining stage logs and exit status. Quiet compilation is not evidence of failure. Configure equivalent resource supervision for a long run; the shell recipe alone does not enforce memory/disk/time ceilings. No GPU or cluster scheduler was required.

## Validated outputs and timing

The native intermediate tables contain 6,582 original reads, 5,487 unique RNA patterns, 5,487 translated RNA-pattern records and 5,339 unique protein patterns. Original multiplicities sum to 6,582 in each aggregated table.

The final file is `distinct/PCLO-CNS-RESUB.DistinctUniqueProteins.<timestamp>.csv` (tab-separated): `Fraction`, `FractionRepetition`, `Algorithm`, `AlgorithmRepetition`, `NumUniqueSamples`, `UniqueSamples`, `AvailableReads`. It has **80 rows** (5 × 4 × 2 × 2). The final full-fraction counts in the first validated run were 3,093–3,097 for Ascending and 3,099–3,103 for Descending. Every selected identifier/cardinality, sampled-read multiplicity and independent-set property passed an independent native-table check. These heuristic sets are not certified globally maximum; independence does not itself establish maximality.

| Measured stage | Seconds |
|---|---:|
| Alignment | 20.13 |
| Pileup and translation | 108.57 |
| Julia selection | 120.91 |
| End to end, including supervisor overhead | 249.75 |
| Second selection run, same input and existing cache | 115.74 |

The second selection also produced 80 valid rows; its full-fraction counts span 3,093–3,104. Both runs used four allowed cores on an Intel Xeon E5-2650 v4 at 2.20 GHz, Rocky Linux 9.8, Julia four threads and one BLAS thread. Peak sampled process-tree PSS was 2.27 GiB for the full pipeline and 1.34 GiB for the second selection. These are measurements, not certified minimum hardware requirements. Existing installed Julia caches were available in both fresh processes; the private overlay stayed empty, so this is not a cold-depot benchmark. Startup/loading and JIT overhead are included in stage timings.

Selection is stochastic: default `consistentfracsampling=false` samples without a fixed global RNG, and heuristic tie breaking also uses randomness. The CLI's default `randseed=1892` does not seed all these operations. Repeated count/identifier equality is therefore not an acceptance criterion. Do not reduce fractions/repetitions, reseed the method or force agreement with historical counts.

## Historical example, retained for provenance

The following six rows are the historical example's last six full-fraction algorithm results, not six genes or six biological replicates. Their exact values are not a deterministic expected-output test. The corresponding first validated run gave **3100, 3102, 3094, 3095, 3099, 3103**. Count differences are consistent with the inspected stochastic sampling/tie-breaking implementation; the historical environment and RNG state are not recoverable from this example, so their exact contribution cannot be isolated.

| Fraction | FractionRepetition | Algorithm | AlgorithmRepetition |	NumUniqueSamples |
| -------- | ------- | -------- | ------- | -------- | 
| 1.0 | 3 | Descending | 1 |	3102 |
| 1.0 | 3 | Descending | 2 |	3102 |
| 1.0 | 4 | Ascending | 1 |	3093 |
| 1.0 | 4 | Ascending | 2 |	3094 |
| 1.0 | 4 | Descending | 1 |	3103 |
| 1.0 | 4 | Descending | 2 |	3102 |

# BTS: setting the demo

```bash
mkdir -p Demo/Data/CCS

samtools view \
--with-header \
--subsample 0.1 \
--subsample-seed 1892 \
--threads 20 \
--bam \
--output Demo/PCLO-CNS-RESUB.C0x1291.ccs.bam \
D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam
```