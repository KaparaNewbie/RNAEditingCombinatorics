# Publication figure and supplementary-table reproduction

## Purpose and prerequisites

This guide reconstructs the manuscript figures and supplementary workbooks from retained scientific results. For sequencing-data processing and new analyses, follow the [main workflow](../README.md#analyses); the [scientific notebook table](../README.md#scientific-notebooks-and-workflow) and [notebook reconstruction command](../README.md#reconstruct-the-notebook-representations) identify the relevant sources.

Use the Python analysis environment and Julia 1.7.2 described in [Installation and software requirements](../README.md#installation-and-software-requirements). Installed versions and clean-install limitations are recorded in [software_versions/README.md](../software_versions/README.md) and the [stage-specific version table](../software_versions/historical_stages.2026-11-23.md). Figure export also needs Chrome for Kaleido, `rsvg-convert`, Poppler, permitted Open Sans/Arial fonts, and the Docker image containing PyMuPDF described under [external assets](#special-external-assets-and-figure-specific-requirements).

Run shell commands from the repository root unless stated otherwise. Python snippets below run from `Code/Notebooks` and share the setup variables introduced in the first snippet. Prepare scientific tables before running their figure sections; the renderer does not repeat all upstream calculations.

The computational figures and supplementary tables have been reproduced from retained scientific inputs in the documented installed environments. This does not cover a fresh raw-sequencing-to-figure run or a complete clean-environment installation. The external inputs and assets listed below remain necessary; numerical reductions can differ at machine precision from historical exports.

## Required retained scientific inputs

A fresh clone supplies code and selection records, but not the large native sequencing, solution or assigned-expression files. Obtain the matching inputs from the authors and retain the relative paths used by the scientific notebooks:

- Original and additional squid PacBio, Illumina, pooled Octopus and single-cell position/read/unique-RNA/protein tables, distinct solutions and assigned-expression files for the documented populations.
- The original and Jan2025 references and annotations, single-cell pseudo-read/cell metadata, and neural-expression annotations specified in the [main workflow](../README.md#reference-data-and-annotations).
- Scientific exports from the dataset notebooks, retained depth solutions, simulation results and the original SVG geometry listed below.
- Workbook templates, the ScanProsite response, manual artwork and permitted fonts listed in the relevant sections below.

[selected_solutions.tsv](../Code/Notebooks/manuscript_input_data/selected_solutions.tsv) records 25 publication selections: ordinary assignments for two PB1 genes and 19 Illumina genes, plus the two PB1 genes under each of the Miyata and Grantham100 comparisons. It records the native source paths and hashes, solution parameters, selected-protein counts and ordered protein identities; it does not contain those native data. Fresh runs of `maximal_independent_set_5.jl` and `expressionlevels.jl` may select different heuristic solutions and are not substitutes for these retained assignments.

The scientific notebooks' **Reproducing manuscript results** tables specify each branch's inputs, cohort, schema, prerequisites and exports. Follow their indicated sections rather than running every exploratory cell. Use the same input versions and population definitions for cached and newly exported tables.

## Configure output locations

Choose a new table root and a separate figure-output root for each run. Producers and the renderer share `MANUSCRIPT_FIGURE_SOURCES` through `manuscript_paths.py`; assignment caches and RNA matrices go in `MANUSCRIPT_WORK`. The default table root is `Code/Notebooks/figure_sources/new_run`; the example below uses a unique run ID to keep new exports separate from existing results.

```bash
# Run from the repository root; choose a unique RUN_ID.
RUN_ID=your_run_id
export MANUSCRIPT_DATA_ROOT="$PWD"  # or the repository containing retained native data
export MANUSCRIPT_FIGURE_SOURCES="$PWD/Code/Notebooks/figure_sources/runs/$RUN_ID"
export MANUSCRIPT_FINAL_FIGURES="$PWD/Code/Notebooks/final_figures/runs/$RUN_ID"
export MANUSCRIPT_WORK="$PWD/Code/Notebooks/manuscript_intermediates/$RUN_ID"
export MANUSCRIPT_ASSIGNMENTS="$MANUSCRIPT_WORK/assignments"
export MANUSCRIPT_FONTS=/path/to/permitted/OpenSans/fonts
export MANUSCRIPT_ARIAL=/path/to/permitted/Arial.TTF
export NUMBA_CACHE_DIR="$MANUSCRIPT_WORK/numba_cache"
export MANUSCRIPT_FIGURE2_ASSETS=/path/to/approved/Figure2/assets
# Do not precreate producer output subdirectories; the producers guard creation.
```

Start Jupyter from the configured shell, or set the same variables in the notebook
before importing the manuscript modules. A running kernel does not inherit later
shell-variable changes. Select the intended Python environment explicitly.

## Preparing publication-selected data

### Extract the publication-selected assignments

[manuscript_inputs.py](../Code/Notebooks/manuscript_inputs.py) reads the retained native assigned-expression files, verifies their SHA-256 hashes, solution parameters, selected-protein counts and ordered identities, then exports standardized tables for the 25 recorded selections and combined ordinary-assignment tables. It does not rerun distinct-isoform selection or recalculate read-to-isoform assignments.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_inputs.py \
  --project-root "$MANUSCRIPT_DATA_ROOT" --output "$MANUSCRIPT_ASSIGNMENTS"
```

### Export expression summaries

[manuscript_expression_tables.py](../Code/Notebooks/manuscript_expression_tables.py) reads those selected assignments and produces plot-ready expression, rank-fit, cumulative-expression and common/rare-isoform tables for Figures 4a, S8, S9 and S17.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_expression_tables.py \
  --assignments "$MANUSCRIPT_ASSIGNMENTS" --output "$MANUSCRIPT_FIGURE_SOURCES/shared_expression"
```

## Generating figure input tables

### Export shared and cross-dataset tables

[various_plots.BQ30.ipynb](../Code/Notebooks/various_plots.BQ30.py) provides cross-platform table exports. Alternatively, call the adapters in `manuscript_plot_tables.py` from `Code/Notebooks` as below. Follow the scientific notebooks' execution tables for upstream producers; use one export route per table and a destination that does not yet exist.

```python
from pathlib import Path
from manuscript_paths import FIGURE_SOURCES, ASSIGNMENTS_DIR
from manuscript_inputs import PROJECT
from manuscript_plot_tables import (
    figure2_table, figure3_tables, figure6_tables, export_native_s3, figure_s3_tables,
    figure_s1_tables, figure_s6_tables, figure_s10_tables, figure_s12_tables,
    figure_s13_tables, figure_s16_table, dependency_matrix_tables, assemble_fisher_tables,
)
native_notebooks = PROJECT / 'Code/Notebooks'
work = ASSIGNMENTS_DIR.parent
figure2_table(PROJECT, FIGURE_SOURCES / 'Figure2/PCLO_top100_unique_reads.tsv')
figure3_tables(native_notebooks, ASSIGNMENTS_DIR, FIGURE_SOURCES / 'Figure3', PROJECT)
figure6_tables(native_notebooks, FIGURE_SOURCES / 'Figure6')
figure_s1_tables(native_notebooks, FIGURE_SOURCES / 'FigureS1')
export_native_s3(PROJECT, work / 's3')
figure_s3_tables(work / 's3', FIGURE_SOURCES / 'FigureS3')
figure_s6_tables(native_notebooks, ASSIGNMENTS_DIR, FIGURE_SOURCES / 'FigureS6', PROJECT)
figure_s10_tables(native_notebooks, ASSIGNMENTS_DIR, FIGURE_SOURCES / 'FigureS10')
figure_s12_tables(native_notebooks, FIGURE_SOURCES / 'FigureS12')
figure_s13_tables(native_notebooks, PROJECT, FIGURE_SOURCES / 'FigureS13')
# Run the documented BY, S18 and S16 notebook branches before these handoffs:
dependency_matrix_tables(native_notebooks, work / 'by', FIGURE_SOURCES)
assemble_fisher_tables(work / 's18', FIGURE_SOURCES / 'FigureS18')
(FIGURE_SOURCES / 'FigureS16').mkdir(parents=True, exist_ok=True)
figure_s16_table(work / 's16', FIGURE_SOURCES / 'FigureS16/figure_s16_summaries.tsv')
```

### Dataset-specific scientific exports

The Figure 6 cell-class producer is [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.ipynb](../Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.NeuronalVsNonneuronal.py), including its `FixedTMR50` export.

The separate [o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.ipynb](../Code/Notebooks/o.vul_SC_tmr50_vs_tmr1000.BQ30.BHAfterNoise.3.py) uses the entropy-considered assignments and original pseudo-read/cell-annotation mappings. Its full-cohort Fisher/BH section produces the manuscript's 856-gene association analysis; follow its execution table for selection, annotation, original-support expansion and common filters. Top100 analyses are separate diagnostics.

[noise_positions_and_reads_squid.ipynb](../Code/Notebooks/noise_positions_and_reads_squid.py) joins noise positions to original read mappings. The squid and pooled-octopus notebooks produce 12-mismatch position tables under `Code/Notebooks/12mm/v2/positions/`. The “12-mismatch analysis” section of [various_plots.BQ30.ipynb](../Code/Notebooks/various_plots.BQ30.py) produces sample/panel counts for Figures S2 and S15. The analysis unit is a mismatch site within its sample/transcript population, with population-specific denominators and undefined-proportion handling.

The PB1 domain analysis uses the `revision-20261109-domains-*` cells and the retained ScanProsite response. Additional-replicate recurrence and relative abundance are calculated within each gene and biological replicate. The Illumina notebook uses its 19-gene `conditions`/`chroms` selection, excluding `VGA_BPS13`; preserve final editing flags and the cohort filters in the pooled and cell-class notebooks.

### Producer-to-renderer mapping

`F` below denotes the selected `MANUSCRIPT_FIGURE_SOURCES` root. Each computational
row is consumed by the corresponding named section of [final_manuscript_figures.ipynb](../Code/Notebooks/final_manuscript_figures.py).
See each scientific notebook's execution table for prerequisites and native producers.

| Figure/table | Scientific producer | Generated handoff under F | Final consumer |
|---|---|---|---|
| 1 | Author manual artwork | External asset, not generated data | Author assembly |
| 2 | PB1 unique RNA; figure2_table | Figure2/PCLO_top100_unique_reads.tsv | figure2-asset / figure2-assembly |
| 3 | PB1/Illumina; figure3_tables | Figure3/squid_original_figure3{a,b,c,d}_plot.tsv | figure3-* |
| 4 | expression_tables; additional squid replicates, batch-export-figure4b | shared_expression/figure4a_expression.tsv; Figure4/figure4b_recurrence.tsv | figure4-* |
| 5 | PB1/Illumina All-reads entropy; BY cells + dependency_matrix_tables | Figure5/{figure5c_entropy,figure5_pclo_matrices.BY}.tsv | figure5-* |
| 6 | Pooled/cell-class notebooks; figure6_tables | Figure6/figure6_{cohorts,distributions,mean_markers}.tsv | figure6-* |
| S1 | Native read/noise exports; figure_s1_tables; violin_inputs | FigureS1/figure_s1_* | figure-s1-* |
| S2 | manuscript_s2; mismatch producers | FigureS2/s2_*; shared_mismatches/* | figure-s2-* |
| S3 | Original-read recoding; export_native_s3 → figure_s3_tables | FigureS3/figure_s3_{counts,statistics}.tsv | figure-s3-* |
| S4 | Additional squid replicates, batch-export-figure-s4 | FigureS4/figure_s4_plot.tsv | figure-s4-* |
| S5 | Pooled batch-export-figure-s5 | FigureS5/figure_s5_plot.tsv | figure-s5-* |
| S6 | Native depth/AA/dispersion; figure_s6_tables | FigureS6/figure_s6* | figure-s6-* |
| S7 | manuscript_s7 matrices → scores → combine | FigureS7/s7_scores.tsv | s7-complete-updated-20261116-* |
| S8/S9 | manuscript_expression_tables | shared_expression/figure_s{8,9}* | figure-s8 / figure-s9 |
| S10 | PB1 fractions; additional squid replicates, separate barcode/UMI exports; figure_s10_tables | FigureS10/figure_s10_* | figure-s10-* |
| S11 | PB1 BY + native MI; dependency_matrix_tables | FigureS11/figure_s11_matrices.BY.tsv | figure-s11-* |
| S12 | Native Corrs.PB1/IL1; figure_s12_tables | FigureS12/figure_s12_{pairs,boxes,outliers}.tsv | figure-s12-* |
| S13 | Pooled selection/read membership; figure_s13_tables | FigureS13/figure_s13_* | figure-s13-* |
| S14 | Pooled s14-native-handoff / batch-export-figure-s14 | FigureS14/figure_s14_{cohorts,curves,test}.tsv | figure-s14-* |
| S15 | Noise and VP mm12 cells | shared_mismatches/combined_s2_s15_mm12_* | s15-* |
| S16 | empirical_graph_assessment; figure_s16_table | FigureS16/figure_s16_summaries.tsv | figure-s16-* |
| S17 | expression_tables; violin_inputs.s17_statistics_geometry | shared_expression/figure_s17_observations.tsv; FigureS17/* | figure-s17-* |
| S18 | PB1/Illumina native S18 Fisher/BH cells; assemble_fisher_tables | FigureS18/figure_s18_fisher.tsv and summaries | figure-s18-* |
| S19 | Code/Simulations/manuscript_s19.jl | FigureS19/s19_{correlations,histograms}.tsv | figure-s19-* |
| S20 | manuscript_s20, original RNA and fixed assignments | FigureS20/figure_s20_{expression,boxes,chimeras}.tsv | figure-s20-* |
| Tables S1/S3 | manuscript_metrics; original/additional squid PacBio native stages | New metrics/workbooks outside F | manuscript_workbooks.export_pacbio_metrics_workbook |
| Table S2 | Illumina original-read sample SD cell | New TableS2_recoding_sample_SD.tsv outside F | export_sample_sd_workbook |
| Table S4 | S18 gene Fisher/BH summary | FigureS18/s18_fisher_gene_summary.tsv | export_fisher_workbook |

For S6, `figure_s6_tables()` in [manuscript_plot_tables.py](../Code/Notebooks/manuscript_plot_tables.py) reads the retained depth solutions in [s6_depth_inputs.tsv](../Code/Notebooks/manuscript_input_data/s6_depth_inputs.tsv), dataset-notebook dispersion exports, and PCLO assignment/amino-acid inputs. The related depth analysis is in [various_plots.BQ30.ipynb](../Code/Notebooks/various_plots.BQ30.py). For S10, `figure_s10_tables()` combines the PB1 fraction comparison with the separate barcode-filtered and deduplicated notebook exports. Both adapters write tables consumed by [final_manuscript_figures.ipynb](../Code/Notebooks/final_manuscript_figures.py).

### Figure S7

S7 uses seven initialized OpenBLAS threads for distance calculations and one internal fitting thread to retain the numerical reduction order. The documented runtime uses libopenblas 0.3.28 (`pthreads_h94d23a6_1`).

#### Prepare amino-acid matrices

```bash
OPENBLAS_NUM_THREADS=7 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s7.py prepare \
  --project-root "$MANUSCRIPT_DATA_ROOT" --assignments "$MANUSCRIPT_ASSIGNMENTS" --output "$MANUSCRIPT_WORK/s7_matrices"
```

#### Calculate dimensionality-reduction scores

```bash
for branch in PacBio_ordinary Illumina_ordinary PacBio_dissimilar; do
  OPENBLAS_NUM_THREADS=7 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s7.py scores \
    --matrices "$MANUSCRIPT_WORK/s7_matrices" --branch "$branch" --output "$MANUSCRIPT_WORK/$branch.tsv"
done
```

#### Combine platform scores

```bash
mkdir -p "$MANUSCRIPT_FIGURE_SOURCES/FigureS7"
python Code/Notebooks/manuscript_s7.py combine --branches \
  "$MANUSCRIPT_WORK/PacBio_ordinary.tsv" "$MANUSCRIPT_WORK/Illumina_ordinary.tsv" "$MANUSCRIPT_WORK/PacBio_dissimilar.tsv" \
  --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS7/s7_scores.tsv"
```

### Figure S19: Simulate editing and export tables

The Figure S19 producer [manuscript_s19.jl](../Code/Simulations/manuscript_s19.jl) implements the simulation from [unifinished_editing_isoforms.jl](../Code/Simulations/unifinished_editing_isoforms.jl), preserving its loop and random-draw order. Use Julia 1.7.2, StatsBase 0.34.3, one thread and seed 1892, with N = 10,000 reads, K = 10 sites and the nine ordered editing/removal probability pairs from {0.001, 0.002, 0.003}. Each panel contains 45 unique off-diagonal correlations and 15 right-open histogram bins; `nextfloat(maximum)` includes the maximum in the upper bin. The renderer reads `s19_correlations.tsv` and `s19_histograms.tsv`.

```bash
OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 Code/Simulations/manuscript_s19.jl \
  "$MANUSCRIPT_FIGURE_SOURCES/FigureS19" "$MANUSCRIPT_WORK/s19"
```

### Figure S20: Analyze chimeric reads and export tables

[visualize_chimeric_reads_origin_squid_2.ipynb](../Code/Notebooks/visualize_chimeric_reads_origin_squid_2.py) describes the analysis. For the manuscript export, use [manuscript_s20.py](../Code/Notebooks/manuscript_s20.py) with the directories configured above. It joins original unique RNA patterns to fixed selected proteins, using common/rare fractions X = 0.01 and Y = 0.30, original multiplicities, ties-to-even rounding and prefix/suffix witnesses. Any eligible original RNA can support a rare-protein chimera call. Expression is normalized in native source-row order before sorting and alignment to canonical selected IDs.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s20.py \
  --project-root "$MANUSCRIPT_DATA_ROOT" --assignments "$MANUSCRIPT_ASSIGNMENTS" \
  --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS20" --work "$MANUSCRIPT_WORK/s20"
```

If the retained data are in the repository selected by `MANUSCRIPT_DATA_ROOT` (or the code repository when that variable is unset), the S20 command can use the default data root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s20.py \
  --assignments "$MANUSCRIPT_ASSIGNMENTS" \
  --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS20" --work "$MANUSCRIPT_WORK/s20"
```

### Export plot geometry and expression statistics

S1/S2/S17 use SVG contours and jitter from the original plotting cells alongside numerical data from scientific tables. Supply the four matching SVG files in an external `GEOMETRY` directory:

| Figure | SVG basename | SHA-256 |
|---|---|---|
| S1 | Per chrom noise levels - PacBio1 vs. Illumina vs. octopus pooled.svg | `93d6706671f4782ef3647132d2862fd8728e7c4698916519a5d484c9dc32c744` |
| S2 | Distribution of pooled editing levels according to their editing status in current vs. previous study - Squid.svg | `a4bfa94291863f51e0e5a234b22fe8b56a2c3c4c9f59ad50e9bbc0103fc49b6c` |
| S17 | Mean edited positions vs expression levels - PacBio1.svg | `0998b5b1538251b007aba8915805f5e25ce3bd9a44fb94ea47a8ba2b2d7f4dbe` |
| S17 | Mean edited positions vs expression levels - Illumina.svg | `50e4e39430b76cfb2867759b495e09b686fa75a47a35fdc932684abd991e4b15` |

From `Code/Notebooks`, after the native noise and expression exports:

```python
from manuscript_figures.violin_inputs import s1_noise_geometry, s17_statistics_geometry
geometry = Path('/path/to/permitted/historical/SVGs')
s1_noise_geometry(native_notebooks,
    geometry / 'Per chrom noise levels - PacBio1 vs. Illumina vs. octopus pooled.svg',
    FIGURE_SOURCES / 'FigureS1')
s17_statistics_geometry(FIGURE_SOURCES / 'shared_expression/figure_s17_observations.tsv',
    [geometry / f'Mean edited positions vs expression levels - {platform}.svg'
     for platform in ['PacBio1', 'Illumina']], FIGURE_SOURCES / 'FigureS17')
```

### Figure S2

For S2, run from the repository root. `sites` projects native positions for three distinct populations: squid, octopus and octopus_motif. `tables` combines those projections with the reference/known-site annotations and distribution SVG. The symbolic links below collect this run's projections:

#### Project native editing sites

```bash
GEOMETRY=/path/to/permitted/historical/SVGs
for branch in squid octopus octopus_motif; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s2.py sites \
    --branch "$branch" --project-root "$MANUSCRIPT_DATA_ROOT" \
    --output "$MANUSCRIPT_WORK/s2_sites/$branch"
done
```

#### Combine projections and export tables

```bash
mkdir "$MANUSCRIPT_WORK/s2_sites/combined"
ln -s "$MANUSCRIPT_WORK"/s2_sites/{squid,octopus,octopus_motif}/*.tsv "$MANUSCRIPT_WORK/s2_sites/combined/"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python Code/Notebooks/manuscript_s2.py tables \
  --sites "$MANUSCRIPT_WORK/s2_sites/combined" --output "$MANUSCRIPT_FIGURE_SOURCES/FigureS2" \
  --project-root "$MANUSCRIPT_DATA_ROOT" \
  --distribution-svg "$GEOMETRY/Distribution of pooled editing levels according to their editing status in current vs. previous study - Squid.svg"
```

## Rendering the figures

Open [final_manuscript_figures.ipynb](../Code/Notebooks/final_manuscript_figures.py) with `Code/Notebooks` as its working directory. Follow its execution table for shared setup, selected sections and cropping. Set `FIGURES_TO_BUILD`; Figure 3 also requires `RENDER_FIGURE3`. Run the selected input, drawing and export cells, leaving archived comparisons and collection-refresh cells unused. Figures 5 and S7 have dedicated crops; Figure 2 requires assembly. Intermediate canvases go in `render_work`, with final outputs in `svg`, `pdf` and `previews`. S11/S18 use Matplotlib physical cropping and Poppler rasterization at 600 dpi; other general page exports use MuPDF.

## Supplementary tables and workbooks

[manuscript_workbooks.py](../Code/Notebooks/manuscript_workbooks.py) writes new supplementary workbooks using the supplied templates, retaining manual fields, labels, primers and author metadata. Supply `review_inputs/9.11.2026/Supplementary_tables_12_revised_completed.xlsx` and the S2 template where required. S1/S3 use `manuscript_metrics.py --project-root ... --output NEW_DIR`; S2 uses original-read sample SD (`ddof=1`); S4 uses the S18 Fisher/BH summary. Keep numerical outputs at full precision.

## Special external assets and figure-specific requirements

**Figure 1** is supplied as assembled artwork: `/Dropbox COMB paper/Figures/Figures 12/Final figures/Figure1.pdf` and `Figure1.png` (13 October 2026). Use the PDF as the separate figure and the PNG as its preview. Their SHA-256 hashes are `be386af49114b0adb2335b53fea1be23d90619bfae028d1e2b82854f2b82dd52` and `cf84ab4940f841ccd8cf947721b40fec02cf5f9245c94fee707adcaa91bdb584`, respectively.

**Figure 2** requires `BioRender_base.pdf` and `BioRender_reference.pdf`, with hashes specified in [docker_export.py](../Code/Notebooks/manuscript_figures/docker_export.py), and matching permitted Arial fonts. Set `MANUSCRIPT_ARIAL` to the Arial TTF file and `MANUSCRIPT_FIGURE2_ASSETS` to the directory containing the two PDFs, outside the generated-table root. The `figure2-asset-*` cells generate panel b from the PCLO rank table before assembly. Proprietary fonts and manual artwork are not distributed with the repository.

The rendering image is `combinatorics-figure2:2026-10-13`, ID `sha256:e1fcbeb34fbe87a26e76d553e8987277b29e3792eba52b7555f96fc1d5a47c97`. Its [Dockerfile](../Code/Notebooks/manuscript_figures/Dockerfile) depends on the external `python-stats:1.0` base, whose public availability is unconfirmed. Obtain the matching image from its maintainer. Rendering containers use no network, read-only inputs, two CPUs and 4 GiB RAM.

The frozen [Jan2025 squid FASTA](../reference_data/Jan2025/orfs_squ.fa) is distributed with the repository; see [reference data and annotations](../README.md#reference-data-and-annotations) for its identity and use. Native sequencing and processed inputs, other required references/annotations, single-cell metadata and fixed solution files must be obtained separately. Domain annotation uses `Code/Notebooks/figure_sources/scanprosite_accepted.json` (SHA-256 `c02bc0b3e0563aab6590a06efc405e1a6fae0f013885ec5e2d8cd82513ad4ed5`). Its original query date/database release is unknown; provide this response at the documented path, separately from newly generated tables.
