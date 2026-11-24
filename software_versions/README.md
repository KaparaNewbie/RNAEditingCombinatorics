# Software environments and provenance

Start with the [repository-root README](../README.md) for maintained installation instructions, input acquisition, producer/consumer order and external assets. Current installed versions and the historical versions used to generate scientific results are not necessarily identical.

## A. Current installation specifications

These files remain at the repository root:

- [combinatorics2_simplified.yml](../combinatorics2_simplified.yml): current broad Python analysis specification.
- [pacbiocomb.yml](../pacbiocomb.yml): PacBio processing tools.
- [Project.toml](../Project.toml) and [Manifest.toml](../Manifest.toml): Julia project and recorded dependency pins.

Follow the stage-specific version requirements in the root README. The Conda specifications are not immutable lockfiles, and their presence does not establish that a clean installation has been tested.

## B. Dated inventories captured on 23 November 2026

These files describe the **installed environments on alu15**, captured without installing, updating or replacing an environment. They are not a claim that every listed package was used at every historical processing stage. See [the historical stage table](historical_stages.2026-11-23.md) for that distinction.

For each of `combinatorics`, `combinatorics2` and `pacbiocomb`, the dated files contain:

- `explicit.txt`: `conda list --prefix <actual-prefix> --explicit`, including conda build URLs.
- `json`: `conda list --prefix <actual-prefix> --json`.
- `environment.yml`: `conda env export --prefix <actual-prefix>`.
- `python-metadata.json`: distribution names, versions and installer metadata read with that prefix's Python. `pacbiocomb` has no Python executable, so this item is explicitly not applicable.

The prefixes were `/home/alu/kobish/anaconda3/envs/<name>`. Explicit conda lists alone do not capture pip-installed distributions; use the JSON/export/Python metadata together. Absolute prefixes identify the inspected installation, not a required destination for a new user. The earlier filename-only snapshots are now [old_versions/combinatorics.explicit.txt](old_versions/combinatorics.explicit.txt) and [old_versions/pacbiocomb.explicit.txt](old_versions/pacbiocomb.explicit.txt), with their contents preserved. Package metadata does not prove that binary artifacts can still be downloaded or that a clean installation will solve successfully.

`executables.2026-11-23.json` records direct executable version calls. Prefix `bin` was placed first in PATH for FastQC so its Perl interpreter came from the intended environment. `host.2026-11-23.json` records Rocky Linux 9.8, Linux x86_64 and glibc 2.34. The CentOS 7 description in [Logs.Final.md](../Logs.Final.md) concerns the historical installation, not this validation host. macOS and Windows were not tested.

`julia.2026-11-23.json` records Julia 1.7.2, the tracked Project/Manifest hashes, manifest pins and the corresponding installed source paths for 397 entries. Installed source paths were present. It records both manifest and source Project versions without forcing them to agree: JLL build revisions and some package version spellings differ. Source trees and artifact bytes were not independently rehashed. The custom BioAlignments repository/revision is recorded in the inventory. This was read-only TOML/path inspection, not `Pkg.instantiate`, package resolution or a clean installation. A first `Pkg.dependencies()` inspection encountered a read-only depot log and was replaced by the non-writing inspection.

The November scientific runtime includes Python 3.13.14, NumPy 2.4.6, pandas 3.0.3 and scikit-learn 1.9.0. Current executable versions can differ from historical requirements: current `combinatorics2` samtools is 1.23.1 and current `pacbiocomb` CCS is 6.4.0. Neither replaces the historical 1.14/1.18/1.22.1 or CCS 6.2.0/8.2.0 records.

The dated inventories supplement the [repository-root README](../README.md). No clean-install timing was measured. Full manuscript reproduction requires native data, historical stage versions, retained inputs and permitted fonts/manual artwork; a small demo is a separate test. See [the current demo validation note](demo_validation.2026-11-23.md).

The [dated SHA-256 manifest](SHA256SUMS.2026-11-23.json) records the original November 23 snapshots. Later documentation edits, including the November 24 demo follow-up and this directory guide, do not update those historical hashes; the dated machine-readable inventory files remain unchanged. Full commands, receipts, native headers and validation limitations are included in `NATURE_FINAL_REPORTING_AND_S16_REVIEW.zip`. The archive is review evidence; the scientific producer code remains tracked under `Code/`.

## C. Historical environments

[old_versions/README.md](old_versions/README.md) describes the five archived YAML specifications and package-filename snapshots formerly stored at the repository root. They are preserved for provenance and historical reproducibility; they are not automatically recommended installation environments. The historical development log is [Logs.Final.md](../Logs.Final.md), whose original command body is preserved with a notice explaining relocated paths. The dated November 2026 inventories remain directly in this directory.

## D. Stage-specific software provenance

[historical_stages.2026-11-23.md](historical_stages.2026-11-23.md) distinguishes versions supported by native execution records from versions merely present in the November 2026 installations. Neither a historical specification nor a current inventory proves which packages executed a particular stage. Historical workflows were not rerun as part of this documentation reorganization.
