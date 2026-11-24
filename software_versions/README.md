# Software provenance captured on 23 November 2026

These files describe the **installed environments on alu15**, captured without installing, updating or replacing an environment. They are not a claim that every listed package was used at every historical processing stage. See [the historical stage table](historical_stages.2026-11-23.md) for that distinction.

For each of `combinatorics`, `combinatorics2` and `pacbiocomb`, the dated files contain:

- `explicit.txt`: `conda list --prefix <actual-prefix> --explicit`, including conda build URLs.
- `json`: `conda list --prefix <actual-prefix> --json`.
- `environment.yml`: `conda env export --prefix <actual-prefix>`.
- `python-metadata.json`: distribution names, versions and installer metadata read with that prefix's Python. `pacbiocomb` has no Python executable, so this item is explicitly not applicable.

The prefixes were `/home/alu/kobish/anaconda3/envs/<name>`. Explicit conda lists alone do not capture pip-installed distributions; use the JSON/export/Python metadata together. Absolute prefixes identify the inspected installation, not a required destination for a new user. The historical root files `combinatorics.explicit.txt` and `pacbiocomb.explicit.txt` remain unchanged. Package metadata does not prove that binary artifacts can still be downloaded or that a clean installation will solve successfully.

`executables.2026-11-23.json` records direct executable version calls. Prefix `bin` was placed first in PATH for FastQC so its Perl interpreter came from the intended environment. `host.2026-11-23.json` records Rocky Linux 9.8, Linux x86_64 and glibc 2.34. The older README's CentOS 7 description concerns the historical installation, not this validation host. macOS and Windows were not tested.

`julia.2026-11-23.json` records Julia 1.7.2, the tracked Project/Manifest hashes, manifest pins and the corresponding installed source paths for 397 entries. Installed source paths were present. It records both manifest and source Project versions without forcing them to agree: JLL build revisions and some package version spellings differ. Source trees and artifact bytes were not independently rehashed. The custom BioAlignments repository/revision is recorded in the inventory. This was read-only TOML/path inspection, not `Pkg.instantiate`, package resolution or a clean installation. A first `Pkg.dependencies()` inspection encountered a read-only depot log and was replaced by the non-writing inspection.

The November scientific runtime includes Python 3.13.14, NumPy 2.4.6, pandas 3.0.3 and scikit-learn 1.9.0. Current executable versions can differ from historical requirements: current `combinatorics2` samtools is 1.23.1 and current `pacbiocomb` CCS is 6.4.0. Neither replaces the historical 1.14/1.18/1.22.1 or CCS 6.2.0/8.2.0 records.

Start with [README.v2.md](../README.v2.md) for input acquisition, environment specifications, producer/consumer order and external assets. The dated inventories supplement those instructions. No clean-install timing was measured. Full manuscript reproduction requires native data, historical stage versions, retained inputs and permitted fonts/manual artwork; a small demo is a separate test. See [the current demo validation note](demo_validation.2026-11-23.md).

The dated SHA-256 manifest identifies these new files. Full commands, receipts, native headers and validation limitations are included in `NATURE_FINAL_REPORTING_AND_S16_REVIEW.zip`. The archive is review evidence; the scientific producer code remains tracked under `Code/`.
