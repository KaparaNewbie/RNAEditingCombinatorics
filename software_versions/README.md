# Software environments and provenance

Start with the [repository-root README](../README.md) for installation instructions, input acquisition, workflow order and external assets. Installation specifications, dated inventories and historical execution versions serve different purposes.

## A. Current installation specifications

These files remain at the repository root:

- [combinatorics2_simplified.yml](../combinatorics2_simplified.yml): current broad Python analysis specification.
- [pacbiocomb.yml](../pacbiocomb.yml): PacBio processing tools.
- [Project.toml](../Project.toml) and [Manifest.toml](../Manifest.toml): Julia project and recorded dependency pins.

Use the [stage-specific software table](historical_stages.2026-11-23.md) for historical execution requirements. The Conda specifications are not immutable lockfiles, and their presence does not establish that a clean installation has been tested.

## B. Dated inventories captured on 23 November 2026

The dated inventories record the software environments installed on the researchers' computing system on **23 November 2026**. They are installation snapshots, not necessarily the versions used in every historical analysis or the latest installations after that date. See the [historical stage table](historical_stages.2026-11-23.md) for stage-specific execution versions.

For each of `combinatorics`, `combinatorics2` and `pacbiocomb`, the dated files contain:

- `explicit.txt`: `conda list --prefix <actual-prefix> --explicit`, including conda build URLs.
- `json`: `conda list --prefix <actual-prefix> --json`.
- `environment.yml`: `conda env export --prefix <actual-prefix>`.
- `python-metadata.json`: distribution names, versions and installer metadata read with that prefix's Python. `pacbiocomb` has no Python executable, so this item is explicitly not applicable.

Explicit Conda lists alone do not capture pip-installed distributions; use the JSON, environment export and Python metadata together. Choose installation locations appropriate to your system. Earlier package-filename snapshots are archived as [old_versions/combinatorics.explicit.txt](old_versions/combinatorics.explicit.txt) and [old_versions/pacbiocomb.explicit.txt](old_versions/pacbiocomb.explicit.txt). Package metadata does not guarantee that binary artifacts remain downloadable or that a clean installation will resolve successfully.

[executables.2026-11-23.json](executables.2026-11-23.json) records executable versions; [host.2026-11-23.json](host.2026-11-23.json) records the operating-system environment: Rocky Linux 9.8, Linux x86_64 and glibc 2.34. The CentOS 7 description in [Logs.Final.md](../Logs.Final.md) concerns an earlier installation. macOS and Windows were not tested.

[julia.2026-11-23.json](julia.2026-11-23.json) records Julia 1.7.2, Project/Manifest hashes, dependency pins, installed source locations and the custom BioAlignments repository/revision. Manifest and source Project versions are recorded separately where their spellings or build revisions differ. This inventory does not verify source/artifact byte identity or establish that a clean Julia installation succeeds.

The November 23 snapshot includes Python 3.13.14, NumPy 2.4.6, pandas 3.0.3 and scikit-learn 1.9.0. Its `combinatorics2` samtools and `pacbiocomb` CCS versions differ from some historical requirements; use the [stage-specific table](historical_stages.2026-11-23.md) when reproducing those stages.

No clean-install timing was measured. Full manuscript reproduction requires native data, historical stage versions, retained inputs and permitted fonts/manual artwork. See [Demo/DEMO.md](../Demo/DEMO.md) for the maintained demo instructions and complete validation.

The [dated SHA-256 manifest](SHA256SUMS.2026-11-23.json) records the original November 23 snapshots. Its historical hashes are preserved, including the entry for the superseded demo report removed from the current tree. Documentation revised after that date may differ from its recorded versions; this manifest does not validate the current tree. The dated machine-readable inventories remain unchanged.

## C. Historical environments

[old_versions/README.md](old_versions/README.md) describes archived YAML specifications and package-filename snapshots. These preserve historical provenance and are not automatically recommended installation environments. Earlier development commands are retained in [Logs.Final.md](../Logs.Final.md).

## D. Stage-specific software provenance

[historical_stages.2026-11-23.md](historical_stages.2026-11-23.md) distinguishes versions supported by native execution records from versions present in the November 23 installations. Neither an installation specification nor a dated inventory alone establishes which packages executed a particular stage.

### Additional documented tool requirements

The workflow also specifies the following tools, which are not listed in the dated stage table:

| Stage | Documented tool versions |
|---|---|
| Octopus BED/GFF/GTF conversion | AGAT container `quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0` |
| Neural-expression annotation | Salmon 1.7.0; OrthoFinder 2.5.5 |

These requirements are retained from the workflow documentation and [historical commands](../Logs.Final.md); they do not establish the versions used in every historical execution.
