# Historical environment specifications and package snapshots

These files were moved from the repository root with their contents unchanged. They are preserved for software provenance and historical reproducibility, and must not automatically be interpreted as the current recommended installation environments.

Use the [repository-root README](../../README.md) for maintained installation instructions. The active specifications remain at the root: [combinatorics2_simplified.yml](../../combinatorics2_simplified.yml), [pacbiocomb.yml](../../pacbiocomb.yml), [Project.toml](../../Project.toml) and [Manifest.toml](../../Manifest.toml).

| Archived file | Purpose |
|---|---|
| [combinatorics.yml](combinatorics.yml) | Earlier broad Conda specification for the original `combinatorics` Python analysis and sequencing-tool environment. |
| [combinatorics2.yml](combinatorics2.yml) | Earlier `combinatorics2` Conda specification, retained separately from the current simplified specification. |
| [l_giremi.yml](l_giremi.yml) | Historical environment specification for the L-GIREMI workflow, including its pip dependency and sequencing utilities. |
| [combinatorics.explicit.txt](combinatorics.explicit.txt) | Historical package-filename/build snapshot for `combinatorics`; it contains filenames rather than a complete Conda URL lockfile. |
| [pacbiocomb.explicit.txt](pacbiocomb.explicit.txt) | Historical package-filename/build snapshot for PacBio processing tools; it is not a complete Conda URL lockfile. |

The dated inventories captured on 23 November 2026 remain directly under [software_versions/](../README.md). Historical execution versions can differ from those November 2026 installed versions. See the [stage-specific provenance table](../historical_stages.2026-11-23.md) for evidence and limits; a package's inclusion in a specification or snapshot alone does not prove that it generated a particular scientific result.

[Logs.Final.md](../../Logs.Final.md) preserves the original historical command body. Commands there that name `combinatorics.yml` or `l_giremi.yml` at the old root location are historical quotations; their files are now at `software_versions/old_versions/combinatorics.yml` and `software_versions/old_versions/l_giremi.yml`, relative to the repository root. No historical workflow or environment installation was rerun during this move.
