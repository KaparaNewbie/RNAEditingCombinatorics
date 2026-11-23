"""One selected generated-table root shared by scientific producers and figures.

Set MANUSCRIPT_FIGURE_SOURCES to a new absolute directory before importing.
Native data, assignment caches and external artwork stay outside this root.
"""
import os
from pathlib import Path

NOTEBOOK_DIR = Path(__file__).resolve().parent
FIGURE_SOURCES = Path(os.environ.get(
    'MANUSCRIPT_FIGURE_SOURCES', NOTEBOOK_DIR / 'figure_sources/new_run')).resolve()
ASSIGNMENTS_DIR = Path(os.environ.get(
    'MANUSCRIPT_ASSIGNMENTS', NOTEBOOK_DIR / 'manuscript_intermediates/new_run/assignments')).resolve()


def figure_source(name, suffix='.tsv'):
    """Return the documented final-table path; never search historical directories."""
    if name in {'figure4a_expression', 'figure_s8_s9_expression', 'figure_s8_fits',
                'figure_s8_equations', 'figure_s9_markers', 'figure_s17_observations'}:
        folder = 'shared_expression'
    elif name.startswith('combined_s2_s15_'):
        folder = 'shared_mismatches'
    elif name.startswith('squid_original_figure3'):
        folder = 'Figure3'
    elif name == 'PCLO_top100_unique_reads':
        folder = 'Figure2'
    elif name.startswith('s2_'):
        folder = 'FigureS2'
    elif name.startswith('s7_'):
        folder = 'FigureS7'
    elif name.startswith('s19_'):
        folder = 'FigureS19'
    elif name.startswith('figure_s'):
        number = name[len('figure_s'):].split('_')[0].rstrip('abc')
        assert number.isdigit(), name
        folder = 'FigureS' + number
    elif name in {'figure4b_recurrence', 'figure5_pclo_matrices.BY', 'figure5c_entropy',
                  'figure6_cohorts', 'figure6_distributions', 'figure6_mean_markers'}:
        folder = 'Figure' + name[len('figure')]
    else:
        raise ValueError(f'No manuscript table contract for {name!r}')
    return FIGURE_SOURCES / folder / (name + suffix)
