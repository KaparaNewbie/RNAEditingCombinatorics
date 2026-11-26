# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Reproducing manuscript results
#
# Use the retained manuscript inputs and selections described in `../../README.md`.
# These are dependency branches, not an instruction to Run All. Run the named setup
# and the required upstream sections before an export; skip unrelated exploratory plots.
# Set `MANUSCRIPT_FIGURE_SOURCES` to the same new absolute directory for every producer
# and `final_manuscript_figures`. Keep assignment caches outside that directory and set
# `MANUSCRIPT_ASSIGNMENTS` explicitly. `MANUSCRIPT_DATA_ROOT` selects the native data
# repository when source and data are separate. Existing saved displays are historical.
#
# Cell IDs below identify the saved author notebook. Jupytext reconstruction may assign new IDs; locate the same cells by their section headings and opening source comments.
#
# | Order / branch | Existing section and required cells | Upstream inputs | Generated outputs | Downstream consumer |
# |---|---|---|---|---|
# | 1 | `batch-selection`, `figure3-imports`, `figure3-paths-typography` | New MANUSCRIPT_FIGURE_SOURCES and MANUSCRIPT_FINAL_FIGURES; permitted MANUSCRIPT_FONTS | Selected figure gates and output directories | Only explicitly selected figures |
# | 2a | Figure3 input/layout/draw/export cells (`figure3-*`); then `final-collection-export` | New Figure3 tables | Figure3 SVG/PDF/PNG | Publication review |
# | 2b | `batch-setup`; selected figure load/draw/export cells; `final-collection-export` | New per-figure/shared tables; generated S1/S17 geometry | F4–F6 and S1–S6/S8–S14/S16–S20 | F5 has a dedicated exporter/crop |
# | 2c | `figure3-layout-helpers`, `supplement-layout`, `batch-setup`, `s15-*`, `final-collection-export` | New shared_mismatches tables | S15 SVG/PDF/PNG | Publication review |
# | 2d | `s7-complete-updated-20261116-*` | FigureS7/s7_scores.tsv from exact manuscript_s7 producer | S7 SVG/PDF/PNG | Publication review |
# | 2e | `batch-setup`; `figure2-asset-*`; `figure2-assembly-*` | New PCLO rank table; MANUSCRIPT_FIGURE2_ASSETS; permitted fonts and exact Docker image | Figure2b asset and assembled F2 | Manual artwork remains an external input |
#
# For reconstructed notebooks, these opening comments identify the shared setup:
# `Select only the figure sections to execute` (selection), `Plotting dependencies only`
# (imports), `One path/physical-unit configuration` (paths), `New figures share the
# established local fonts` (batch setup), `Axes rectangles are specified` (Figure 3/S15
# layout), and `Crop only newly rendered` (final crop). S15 additionally needs the
# `from matplotlib.transforms import Bbox` layout cell and its `Read prepared counts only`
# load cell. The S7 branch starts `Independent plotting setup for the sole canonical S7
# producer`. Run a selected figure's loading, drawing and export cells in source order.

# %%
# Select only the figure sections to execute; scientific producers run separately.
FIGURES_TO_BUILD = []
RENDER_FIGURE3 = False


# %% [markdown]
# The publication outputs are under `final_figures/`. All render gates default to off. Select explicit figure sections and set `FINAL_FIGURES_DIR` to a new subdirectory before execution. Generate scientific tables first using the producers in README.md. The renderer does not run scientific inference or refresh the historical collection. Figure 1 and the manual portions of Figure 2 require external artwork; matching fonts and the Figure 2 Docker image are external dependencies.
#

# %%
# Plotting dependencies only; the selected environment is explicit and checked.
import os
import sys
import json
import hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/combinatorics2-figure3-matplotlib')
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.text import Text
from IPython.display import display, Image
assert hasattr(pd.DataFrame, 'map')
print({'executable': sys.executable, 'prefix': sys.prefix,
       'pandas_version': pd.__version__, 'pandas_file': pd.__file__})

# %%
# One path/physical-unit configuration and one typography definition.
NOTEBOOK_DIR = Path.cwd()
from manuscript_paths import FIGURE_SOURCES, figure_source
FINAL_FIGURES_DIR = Path(os.environ.get('MANUSCRIPT_FINAL_FIGURES', NOTEBOOK_DIR / 'final_figures/new_run')).resolve()
INPUT_DIR = FIGURE_SOURCES / 'Figure3'
OUTPUT_DIR = FINAL_FIGURES_DIR / 'render_work/Figure3'
# Refuse to replace any already published selected figure.
selected_outputs = set(FIGURES_TO_BUILD) | ({'Figure3'} if RENDER_FIGURE3 else set())
for stem in selected_outputs:
    for folder, suffix in [('svg', 'svg'), ('pdf', 'pdf'), ('previews', 'png')]:
        assert not (FINAL_FIGURES_DIR / folder / f'{stem}.{suffix}').exists(), 'Choose a new FINAL_FIGURES_DIR'
for folder in ['svg', 'pdf', 'previews']:
    (FINAL_FIGURES_DIR / folder).mkdir(parents=True, exist_ok=True)
FONT_DIR = Path(os.environ['MANUSCRIPT_FONTS']).resolve()
WIDTH_MM, HEIGHT_MM = 180.0, 180.0
PT_PER_MM = 72 / 25.4
TEXT_PT = {'panel': 7, 'axis': 7, 'heading': 7, 'legend_heading': 7,
           'tick': 6, 'legend': 6, 'annotation': 6}
for font in ['OpenSans-Regular.ttf', 'OpenSans-Bold.ttf']:
    font_manager.fontManager.addfont(FONT_DIR / font)
assert Path(font_manager.findfont('Open Sans', fallback_to_default=False)).name == 'OpenSans-Regular.ttf'
mpl.rcParams.update({'font.family': 'Open Sans', 'font.size': 6,
    'axes.labelsize': 7, 'axes.titlesize': 7, 'xtick.labelsize': 6,
    'ytick.labelsize': 6, 'legend.fontsize': 6, 'legend.title_fontsize': 7,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.linewidth': 0.35, 'grid.linewidth': 0.35, 'grid.color': '#dbe4ef',
    'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'svg.hashsalt': 'combinatorics-figure3', 'savefig.facecolor': 'white'})
for folder in ['previews', 'evidence']:
    (OUTPUT_DIR / folder).mkdir(parents=True, exist_ok=True)

# %%
# Read exact stored values and preserve the declared scientific order.
tables = {}
metadata = {'a': {'keys': ['Gene', 'Protein']},
 'b': {'keys': ['Platform', 'Gene'],
       'gene_order': {'Long-reads': ['GRIA2', 'PCLO'],
                      'Short-reads': ['ABL',
                                      'ACHA4',
                                      'ANR17',
                                      'CA2D3',
                                      'CACB2',
                                      'CSKI1',
                                      'DGLA',
                                      'DOP1',
                                      'IQEC1',
                                      'K0513',
                                      'KCNAS',
                                      'MTUS2',
                                      'PCLO',
                                      'RIMS2',
                                      'ROBO2',
                                      'RUSC2',
                                      'SCN1',
                                      'TRIM2',
                                      'TWK7']},
       'gene_colors': {'GRIA2': '#00A08B',
                       'PCLO': '#AF0038',
                       'ABL': '#0D2A63',
                       'ACHA4': '#DA16FF',
                       'ANR17': '#FC0080',
                       'CA2D3': '#6C4516',
                       'CACB2': '#1616A7',
                       'CSKI1': '#222A2A',
                       'DGLA': '#620042',
                       'DOP1': '#A777F1',
                       'IQEC1': '#EB663B',
                       'K0513': '#750D86',
                       'KCNAS': '#FB0D0D',
                       'MTUS2': '#1CA71C',
                       'RIMS2': '#FB00D1',
                       'ROBO2': '#2E91E5',
                       'RUSC2': '#B68100',
                       'SCN1': '#6C7C32',
                       'TRIM2': '#B2828D',
                       'TWK7': '#862A16'}},
 'c': {'keys': ['Platform', 'Gene', 'Fraction'],
       'gene_order': {'Long-reads': ['GRIA2', 'PCLO'],
                      'Short-reads': ['ABL',
                                      'ACHA4',
                                      'ANR17',
                                      'CA2D3',
                                      'CACB2',
                                      'CSKI1',
                                      'DGLA',
                                      'DOP1',
                                      'IQEC1',
                                      'K0513',
                                      'KCNAS',
                                      'MTUS2',
                                      'PCLO',
                                      'RIMS2',
                                      'ROBO2',
                                      'RUSC2',
                                      'SCN1',
                                      'TRIM2',
                                      'TWK7']},
       'gene_colors': {'GRIA2': '#00A08B',
                       'PCLO': '#AF0038',
                       'ABL': '#0D2A63',
                       'ACHA4': '#DA16FF',
                       'ANR17': '#FC0080',
                       'CA2D3': '#6C4516',
                       'CACB2': '#1616A7',
                       'CSKI1': '#222A2A',
                       'DGLA': '#620042',
                       'DOP1': '#A777F1',
                       'IQEC1': '#EB663B',
                       'K0513': '#750D86',
                       'KCNAS': '#FB0D0D',
                       'MTUS2': '#1CA71C',
                       'RIMS2': '#FB00D1',
                       'ROBO2': '#2E91E5',
                       'RUSC2': '#B68100',
                       'SCN1': '#6C7C32',
                       'TRIM2': '#B2828D',
                       'TWK7': '#862A16'},
       'reference_line_coordinates': [{'name': 'Theoretical maximum    ',
                                       'x': [3571.7000000000003, 75005.7],
                                       'y': [36016, 36016],
                                       'platform': 'Long-reads'},
                                      {'name': 'Observed',
                                       'x': [3571.7000000000003, 75005.7],
                                       'y': [18496, 18496],
                                       'platform': 'Long-reads'},
                                      {'name': 'Theoretical maximum    ',
                                       'x': [72753.95, 1527832.95],
                                       'y': [36016, 36016],
                                       'platform': 'Short-reads'},
                                      {'name': 'Observed',
                                       'x': [72753.95, 1527832.95],
                                       'y': [18496, 18496],
                                       'platform': 'Short-reads'}]},
 'd': {'keys': ['Dissimilarity', 'Gene', 'Fraction']}}
for panel in 'abcd':
    path = INPUT_DIR / f'squid_original_figure3{panel}_plot.tsv'
    table = pd.read_table(path, keep_default_na=False, na_values=[''], float_precision='round_trip')
    keys = metadata[panel]['keys']
    assert not table.duplicated(keys).any() and table[keys].notna().all().all()
    tables[panel] = table
metadata['a']['codon_order'] = tables['a'].columns[4:].tolist()
matrix = tables['a'].loc[:, metadata['a']['codon_order']].to_numpy(dtype=float)
assert matrix.shape == (100, 94) and np.isnan(matrix).sum() == 623
assert set(np.unique(matrix[~np.isnan(matrix)])) == {0, 1}
assert tables['a']['IsoformRank'].tolist() == list(range(1, 101))
assert tables['a']['Gene'].eq('PCLO').all()
assert len(tables['b']) == 21 and len(tables['c']) == 105 and len(tables['d']) == 20
assert tables['d']['Dissimilarity'].drop_duplicates().tolist() == ['Miyata', 'Grantham 100']
gene_colors = metadata['b']['gene_colors']


# %%
# Axes rectangles are specified as x, top, width, height in millimetres.
def axes_mm(figure, x, top, width, height):
    return figure.add_axes([x / WIDTH_MM, 1 - (top + height) / HEIGHT_MM,
                            width / WIDTH_MM, height / HEIGHT_MM])

def style_axes(axis, hide_y=False):
    axis.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
    axis.grid(True, axis='y', zorder=0)
    axis.tick_params(pad=1.5, labelleft=not hide_y)
    axis.tick_params(axis='y', which='both', left=False, right=False)
    axis.tick_params(axis='x', which='both', bottom=False, top=False)
    axis.set_axisbelow(True)
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f'{value / 1000:g}k' if value >= 1000 else f'{value:g}'))

# Every text role is retained for a post-draw point-size and bounds audit.
role_text = []
def figure_text(figure, x, top, text, role, **kwargs):
    artist = figure.text(x / WIDTH_MM, 1 - top / HEIGHT_MM, text,
                        fontsize=TEXT_PT[role], **kwargs)
    artist.set_gid(f'role-{role}-{len(role_text)}')
    role_text.append((role, artist))
    return artist


# %%
# Panel a: native vector cells; scientific NA is separate from display encoding.
def draw_panel_a(figure):
    axis = axes_mm(figure, 13, 10, 32, 39)
    axis.set_facecolor('#000000')  # Prevent white seams between vector cells.
    states = np.where(np.isnan(matrix), -1, matrix)
    colors = ListedColormap(['#c0c0c0', '#000000', '#4af626'])
    mesh = axis.pcolormesh(states, cmap=colors,
        norm=BoundaryNorm([-1.5, -0.5, 0.5, 1.5], 3),
        edgecolors='face', linewidth=0.03, antialiased=False, rasterized=False)
    np.testing.assert_array_equal(np.asarray(mesh.get_array()), states)
    axis.set_ylim(100, 0)
    axis.set_xticks([])
    axis.set_yticks([19.5, 39.5, 59.5, 79.5, 99.5], [20, 40, 60, 80, 100])
    axis.set_ylabel('Most expressed PCLO isoforms', labelpad=2)
    axis.set_xlabel('Editable amino acids', labelpad=3)
    axis.tick_params(axis='y', pad=1, length=1.5)
    for top, color, label in [(24.5, '#4af626', 'Recoded'),
                               (30, '#000000', 'Not-recoded'),
                               (35.5, '#c0c0c0', 'NA')]:
        figure.add_artist(Rectangle((48 / WIDTH_MM, 1 - (top + 2.8) / HEIGHT_MM),
             2.8 / WIDTH_MM, 2.8 / HEIGHT_MM, transform=figure.transFigure,
             facecolor=color, edgecolor='none'))
        figure_text(figure, 51.5, top + 1.4, label, 'legend', ha='left', va='center')
    figure_text(figure, 3, 7, 'a', 'panel', weight='bold')
    return axis


# %%
# Panel b: one maximal isoform count per gene and sequencing platform.
def draw_panel_b(figure):
    axes = []
    for platform, x, width in [('Long-reads', 80, 12), ('Short-reads', 96, 81)]:
        axis = axes_mm(figure, x, 13, width, 32)
        table = tables['b'].loc[tables['b']['Platform'].eq(platform)]
        assert table['Gene'].tolist() == metadata['b']['gene_order'][platform]
        bars = axis.bar(range(len(table)), table['Max'],
                        color=[gene_colors[gene] for gene in table['Gene']], width=0.78, zorder=3)
        np.testing.assert_array_equal([bar.get_height() for bar in bars], table['Max'])
        axis.set_xticks(range(len(table)), table['Gene'], rotation=60, ha='right', rotation_mode='anchor')
        axis.set_xlim(-0.65, len(table) - 0.3)
        axis.set_ylim(0, 70000)
        axis.yaxis.set_major_locator(MultipleLocator(10000))
        axis.set_title(platform, pad=4)
        style_axes(axis, hide_y=platform == 'Short-reads')
        axes.append(axis)
    axes[0].set_ylabel('Max distinct protein\nisoforms observed', labelpad=2)
    figure_text(figure, 129, 57, 'Gene', 'axis', ha='center')
    figure_text(figure, 69, 7, 'b', 'panel', weight='bold')
    return axes


# %%
# Panel c: keep maximal-per-fraction curves and original DSCAM line endpoints.
def draw_panel_c(figure):
    axes = []
    for platform, x, width in [('Long-reads', 14, 32), ('Short-reads', 50, 77)]:
        axis = axes_mm(figure, x, 70, width, 55)
        table = tables['c'].loc[tables['c']['Platform'].eq(platform)]
        for gene in metadata['c']['gene_order'][platform]:
            curve = table.loc[table['Gene'].eq(gene)]
            assert curve['Fraction'].is_monotonic_decreasing and len(curve) == 5
            line, = axis.plot(curve['NumOfReads'], curve['NumOfProteins'],
                color=gene_colors[gene], marker='D', markersize=1.5, linewidth=0.6)
            np.testing.assert_array_equal(line.get_xydata(), curve[['NumOfReads', 'NumOfProteins']].to_numpy())
        references = [line for line in metadata['c']['reference_line_coordinates'] if line['platform'] == platform]
        for reference, color in zip(references, ['#aaaaaa', '#686868']):
            axis.plot(reference['x'], reference['y'], color=color, linewidth=0.9, dashes=(3, 2))
        axis.set_xlim(0, table['NumOfReads'].max() * 1.08)
        axis.set_ylim(0, 71000)
        axis.yaxis.set_major_locator(MultipleLocator(10000))
        axis.xaxis.set_major_locator(MultipleLocator(20000 if platform == 'Long-reads' else 200000))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, pos, platform=platform: '0' if value == 0 else (f'{value/1000:g}k' if platform == 'Long-reads' else f'{value/1e6:g}M')))
        axis.tick_params(axis='x', labelrotation=25)
        axis.set_title(platform, pad=4)
        style_axes(axis, hide_y=platform == 'Short-reads')
        axis.grid(True, axis='x')
        axes.append(axis)
    axes[0].set_ylabel('Distinct protein isoforms observed', labelpad=3)
    figure_text(figure, 70.5, 134, 'Coverage', 'axis', ha='center')
    figure_text(figure, 3, 63, 'c', 'panel', weight='bold')
    return axes


# %%
# Legend groups occupy their own area, outside the data axes.
def draw_panel_c_legend(figure):
    top = 68.0
    for platform in ['Long-reads', 'Short-reads']:
        figure_text(figure, 132, top, platform, 'legend_heading')
        top += 3.2
        for gene in metadata['c']['gene_order'][platform]:
            figure.add_artist(Line2D(np.array([132, 133.5, 135]) / WIDTH_MM,
                [1 - top / HEIGHT_MM] * 3, transform=figure.transFigure,
                color=gene_colors[gene], linewidth=0.6, marker='D', markersize=2, markevery=[1]))
            figure_text(figure, 136.5, top, gene, 'legend', va='center')
            top += 2.35
        top += 3.4  # Separate each heading from the preceding group's items.
    # DSCAM occupies the upper-right column beside the Long-reads legend.
    top = 68.0
    figure_text(figure, 151, top, 'Drosophila’s DSCAM', 'legend_heading')
    top += 3.3
    for label, color in [('Theoretical maximum', '#aaaaaa'), ('Observed', '#686868')]:
        figure.add_artist(Line2D([151 / WIDTH_MM, 154 / WIDTH_MM],
            [1 - top / HEIGHT_MM] * 2, transform=figure.transFigure,
            color=color, linewidth=0.9, dashes=(3, 2)))
        figure_text(figure, 155.5, top, label, 'legend', va='center')
        top += 2.7


# %%
# Panel d: original Miyata / Grantham 100 maximal curves and separate y axes.
def draw_panel_d(figure):
    axes = []
    for category, x, ymax, step, title in [
        ('Miyata', 20, 4300, 600, 'Miyata'),
        ('Grantham 100', 61.4, 360, 50, 'Grantham 100')]:
        axis = axes_mm(figure, x, 149, 31.8, 22)
        table = tables['d'].loc[tables['d']['Dissimilarity'].eq(category)]
        for gene in ['GRIA2', 'PCLO']:
            curve = table.loc[table['Gene'].eq(gene)]
            line, = axis.plot(curve['NumOfReads'], curve['NumOfProteins'],
                color=gene_colors[gene], marker='D', markersize=1.5, linewidth=0.6)
            np.testing.assert_array_equal(line.get_xydata(), curve[['NumOfReads', 'NumOfProteins']].to_numpy())
        axis.set_xlim(0, tables['d']['NumOfReads'].max() * 1.05)
        axis.set_ylim(0, ymax)
        axis.yaxis.set_major_locator(MultipleLocator(step))
        axis.xaxis.set_major_locator(MultipleLocator(20000))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: '0' if value == 0 else f'{value/1000:g}k'))
        axis.set_title(title, pad=4)
        style_axes(axis)
        axis.grid(True, axis='x')
        axes.append(axis)
    axes[0].set_ylabel('Dissimilar protein\nisoforms observed', labelpad=3)
    figure_text(figure, 20, 142, 'Long-reads', 'heading')
    figure_text(figure, 56.6, 178, 'Coverage', 'axis', ha='center')
    handles = [Line2D([], [], color=gene_colors[g], marker='D', markersize=2, linewidth=0.6, label=g) for g in ['GRIA2', 'PCLO']]
    figure.legend(handles=handles, loc='upper left', bbox_to_anchor=(100/WIDTH_MM, 1-149/HEIGHT_MM), borderaxespad=0, frameon=False, handlelength=1.2)
    figure_text(figure, 3, 142, 'd', 'panel', weight='bold')
    return axes


# %%
if RENDER_FIGURE3:
    # Compose directly at the final physical size; no panel scaling after export.
    role_text.clear()
    figure = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=150)
    panel_axes = {'a': [draw_panel_a(figure)], 'b': draw_panel_b(figure),
                  'c': draw_panel_c(figure), 'd': draw_panel_d(figure)}
    draw_panel_c_legend(figure)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    # Axis titles/headings use 7 pt; every tick and legend entry uses 6 pt.
    for axes in panel_axes.values():
        for axis in axes:
            for text in [axis.xaxis.label, axis.yaxis.label, axis.title]:
                if text.get_text():
                    assert text.get_fontsize() == 7
            for text in axis.get_xticklabels() + axis.get_yticklabels():
                assert text.get_fontsize() == 6
    for role, text in role_text:
        assert text.get_fontsize() == TEXT_PT[role]
    print('Figure 3 scientific artists and requested font roles checked.')

# %%
if RENDER_FIGURE3:
    # Export the native full canvas; the selected-page cell applies the final crop.
    figure3_raw = FINAL_FIGURES_DIR / 'render_work/raw'
    for folder in ['svg', 'pdf', 'previews']:
        (figure3_raw / folder).mkdir(parents=True, exist_ok=True)
    for extension in ['svg', 'pdf']:
        figure.savefig(figure3_raw / extension / f'Figure3.{extension}', format=extension,
                       metadata={'Creator': 'Combinatorics — reproducible Figure 3'})
    figure.savefig(figure3_raw / 'previews/Figure3.png', dpi=300)
    figure.savefig(OUTPUT_DIR / 'previews/Figure3_print_scale.png', dpi=96)
    # These panel crops are diagnostic previews; the final physical crop is separate.
    from matplotlib.transforms import Bbox
    preview_regions_mm = {'a': (0, 4, 64.5, 58), 'b': (64.5, 4, 180, 60),
                          'c': (0, 59, 180, 139), 'd': (0, 138, 180, 180)}
    for panel, (left, top, right, bottom) in preview_regions_mm.items():
        bounds = Bbox.from_extents(left / 25.4, (HEIGHT_MM-bottom) / 25.4,
                                  right / 25.4, (HEIGHT_MM-top) / 25.4)
        figure.savefig(OUTPUT_DIR / f'previews/Figure3{panel}.png', dpi=300, bbox_inches=bounds)

# %%
if RENDER_FIGURE3:
    # Save measured text bounds and sizes; external SVG/PDF inspection follows.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    text_audit = []
    for text in figure.findobj(Text):
        if text.get_visible() and text.get_text():
            bounds = text.get_window_extent(renderer)
            text_audit.append({'text': text.get_text(), 'size_pt': text.get_fontsize(),
                'font_family': text.get_fontfamily(), 'weight': text.get_fontweight(),
                'bounds_mm': [value / figure.dpi * 25.4 for value in bounds.extents],
                'gid': text.get_gid(), 'rotation': text.get_rotation()})
    canvas_mm = [WIDTH_MM, HEIGHT_MM]
    axes_audit = []
    for panel, axes in panel_axes.items():
        for axis in axes:
            x, y, width, height = axis.get_position().bounds
            axes_audit.append({'panel': panel, 'x_mm': x * WIDTH_MM,
                'top_mm': (1-y-height) * HEIGHT_MM, 'width_mm': width * WIDTH_MM,
                'height_mm': height * HEIGHT_MM,
                'visible_spines': [name for name, spine in axis.spines.items() if spine.get_visible()],
                'visible_y_tick_marks': sum(line.get_visible() for tick in axis.yaxis.get_major_ticks()
                                            for line in [tick.tick1line, tick.tick2line])})
            if panel != 'a':
                assert not any(spine.get_visible() for spine in axis.spines.values())
                assert axes_audit[-1]['visible_y_tick_marks'] == 0
    clipped = [item for item in text_audit if item['bounds_mm'][0] < 0 or item['bounds_mm'][1] < 0
               or item['bounds_mm'][2] > WIDTH_MM or item['bounds_mm'][3] > HEIGHT_MM]
    (OUTPUT_DIR / 'evidence/figure3_artist_audit.json').write_text(json.dumps(
        {'canvas_mm': canvas_mm, 'axes': axes_audit, 'text': text_audit, 'outside_canvas': clipped,
         'matrix_cells': int(matrix.size), 'matrix_na': int(np.isnan(matrix).sum()),
         'scientific_artist_arrays_exact': True, 'raw_analyses_executed': False}, indent=2))
    assert not clipped, [(item['text'], item['bounds_mm']) for item in clipped]
    print('No visible text extends beyond the final canvas.')
    display(Image(filename=str(OUTPUT_DIR / 'previews/Figure3_print_scale.png')))

# %% [markdown]
# # Supplemental figure S15
# Continue after the Figure 3 cells above. Figure 3 **revision8** is the approved
# baseline, including the absence of X/Y tick marks in panels b–d. Its drawing
# cells and input tables are unchanged by this extension.
#
# Only verified Phase 1 TSVs are loaded. No source notebook, clustering, raw
# positions, or upstream analysis is executed. The saved source traces and the
# designated reference SVGs were checked separately; their small layout manifest
# is defined explicitly in the following plotting cells.
#
# The canonical S7 producer is at the end of this notebook.
# S15 uses a 180 × 250 mm working canvas. The selected-page export applies the approved final crop at physical scale one; temporary canvases and evidence stay under the new `FINAL_FIGURES_DIR/render_work` directory.
#
# Caption follow-up: Figure 3a shows **94 positions**, not 96. Counts of observed
# amino-acid alternatives do not necessarily equal theoretical possibilities.
# Manuscript integration and caption edits are deferred to the current Word
# manuscript, using Track Changes and preserving the agreed display dimensions.

# %%
from matplotlib.transforms import Bbox
# Physical placement helpers use each supplemental figure's own canvas size.
SUPPLEMENT_DIR = FINAL_FIGURES_DIR / 'render_work/S15'
reference_layout = {'s15': {'rows': [{'panel': 'a',
                   'figure': 1,
                   'eligibility': 'All genes',
                   'ylim': [0, 22.307084339146076],
                   'dtick': 5,
                   'labels': {'SquidLongReads': {},
                              'SquidShortReads': {},
                              'Octopus': {'A>C': {'inside': True, 'color': '#ffffff'},
                                          'A>G': {'inside': True, 'color': '#444444'},
                                          'A>T': {'inside': True, 'color': '#ffffff'},
                                          'C>A': {'inside': True, 'color': '#ffffff'},
                                          'C>G': {'inside': True, 'color': '#ffffff'},
                                          'C>T': {'inside': True, 'color': '#ffffff'},
                                          'G>A': {'inside': True, 'color': '#444444'},
                                          'G>C': {'inside': True, 'color': '#ffffff'},
                                          'G>T': {'inside': True, 'color': '#444444'},
                                          'T>A': {'inside': True, 'color': '#ffffff'},
                                          'T>C': {'inside': True, 'color': '#ffffff'},
                                          'T>G': {'inside': True, 'color': '#ffffff'}}}},
                  {'panel': 'b_top',
                   'figure': 2,
                   'eligibility': 'Passed',
                   'ylim': [0, 105],
                   'dtick': 10,
                   'labels': {'SquidLongReads': {'A>G': {'inside': True, 'color': '#444444'}},
                              'SquidShortReads': {'A>G': {'inside': True, 'color': '#444444'},
                                                  'C>T': {'inside': False, 'color': '#2a3f5f'},
                                                  'T>C': {'inside': False, 'color': '#2a3f5f'}},
                              'Octopus': {'A>C': {'inside': False, 'color': '#2a3f5f'},
                                          'A>G': {'inside': True, 'color': '#444444'},
                                          'A>T': {'inside': False, 'color': '#2a3f5f'},
                                          'C>A': {'inside': False, 'color': '#2a3f5f'},
                                          'C>G': {'inside': False, 'color': '#2a3f5f'},
                                          'C>T': {'inside': False, 'color': '#2a3f5f'},
                                          'G>A': {'inside': False, 'color': '#2a3f5f'},
                                          'G>C': {'inside': False, 'color': '#2a3f5f'},
                                          'G>T': {'inside': False, 'color': '#2a3f5f'},
                                          'T>A': {'inside': False, 'color': '#2a3f5f'},
                                          'T>C': {'inside': False, 'color': '#2a3f5f'},
                                          'T>G': {'inside': False, 'color': '#2a3f5f'}}}},
                  {'panel': 'b_bottom',
                   'figure': 2,
                   'eligibility': 'Excluded',
                   'ylim': [0, 36.660591492364176],
                   'dtick': 10,
                   'labels': {'SquidLongReads': {},
                              'SquidShortReads': {},
                              'Octopus': {'A>C': {'inside': True, 'color': '#ffffff'},
                                          'A>G': {'inside': True, 'color': '#444444'},
                                          'A>T': {'inside': True, 'color': '#ffffff'},
                                          'C>A': {'inside': True, 'color': '#ffffff'},
                                          'C>G': {'inside': False, 'color': '#2a3f5f'},
                                          'C>T': {'inside': True, 'color': '#ffffff'},
                                          'G>A': {'inside': True, 'color': '#444444'},
                                          'G>C': {'inside': False, 'color': '#2a3f5f'},
                                          'G>T': {'inside': False, 'color': '#2a3f5f'},
                                          'T>A': {'inside': True, 'color': '#ffffff'},
                                          'T>C': {'inside': True, 'color': '#ffffff'},
                                          'T>G': {'inside': False, 'color': '#2a3f5f'}}}},
                  {'panel': 'c',
                   'figure': 3,
                   'eligibility': 'Passed',
                   'ylim': [0, 79.55749128919861],
                   'dtick': 10,
                   'labels': {'SquidLongReads': {'A>C': {'inside': False, 'color': '#2a3f5f'},
                                                 'A>G': {'inside': True, 'color': '#444444'},
                                                 'A>T': {'inside': False, 'color': '#2a3f5f'},
                                                 'C>A': {'inside': False, 'color': '#2a3f5f'},
                                                 'C>G': {'inside': False, 'color': '#2a3f5f'},
                                                 'C>T': {'inside': True, 'color': '#ffffff'},
                                                 'G>A': {'inside': True, 'color': '#444444'},
                                                 'G>C': {'inside': False, 'color': '#2a3f5f'},
                                                 'G>T': {'inside': False, 'color': '#2a3f5f'},
                                                 'T>A': {'inside': False, 'color': '#2a3f5f'},
                                                 'T>C': {'inside': True, 'color': '#ffffff'},
                                                 'T>G': {'inside': False, 'color': '#2a3f5f'}},
                              'SquidShortReads': {'A>C': {'inside': False, 'color': '#2a3f5f'},
                                                  'A>G': {'inside': True, 'color': '#444444'},
                                                  'A>T': {'inside': False, 'color': '#2a3f5f'},
                                                  'C>A': {'inside': True, 'color': '#ffffff'},
                                                  'C>G': {'inside': False, 'color': '#2a3f5f'},
                                                  'C>T': {'inside': True, 'color': '#ffffff'},
                                                  'G>A': {'inside': True, 'color': '#444444'},
                                                  'G>C': {'inside': False, 'color': '#2a3f5f'},
                                                  'G>T': {'inside': True, 'color': '#444444'},
                                                  'T>A': {'inside': False, 'color': '#2a3f5f'},
                                                  'T>C': {'inside': True, 'color': '#ffffff'},
                                                  'T>G': {'inside': False, 'color': '#2a3f5f'}},
                              'Octopus': {'A>C': {'inside': False, 'color': '#2a3f5f'},
                                          'A>G': {'inside': True, 'color': '#444444'},
                                          'A>T': {'inside': False, 'color': '#2a3f5f'},
                                          'C>A': {'inside': False, 'color': '#2a3f5f'},
                                          'C>G': {'inside': False, 'color': '#2a3f5f'},
                                          'C>T': {'inside': False, 'color': '#2a3f5f'},
                                          'G>A': {'inside': False, 'color': '#2a3f5f'},
                                          'G>C': {'inside': False, 'color': '#2a3f5f'},
                                          'G>T': {'inside': False, 'color': '#2a3f5f'},
                                          'T>A': {'inside': False, 'color': '#2a3f5f'},
                                          'T>C': {'inside': False, 'color': '#2a3f5f'},
                                          'T>G': {'inside': False, 'color': '#2a3f5f'}}}},
                  {'panel': 'd',
                   'figure': 4,
                   'eligibility': 'Passed',
                   'ylim': [0, 105],
                   'dtick': 10,
                   'labels': {'SquidLongReads': {'A>G': {'inside': True, 'color': '#444444'}},
                              'SquidShortReads': {'A>C': {'inside': False, 'color': '#2a3f5f'},
                                                  'A>G': {'inside': True, 'color': '#444444'},
                                                  'C>T': {'inside': False, 'color': '#2a3f5f'},
                                                  'G>A': {'inside': False, 'color': '#2a3f5f'},
                                                  'G>T': {'inside': False, 'color': '#2a3f5f'},
                                                  'T>C': {'inside': False, 'color': '#2a3f5f'}},
                              'Octopus': {'A>C': {'inside': False, 'color': '#2a3f5f'},
                                          'A>G': {'inside': True, 'color': '#444444'},
                                          'A>T': {'inside': False, 'color': '#2a3f5f'},
                                          'C>A': {'inside': False, 'color': '#2a3f5f'},
                                          'C>G': {'inside': False, 'color': '#2a3f5f'},
                                          'C>T': {'inside': False, 'color': '#2a3f5f'},
                                          'G>A': {'inside': False, 'color': '#2a3f5f'},
                                          'G>C': {'inside': False, 'color': '#2a3f5f'},
                                          'G>T': {'inside': False, 'color': '#2a3f5f'},
                                          'T>A': {'inside': False, 'color': '#2a3f5f'},
                                          'T>C': {'inside': False, 'color': '#2a3f5f'},
                                          'T>G': {'inside': False, 'color': '#2a3f5f'}}}}]}}
for folder in ['previews', 'evidence']:
    (SUPPLEMENT_DIR / folder).mkdir(parents=True, exist_ok=True)
REFERENCE_INK = '#2a3f5f'

def supplement_axes(figure, x, top, width, height):
    canvas_width, canvas_height = figure.get_size_inches() * 25.4
    return figure.add_axes([x / canvas_width, 1 - (top + height) / canvas_height,
                            width / canvas_width, height / canvas_height])

def supplement_text(figure, x, top, label, role, **kwargs):
    canvas_width, canvas_height = figure.get_size_inches() * 25.4
    artist = figure.text(x / canvas_width, 1 - top / canvas_height, label,
        fontsize=TEXT_PT[role], color=REFERENCE_INK, **kwargs)
    artist.set_gid(f'role-{role}-{len(figure.texts)}')
    return artist

def supplement_style(axis, hide_y=False):
    style_axes(axis, hide_y=hide_y)
    axis.tick_params(colors=REFERENCE_INK)
    for text in [axis.title, axis.xaxis.label, axis.yaxis.label]:
        text.set_color(REFERENCE_INK)


# %%
# Audit full-canvas text, final point sizes, and the absence of black frames/ticks.
def audit_supplement(figure, name):
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    width, height = figure.get_size_inches() * 25.4
    text_records, axes_records = [], []
    for artist in figure.findobj(Text):
        if not artist.get_visible() or not artist.get_text():
            continue
        bounds = artist.get_window_extent(renderer).extents / figure.dpi * 25.4
        assert artist.get_fontsize() in [6, 7]
        assert artist.get_fontfamily() == ['Open Sans']
        text_records.append({'text': artist.get_text(), 'size_pt': artist.get_fontsize(),
            'bounds_mm': bounds.tolist(), 'gid': artist.get_gid(),
            'rotation': artist.get_rotation(), 'weight': artist.get_fontweight()})
    for axis in figure.axes:
        assert not any(spine.get_visible() for spine in axis.spines.values())
        ticks = axis.xaxis.get_major_ticks() + axis.yaxis.get_major_ticks()
        assert not any(tick.tick1line.get_visible() or tick.tick2line.get_visible() for tick in ticks)
        axes_records.append({'bounds_fraction': list(axis.get_position().bounds),
            'xlim': list(axis.get_xlim()), 'ylim': list(axis.get_ylim()),
            'lines': len(axis.lines), 'collections': len(axis.collections),
            'bars': len(axis.patches), 'visible_spines': [], 'visible_tick_marks': 0})
    clipped = [t for t in text_records if t['bounds_mm'][0] < 0 or t['bounds_mm'][1] < 0
               or t['bounds_mm'][2] > width or t['bounds_mm'][3] > height]
    record = {'canvas_mm': [width, height], 'text': text_records, 'axes': axes_records,
              'outside_canvas': clipped, 'scientific_artist_arrays_exact': True}
    (SUPPLEMENT_DIR / f'evidence/{name}_artist_audit.json').write_text(json.dumps(record, indent=2))
    assert not clipped, [(t['text'], t['bounds_mm']) for t in clipped]
    return record


# %%
# Vector exports retain the full canvas. Cropping applies only to panel previews.
def export_supplement(figure, name, panel_regions):
    audit_supplement(figure, name)
    for extension in ['svg', 'pdf']:
        figure.savefig(FINAL_FIGURES_DIR / extension / f'{name}.{extension}', format=extension,
            metadata={'Creator': 'Combinatorics — verified supplemental figures'})
    figure.savefig(FINAL_FIGURES_DIR / f'previews/{name}.png', dpi=300)
    figure.savefig(SUPPLEMENT_DIR / f'previews/{name}_print_scale.png', dpi=96)
    height = figure.get_size_inches()[1] * 25.4
    for panel, (left, top, right, bottom) in panel_regions.items():
        bounds = Bbox.from_extents(left / 25.4, (height - bottom) / 25.4,
                                  right / 25.4, (height - top) / 25.4)
        figure.savefig(SUPPLEMENT_DIR / f'previews/{name}{panel}.png', dpi=300, bbox_inches=bounds)
    print(f'{name}: vector SVG/PDF exported at {figure.get_size_inches() * 25.4} mm.')
    display(Image(filename=str(SUPPLEMENT_DIR / f'previews/{name}_print_scale.png')))


# %% [markdown]
# ## Figure S15 — combined mismatch-frequency panels
# The reference is `Combined.12mm.Figures1-4ForS14.svg`; the final figure number
# is **S15**. Sources: `various_plots.BQ30.ipynb`, cells `af126453` and `1fb847ad`.
# The layout contains a, b (two rows), c and d, each with the original three
# dataset columns and 12 mismatch categories. Family 5 belongs to S2 and is
# excluded. The four zero-denominator facets keep their axes and missing
# percentages, with no invented bars. Percentages retain the separate
# Figure/Dataset/Eligibility denominators. Color, category order, axis ranges,
# and one-decimal annotation text follow the combined reference.

# %%
# Read prepared counts only; sample counts are used to verify aggregation.
s15_meta = {'dataset_order': ['SquidLongReads', 'SquidShortReads', 'Octopus'],
 'dataset_labels': {'SquidLongReads': "Squid's Long-reads",
                    'SquidShortReads': "Squid's Short-reads",
                    'Octopus': 'Whole-transcriptome octopus data'},
 'mismatch_order': ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G'],
 'mismatch_colors': {'A>C': '#2E91E5',
                     'A>G': '#E15F99',
                     'A>T': '#1CA71C',
                     'C>A': '#FB0D0D',
                     'C>G': '#DA16FF',
                     'C>T': '#222A2A',
                     'G>A': '#B68100',
                     'G>C': '#750D86',
                     'G>T': '#EB663B',
                     'T>A': '#511CFB',
                     'T>C': '#00A08B',
                     'T>G': '#FB00D1'},
 'figure_titles': {'1': 'Mismatches at 100% mismatch frequency',
                   '2': 'Mismatches at or above suspected SNP threshold, but below 100% mismatch frequency',
                   '3': 'Mismatches at or below editing threshold',
                   '4': 'Mismatches above editing threshold, but below suspected SNP threshold',
                   '5': 'Mismatches above editing threshold, but below 100% mismatch frequency'},
 'subtitle': 'Only genes with at most 3 suspected SNPs are included (editing detection enabled)',
 'panel_keys': ['Figure', 'Dataset', 'Eligibility', 'Mismatch'],
 'sample_keys': ['Figure', 'Dataset', 'Platform', 'Sample', 'Chrom', 'Eligibility', 'Mismatch'],
 'columns': {'mm12_panel_counts': ['Dataset',
                                   'Mismatch',
                                   'Count',
                                   'PanelTotal',
                                   'Percentage',
                                   'Figure',
                                   'Eligibility'],
             'mm12_sample_counts': ['Dataset',
                                    'Platform',
                                    'Sample',
                                    'Chrom',
                                    'Mismatch',
                                    'Count',
                                    'Figure',
                                    'Eligibility']}}
s15_meta['checksums'] = {}
s15_input_tables = {}
for name in ['panel', 'sample']:
    path = figure_source(f'combined_s2_s15_mm12_{name}_counts')
    s15_meta['checksums'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    table = pd.read_table(path, float_precision='round_trip')
    schema_name = f'mm12_{name}_counts'
    assert table.columns.tolist() == s15_meta['columns'][schema_name]
    assert not table.duplicated(s15_meta[f'{name}_keys']).any()
    s15_input_tables[name] = table
s15_panels = s15_input_tables['panel'].loc[s15_input_tables['panel']['Figure'].between(1, 4)]
s15_samples = s15_input_tables['sample'].loc[s15_input_tables['sample']['Figure'].between(1, 4)]
assert len(s15_panels) == 180
count_keys = s15_meta['panel_keys']
grouped_counts = s15_samples.groupby(count_keys, sort=False)['Count'].sum().sort_index()
panel_counts = s15_panels.set_index(count_keys)['Count'].sort_index()
# Sample counts are sparse: absent categories must have exactly zero panel counts.
assert not s15_samples[count_keys].isna().any().any()
assert grouped_counts.index.difference(panel_counts.index).empty
absent_categories = panel_counts.index.difference(grouped_counts.index)
assert panel_counts.loc[absent_categories].eq(0).all()
pd.testing.assert_series_equal(grouped_counts, panel_counts.loc[grouped_counts.index])

# %%
# Zero totals remain NA; all other percentages use the original facet total.
s15_facet_checks = []
for keys, facet in s15_panels.groupby(['Figure', 'Dataset', 'Eligibility'], sort=False):
    assert facet['Mismatch'].tolist() == s15_meta['mismatch_order']
    total = int(facet['Count'].sum())
    assert facet['PanelTotal'].eq(total).all()
    if total == 0:
        assert facet['Percentage'].isna().all()
    else:
        np.testing.assert_allclose(facet['Percentage'], facet['Count'] / total * 100,
                                   rtol=1e-14, atol=1e-14)
    s15_facet_checks.append({'figure': int(keys[0]), 'dataset': keys[1],
        'eligibility': keys[2], 'total': total, 'missing_percentages': int(facet['Percentage'].isna().sum())})
assert len(s15_facet_checks) == 15
assert sum(item['total'] == 0 for item in s15_facet_checks) == 4
s15_figure = plt.figure(figsize=(180 / 25.4, 250 / 25.4), dpi=150)
s15_artist_checks, s15_label_adjustments = [], []


# %%
# Retain the reference labels; move an inside label out only if 6 pt cannot fit.
def label_s15_bars(axis, bars, values, row, dataset):
    label_style = row['labels'][dataset]
    renderer = s15_figure.canvas.get_renderer()
    font = font_manager.FontProperties(family='Open Sans', size=6)
    height_mm = axis.get_position().height * 250
    ymax = row['ylim'][1]
    for mismatch, bar, value in zip(s15_meta['mismatch_order'], bars, values):
        if value == 0:
            assert mismatch not in label_style
            continue
        spec = label_style[mismatch]
        label = f'{value:.1f}'
        text_width_px = renderer.get_text_width_height_descent(label, font, False)[0]
        text_height_mm = text_width_px / s15_figure.dpi * 25.4
        inside = spec['inside'] and value / ymax * height_mm > text_height_mm + 0.6
        if inside != spec['inside']:
            s15_label_adjustments.append({'row': row['panel'], 'dataset': dataset, 'mismatch': mismatch})
        offset = 0.3 / height_mm * ymax
        text = axis.text(bar.get_x() + bar.get_width() / 2, value + (-offset if inside else offset),
            label, fontsize=6, rotation=270 if inside else 0,
            ha='center', va='top' if inside else 'bottom',
            color=spec['color'] if inside else REFERENCE_INK, clip_on=False)
        text.set_gid(f'S15-value-{row["panel"]}-{dataset}-{mismatch.replace(">", "to")}')


# %%
# Each row has the same category/color order, with separate original denominators.
def draw_s15_row(row, top):
    for column, (dataset, x) in enumerate(zip(s15_meta['dataset_order'], [16, 67, 118])):
        axis = supplement_axes(s15_figure, x, top, 44, S15_ROW_HEIGHT_MM)
        facet = s15_panels.loc[s15_panels['Figure'].eq(row['figure'])
            & s15_panels['Dataset'].eq(dataset) & s15_panels['Eligibility'].eq(row['eligibility'])]
        values = facet['Percentage'].to_numpy()
        assert facet['Mismatch'].tolist() == s15_meta['mismatch_order']
        axis.set(xlim=(-0.5, 11.5), ylim=row['ylim'])
        axis.set_xticks(range(12), s15_meta['mismatch_order'],
                       rotation=-45, ha='left', rotation_mode='anchor')
        axis.set_yticks(np.arange(0, row['ylim'][1], row['dtick']))
        if column == 0:
            axis.set_ylabel('Percentage', labelpad=3)
        supplement_style(axis, hide_y=column != 0)
        if np.isnan(values).all():
            s15_artist_checks.append({'row': row['panel'], 'dataset': dataset, 'bars': 0, 'missing': 12})
            continue
        assert np.isfinite(values).all()
        bars = axis.bar(range(12), values, width=0.8, zorder=3,
            color=[s15_meta['mismatch_colors'][m] for m in s15_meta['mismatch_order']], linewidth=0)
        for index, bar in enumerate(bars):
            bar.set_gid(f'S15-bar-{row["panel"]}-{dataset}-{index}')
        np.testing.assert_array_equal([bar.get_height() for bar in bars], values)
        label_s15_bars(axis, bars, values, row, dataset)
        s15_artist_checks.append({'row': row['panel'], 'dataset': dataset, 'bars': len(bars), 'missing': 0})


# %%
# Reduce plotting heights only; preserve the 8 mm group and 10 mm internal gaps.
S15_ROW_HEIGHT_MM = 28.5  # 5% less than the approved 30 mm rows.
S15_INTERNAL_ROW_GAP_MM = 10.0
S15_GROUP_GAP_MM = 8.0
S15_COLUMN_HEADING_GAP_MM = 6.0  # Column headings 2 mm closer to panel a.
S15_COLUMN_HEADING_INK_BOTTOM_MM = 7.60
S15_TITLE_INK_ASCENT_MM, S15_X_TICK_INK_DEPTH_MM = 1.90, 5.71
s15_panel_tops = []
previous_content_bottom = S15_COLUMN_HEADING_INK_BOTTOM_MM
for family in range(1, 5):
    heading_to_plot = 1.7 if family in [1, 2] else 5.7  # c/d include a subtitle.
    gap = S15_COLUMN_HEADING_GAP_MM if family == 1 else S15_GROUP_GAP_MM
    title_baseline = previous_content_bottom + gap + S15_TITLE_INK_ASCENT_MM
    plot_top = title_baseline + heading_to_plot
    s15_panel_tops.append(plot_top)
    panel_height = (2 * S15_ROW_HEIGHT_MM + S15_INTERNAL_ROW_GAP_MM
                    if family == 2 else S15_ROW_HEIGHT_MM)
    previous_content_bottom = plot_top + panel_height + S15_X_TICK_INK_DEPTH_MM
s15_row_tops = [s15_panel_tops[0], s15_panel_tops[1],
                s15_panel_tops[1] + S15_ROW_HEIGHT_MM + S15_INTERNAL_ROW_GAP_MM,
                s15_panel_tops[2], s15_panel_tops[3]]


# %%
if 'FigureS15' in FIGURES_TO_BUILD:
    # Five plot rows form four panels; headings and row labels remain outside the axes.
    for row, top in zip(reference_layout['s15']['rows'], s15_row_tops):
        draw_s15_row(row, top)
    for dataset, center in zip(s15_meta['dataset_order'], [38, 89, 140]):
        supplement_text(s15_figure, center, 7, s15_meta['dataset_labels'][dataset], 'heading', ha='center')
    # Anchor headings to their own plot; leave about 1 mm below the last heading line.
    for panel, family, plot_top in zip('abcd', range(1, 5), s15_panel_tops):
        top = plot_top - (5.7 if panel in ['c', 'd'] else 1.7)
        supplement_text(s15_figure, 3, top, panel, 'panel', weight='bold')
        supplement_text(s15_figure, 89, top, s15_meta['figure_titles'][str(family)],
                        'heading', ha='center', weight='bold')
        if panel in ['c', 'd']:
            supplement_text(s15_figure, 89, top + 4, s15_meta['subtitle'], 'annotation', ha='center')
    for row, top, label in [
        (reference_layout['s15']['rows'][1], s15_row_tops[1], '≤3 suspected SNPs\n(editing detection enabled)'),
        (reference_layout['s15']['rows'][2], s15_row_tops[2], '4+ suspected SNPs\n(editing detection disabled)')]:
        highest_tick = np.floor(row['ylim'][1] / row['dtick']) * row['dtick']
        center_top = top + S15_ROW_HEIGHT_MM * (1 - highest_tick / 2 / row['ylim'][1])
        supplement_text(s15_figure, 172, center_top, label, 'heading',
                        ha='center', va='center', rotation=270)

# %%
if 'FigureS15' in FIGURES_TO_BUILD:
    # Validate all facets, then use the common native-size export/crop workflow.
    assert len(s15_artist_checks) == 15
    assert sum(item['bars'] for item in s15_artist_checks) == 132
    assert sum(item['missing'] for item in s15_artist_checks) == 48
    s15_figure.canvas.draw()
    renderer = s15_figure.canvas.get_renderer()
    side_labels = [text for text in s15_figure.texts if text.get_rotation() == 270]
    side_bounds = [text.get_window_extent(renderer).extents * 25.4 / s15_figure.dpi
                   for text in side_labels]
    # The unchanged two-line labels fit beside their rows, within the 10 mm gap.
    assert all(box[3] - box[1] < S15_ROW_HEIGHT_MM + 4 for box in side_bounds)
    assert side_bounds[0][1] > side_bounds[1][3]
    s15_checks = {'scientific_artist_arrays_exact': True, 'sample_aggregation_exact': True,
        'facets': s15_facet_checks, 'artists': s15_artist_checks,
        'labels_moved_outside_for_6pt_readability': s15_label_adjustments,
        'input_sha256': s15_meta['checksums'], 'raw_positions_loaded': False,
        'row_height_mm': S15_ROW_HEIGHT_MM, 'row_tops_mm': s15_row_tops,
        'group_gap_mm': S15_GROUP_GAP_MM, 'internal_row_gap_mm': S15_INTERNAL_ROW_GAP_MM,
        'column_heading_gap_mm': S15_COLUMN_HEADING_GAP_MM,
        'side_label_bounds_mm': [box.tolist() for box in side_bounds]}
    # S15 precedes the batch setup; import its exporter without executing other figures.
    from manuscript_figures.word_layout import export_layout_figure
    from manuscript_figures import word_layout
    word_layout.ROOT = FINAL_FIGURES_DIR / 'render_work'
    export_layout_figure(s15_figure, 'FigureS15', s15_checks)


# %% [markdown]
# ## Remaining manuscript figures
#
# Select the required load, draw and export cells explicitly. Shared helpers come from `manuscript_figures`; scientific tables must be generated first. Finish the selected raw canvases with the per-page crop cell.
#

# %%
# New figures share the established local fonts and physical-unit conventions.
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
# Scientific tables share FIGURE_SOURCES; generated geometry preserves original contours.
from manuscript_figures.batch_helpers import load_table, axes_mm, text_mm, check_xy, export_figure, check_layout_spacing
from manuscript_figures.batch_helpers import close_facet_title, check_facet_title_attachment
assert hasattr(pd.DataFrame, 'map'), 'pandas 2.1 or newer is required'
for filename in ['OpenSans-Regular.ttf', 'OpenSans-Bold.ttf']:
    font_manager.fontManager.addfont(FONT_DIR / filename)
mpl.rcParams.update({'font.family': 'Open Sans', 'font.size': 6,
    'axes.labelsize': 7, 'axes.titlesize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6, 'legend.title_fontsize': 7, 'pdf.fonttype': 42,
    'svg.fonttype': 'none', 'svg.hashsalt': 'combinatorics-manuscript-batch',
    'savefig.facecolor': 'white', 'axes.labelcolor': '#2a3f5f', 'text.color': '#2a3f5f',
    'xtick.bottom': False, 'xtick.top': False, 'ytick.left': False, 'ytick.right': False,
    'xtick.minor.visible': False, 'ytick.minor.visible': False})
# Fixed native canvas dimensions (mm).
BATCH_LAYOUT = {'Figure4': {'canvas_mm': [180.0, 225.0]}, 'Figure5': {'canvas_mm': [180.0, 225.0]}, 'Figure6': {'canvas_mm': [180.0, 169]}, 'FigureS1': {'canvas_mm': [180.0, 190.0]}, 'FigureS3': {'canvas_mm': [180.0, 250.0]}, 'FigureS4': {'canvas_mm': [180.0, 129]}, 'FigureS5': {'canvas_mm': [180.0, 110]}, 'FigureS6': {'canvas_mm': [180.0, 260.0]}, 'FigureS8': {'canvas_mm': [180.0, 260]}, 'FigureS9': {'canvas_mm': [180.0, 260]}, 'FigureS10': {'canvas_mm': [180.0, 250]}, 'FigureS11': {'canvas_mm': [180.0, 255.0]}, 'FigureS12': {'canvas_mm': [180.0, 260.0]}, 'FigureS13': {'canvas_mm': [180.0, 150.0]}, 'FigureS14': {'canvas_mm': [180, 112]}, 'FigureS16': {'canvas_mm': [180.0, 235]}, 'FigureS17': {'canvas_mm': [180.0, 260]}, 'FigureS20': {'canvas_mm': [180.0, 205.0]}, 'FigureS6c_gene_maxima_review': {'canvas_mm': [180.0, 87.0]}, 'Figure2b': {'canvas_mm': [49, 45.5]}}
from manuscript_figures.word_layout import export_layout_figure as export_figure, place_clear_legends, balance_middle_panel

# Data roots and temporary/full-canvas export destination are explicit.
from manuscript_figures import batch_helpers, word_layout
batch_helpers.FINAL = FINAL_FIGURES_DIR
batch_helpers.FINALIZATION = FINAL_FIGURES_DIR / 'render_work'
word_layout.ROOT = FINAL_FIGURES_DIR / 'render_work'


# %% [markdown] figure="Figure4"
# ## Figure 4 — ranked expression and TPM recurrence
# Panel a preserves the original PacBio/Illumina ranked expression curves;
# panel b uses the separate PacBio3 cohort, three replicates, five TPM thresholds,
# and original Top-X versus All-isoform denominators. No points are sampled away.
# The 180 × 225 mm canvas provides 6 pt labels for the dense recurrence facets.
#
# This publication section loads the fixed selected-expression handoff documented in README.md; comparison variants remain separate.
#

# %% figure="Figure4"
if 'Figure4' in FIGURES_TO_BUILD:
    figure4a_table, figure4a_metadata = load_table('figure4a_expression', {'styles': [{'Platform': 'Long-reads', 'Gene': 'GRIA2', 'Color': '#00A08B'},
            {'Platform': 'Long-reads', 'Gene': 'PCLO', 'Color': '#AF0038'},
            {'Platform': 'Short-reads', 'Gene': 'RUSC2', 'Color': '#B68100'},
            {'Platform': 'Short-reads', 'Gene': 'TRIM2', 'Color': '#B2828D'},
            {'Platform': 'Short-reads', 'Gene': 'CA2D3', 'Color': '#6C4516'},
            {'Platform': 'Short-reads', 'Gene': 'ABL', 'Color': '#0D2A63'},
            {'Platform': 'Short-reads', 'Gene': 'DGLA', 'Color': '#620042'},
            {'Platform': 'Short-reads', 'Gene': 'K0513', 'Color': '#750D86'},
            {'Platform': 'Short-reads', 'Gene': 'KCNAS', 'Color': '#FB0D0D'},
            {'Platform': 'Short-reads', 'Gene': 'ACHA4', 'Color': '#DA16FF'},
            {'Platform': 'Short-reads', 'Gene': 'ANR17', 'Color': '#FC0080'},
            {'Platform': 'Short-reads', 'Gene': 'TWK7', 'Color': '#862A16'},
            {'Platform': 'Short-reads', 'Gene': 'SCN1', 'Color': '#6C7C32'},
            {'Platform': 'Short-reads', 'Gene': 'CACB2', 'Color': '#1616A7'},
            {'Platform': 'Short-reads', 'Gene': 'RIMS2', 'Color': '#FB00D1'},
            {'Platform': 'Short-reads', 'Gene': 'PCLO', 'Color': '#AF0038'},
            {'Platform': 'Short-reads', 'Gene': 'DOP1', 'Color': '#A777F1'},
            {'Platform': 'Short-reads', 'Gene': 'IQEC1', 'Color': '#EB663B'},
            {'Platform': 'Short-reads', 'Gene': 'CSKI1', 'Color': '#222A2A'},
            {'Platform': 'Short-reads', 'Gene': 'MTUS2', 'Color': '#1CA71C'},
            {'Platform': 'Short-reads', 'Gene': 'ROBO2', 'Color': '#2E91E5'}]})
    figure4b_table, figure4b_metadata = load_table('figure4b_recurrence', {'styles': {'10': {'color': '#DC3912', 'dash': 'longdashdot'},
            '50': {'color': '#FF9900', 'dash': 'solid'},
            '100': {'color': '#109618', 'dash': 'dot'},
            '250': {'color': '#990099', 'dash': 'dash'},
            '500': {'color': '#3366CC', 'dash': 'longdash'}}})
    assert len(figure4a_table) == 400407 and len(figure4b_table) == 220
    figure4 = plt.figure(figsize=np.array(BATCH_LAYOUT['Figure4']['canvas_mm']) / 25.4)
    figure4_points_checked = 0


# %% figure="Figure4"
if 'Figure4' in FIGURES_TO_BUILD:
    # Log rank and relative expression; original cumulative values stay linear.
    figure4a_axes = {}
    for col, platform in enumerate(['Long-reads', 'Short-reads']):
        for row, value in enumerate(['%RelativeExpression', '%CummulativeRelativeExpression']):
            axis = axes_mm(figure4, 20 + col * 63, [12, 57][row], 55, [38, 34][row], grid='both', hide_y=col == 1)
            axis.set_xscale('log')
            if row == 0:
                axis.set_yscale('log')
                axis.set_ylim(4e-5, 100)
                axis.set_yticks([1e-4, 1e-3, 1e-2, 0.1, 1, 10], ['0.0001', '0.001', '0.01', '0.1', '1', '10'])
                axis.yaxis.set_minor_formatter(NullFormatter())
                axis.set_title(platform, pad=4)
                axis.tick_params(labelbottom=False)
            else:
                axis.set_ylim(-2, 106)
                axis.set_yticks([0, 20, 40, 60, 80, 100])
            subset = figure4a_table.loc[figure4a_table['Platform'].eq(platform)]
            max_rank = subset['#Protein'].max()
            axis.set_xlim(0.6, max_rank * 1.35)
            ranks = [v for v in [1, 10, 100, 1000, 10000, 100000] if v <= max_rank * 1.35]
            axis.set_xticks(ranks, [f'{v // 1000}k' if v >= 1000 else str(v) for v in ranks])
            axis.xaxis.set_minor_formatter(NullFormatter())
            for style in figure4a_metadata['styles']:
                if style['Platform'] != platform:
                    continue
                curve = subset.loc[subset['Gene'].eq(style['Gene'])]
                dots = axis.scatter(curve['#Protein'], curve[value], s=0.55**2, color=style['Color'], alpha=0.3, linewidths=0)
                dots.set_gid(f'Figure4a-{row}-{platform}-{style["Gene"]}')
                check_xy(dots, curve['#Protein'], curve[value], 'scatter')
                figure4_points_checked += len(curve)
            if col == 0:
                axis.set_ylabel(['Relative expression [%]', 'Cumulative relative\nexpression [%]'][row], labelpad=4)
            figure4a_axes[(row, col)] = axis
    text_mm(figure4, 3, 8, 'a', 'panel')
    text_mm(figure4, 79, 101, 'Isoform rank', 'axis', ha='center')


# %% figure="Figure4"
if 'Figure4' in FIGURES_TO_BUILD:
    legend_top = 12
    for platform in ['Long-reads', 'Short-reads']:
        text_mm(figure4, 145, legend_top, platform, 'legend_heading')
        legend_top += 5.2
        for style in figure4a_metadata['styles']:
            if style['Platform'] != platform:
                continue
            figure4.add_artist(Line2D([146.5 / 180], [1 - legend_top / 225], transform=figure4.transFigure,
                marker='o', markersize=2.4, linestyle='none', color=style['Color']))
            text_mm(figure4, 150, legend_top, style['Gene'], 'legend', va='center')
            legend_top += 3
        legend_top += 5


# %% figure="Figure4"
if 'Figure4' in FIGURES_TO_BUILD:
    # Paired Top-X/All axes per gene; threshold colors and dash patterns are original.
    dash_styles = {'solid': '-', 'dot': (0, (1, 1.5)), 'dash': (0, (3.5, 2)),
                   'longdash': (0, (6, 3)), 'longdashdot': (0, (6, 2, 1, 2))}
    text_mm(figure4, 3, 113, 'b', 'panel')
    text_mm(figure4, 20, 112, 'Expression per replicate', 'legend_heading')
    for index, threshold in enumerate([10, 50, 100, 250, 500]):
        style = figure4b_metadata['styles'][str(threshold)]
        x = 20 + index * 29
        figure4.add_artist(Line2D([x / 180, (x + 5) / 180], [1 - 118 / 225] * 2,
            transform=figure4.transFigure, color=style['color'], linestyle=dash_styles[style['dash']], linewidth=0.7))
        text_mm(figure4, x + 6.5, 118, f'TPM ≥ {threshold}', 'legend', va='center')
    for row, definition in enumerate(['At least 2 replicates', 'All 3 replicates']):
        top = [132, 175][row]
        for col, gene in enumerate(['GRIA2', 'ADAR1', 'IQEC1']):
            left = 20 + col * 49
            top_values = [5, 10, 50, 100, 200, 500] + ([1000] if gene == 'IQEC1' else [])
            for all_isoforms in [False, True]:
                axis = axes_mm(figure4, left + (38 if all_isoforms else 0), top, 5 if all_isoforms else 36, 27,
                               hide_y=col > 0 or all_isoforms)
                axis.set_ylim(0, 102)
                axis.set_yticks([0, 20, 40, 60, 80, 100])
                axis.set_xlim((-0.32, 0.32) if all_isoforms else (-0.25, len(top_values) - 0.75))
                axis.set_xticks([0] if all_isoforms else range(len(top_values)), ['All'] if all_isoforms else [f'Top {v}' for v in top_values],
                                 rotation=0 if all_isoforms else 55, ha='center' if all_isoforms else 'right')
                for threshold in [10, 50, 100, 250, 500]:
                    data = figure4b_table.loc[figure4b_table['Gene'].eq(gene) & figure4b_table['RecurrenceDefinition'].eq(definition)
                        & figure4b_table['TPMThreshold'].eq(threshold) & figure4b_table['TopXIsoforms'].isna().eq(all_isoforms)]
                    style = figure4b_metadata['styles'][str(threshold)]
                    line, = axis.plot(data['TopXPosition'], data['%NumOfProteins'], color=style['color'],
                        linestyle='none' if all_isoforms else dash_styles[style['dash']], linewidth=0.65,
                        marker='o', markersize=1.7, markeredgewidth=0, alpha=0.85 if all_isoforms else 0.8, clip_on=False)
                    line.set_gid(f'Figure4b-{row}-{gene}-{all_isoforms}-{threshold}')
                    check_xy(line, data['TopXPosition'], data['%NumOfProteins'])
                    figure4_points_checked += len(data)
            if row == 0:
                text_mm(figure4, left + 21.5, 128, gene, ha='center')
        text_mm(figure4, 170, top + 13.5, definition, 'subtitle', rotation=-90, ha='center', va='center')
    text_mm(figure4, 8, 167, '% of isoforms in set', 'axis', rotation=90, ha='center', va='center')
    text_mm(figure4, 90.5, 219, 'Isoforms', 'axis', ha='center')


# %% figure="Figure4"
if 'Figure4' in FIGURES_TO_BUILD:
    assert figure4_points_checked == 801034
    export_figure(figure4, 'Figure4', {'all_xy_values_exact': True, 'points': figure4_points_checked,
        'expression_rows': 400407, 'recurrence_rows': 220, 'thresholds': [10, 50, 100, 250, 500],
        'original_denominators_and_Top1000_visibility_preserved': True})


# %% [markdown] figure="Figure5"
# ## Figure 5 — PCLO site dependencies and entropy
# The original masked correlation/MI caches carry 144 PacBio and 45 Illumina sites. Numeric site order is verified against the source helpers, pairwise caches and original read headers. Four undefined PacBio MI values remain NA. Heatmap cells are native vectors.
#
# Panel c loads `figure_sources/figure5c_entropy.tsv`, produced by the canonical Shannon subsections in the PacBio and Illumina notebooks. Observed entropy uses selected-isoform assigned-expression probabilities; hypothetical entropy sums categorical AA entropies from informative original All reads, with a denominator per codon. The two empirical constructions do not define an exact total-correlation gap. The table carries platform/gene identities and row order; publication colors and crop are declared below.

# %% figure="Figure5"
if 'Figure5' in FIGURES_TO_BUILD:
    from manuscript_figures.batch_helpers import continuous_colorbar
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle, Patch
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.cm import ScalarMappable
    figure5_table, figure5_metadata = load_table('figure5_pclo_matrices.BY', {'color_scales': {'Pearson': [[0, '#67001f'],
                              [0.1, '#b2182b'],
                              [0.2, '#d6604d'],
                              [0.3, '#f4a582'],
                              [0.4, '#fddbc7'],
                              [0.5, '#f7f7f7'],
                              [0.6, '#d1e5f0'],
                              [0.7, '#92c5de'],
                              [0.8, '#4393c3'],
                              [0.9, '#2166ac'],
                              [1, '#053061']],
                  'MI': [[0, '#f7f7f7'],
                         [0.25, '#d1e5f0'],
                         [0.5, '#67a9cf'],
                         [0.75, '#2166ac'],
                         [1, '#053061']]}})
    # Revised BY masks; MI and all original graphical metadata are preserved.
    figure5_metadata['site_order'] = {p: g.loc[g.ColumnIndex.eq(0), 'RowSite'].drop_duplicates().tolist()
        for p, g in figure5_table.groupby('Platform', sort=False)}
    entropy_path = figure_source('figure5c_entropy')
    entropy_summary_df = pd.read_table(entropy_path, float_precision='round_trip')
    assert entropy_summary_df.columns.tolist() == ['Platform', 'Gene', 'ObservedEntropyBits', 'HypotheticalEntropyBits']
    assert entropy_summary_df['Platform'].tolist() == ['PacBio'] * 2 + ['Illumina'] * 19
    assert not entropy_summary_df.duplicated(['Platform', 'Gene']).any()
    assert np.isfinite(entropy_summary_df[['ObservedEntropyBits', 'HypotheticalEntropyBits']]).all().all()
    entropy_gene_order = {
        label: entropy_summary_df.loc[entropy_summary_df.Platform.eq(platform), 'Gene'].tolist()
        for platform, label in [('PacBio', 'Long-reads'), ('Illumina', 'Short-reads')]
    }
    # Fixed publication colors; gene identities and order come from the scientific table.
    entropy_gene_colors = {
        "GRIA2": "#00A08B",
        "PCLO": "#AF0038",
        "ABL": "#0D2A63",
        "ACHA4": "#DA16FF",
        "ANR17": "#FC0080",
        "CA2D3": "#6C4516",
        "CACB2": "#1616A7",
        "CSKI1": "#222A2A",
        "DGLA": "#620042",
        "DOP1": "#A777F1",
        "IQEC1": "#EB663B",
        "K0513": "#750D86",
        "KCNAS": "#FB0D0D",
        "MTUS2": "#1CA71C",
        "RIMS2": "#FB00D1",
        "ROBO2": "#2E91E5",
        "RUSC2": "#B68100",
        "SCN1": "#6C7C32",
        "TRIM2": "#B2828D",
        "TWK7": "#862A16"
    }
    assert set(entropy_summary_df['Gene']) == set(entropy_gene_colors)
    entropy_table = entropy_summary_df.rename(columns={
        'ObservedEntropyBits': 'WeightedExpressionData', 'HypotheticalEntropyBits': 'Hypothetical'}).melt(
        id_vars=['Platform', 'Gene'], var_name='EntropyName', value_name='EntropyValue')
    entropy_table['Platform'] = entropy_table['Platform'].replace({'PacBio': 'Long-reads', 'Illumina': 'Short-reads'})
    assert len(entropy_table) == 42
    figure5_output_dir = FINAL_FIGURES_DIR
    figure5 = plt.figure(figsize=np.array([180.0, 225.0]) / 25.4)
    figure5_cell_checks = []

# %% figure="Figure5"
if 'Figure5' in FIGURES_TO_BUILD:
    # Exact reference RGB interpolation; each finite cell is a vector rectangle.
    for row, measure in enumerate(['Pearson', 'MI']):
        top = [20, 97][row]
        scale = figure5_metadata['color_scales'][measure]
        stops = np.array([point[0] for point in scale])
        channels = np.array([[int(point[1][k:k + 2], 16) for k in [1, 3, 5]] for point in scale])
        for col, platform in enumerate(['PacBio', 'Illumina']):
            selected = figure5_table.loc[figure5_table['Platform'].eq(platform) & figure5_table['Measure'].eq(measure)]
            count = len(figure5_metadata['site_order'][platform])
            values = selected['Value'].to_numpy().reshape(count, count)
            mask = selected['Masked'].to_numpy().reshape(count, count)
            np.testing.assert_array_equal(np.isnan(values), mask)
            axis = axes_mm(figure5, 18 + col * 72, top, 54, 54, grid=None)
            axis.set(xlim=(0, count), ylim=(count, 0), xticks=[], yticks=[])
            rows, columns = np.where(~mask)
            normalized = (values[~mask] + 1) / 2 if measure == 'Pearson' else values[~mask]
            colors = np.round(np.column_stack([np.interp(normalized, stops, channels[:,k]) for k in range(3)])) / 255
            # Preserve the reference's one-pixel gaps relative to each heatmap extent.
            gap_x, gap_y = count / (432 if measure == 'Pearson' else 438), count / 420
            rectangles = [Rectangle((x + gap_x / 2, y + gap_y / 2), 1 - gap_x, 1 - gap_y) for y, x in zip(rows, columns)]
            tiles = PatchCollection(rectangles, facecolors=colors, edgecolors='none', antialiased=False)
            tiles.set_gid(f'Figure5-{measure}-{platform}')
            axis.add_collection(tiles)
            np.testing.assert_array_equal(tiles.get_facecolors()[:, :3], colors)
            assert len(tiles.get_paths()) == (~mask).sum()
            figure5_cell_checks.append({'measure': measure, 'platform': platform, 'sites': count,
                'finite_cells': int((~mask).sum()), 'masked_cells': int(mask.sum()), 'exact_reference_rgb': True})
            text_mm(figure5, 45 + col * 72, top - 6, ['Long-reads', 'Short-reads'][col], ha='center')
        text_mm(figure5, 3, top - 12, ['a', 'b'][row], 'panel')
        text_mm(figure5, 18, top - 12, ['Site-Site correlations (PCLO)', 'Mutual information (PCLO)'][row])
        color_axis = axes_mm(figure5, 159, top + 10, 3, 40, grid=None)
        continuous_colorbar(color_axis, scale, -1 if measure == 'Pearson' else 0, 1,
            np.arange(-1, 1.01, .25) if measure == 'Pearson' else np.arange(0, 1.01, .2))


# %% figure="Figure5"
if 'Figure5' in FIGURES_TO_BUILD:
    text_mm(figure5, 3, 160, 'c', 'panel')
    entropy_bars_checked = 0
    for platform, left, width in [('Long-reads', 18, 17), ('Short-reads', 45, 126)]:
        genes = entropy_gene_order[platform]
        axis = axes_mm(figure5, left, 167, width, 35, hide_y=platform == 'Short-reads')
        axis.set(xlim=(-0.6, len(genes) - 0.4), ylim=(0, 37))
        axis.set_xticks(range(len(genes)), genes, rotation=60, ha='right')
        axis.set_yticks(range(0, 36, 5))
        for index, gene in enumerate(genes):
            color = entropy_gene_colors[gene]
            for shift, category in [(-0.18, 'WeightedExpressionData'), (0.18, 'Hypothetical')]:
                record = entropy_table.loc[entropy_table['Platform'].eq(platform) & entropy_table['Gene'].eq(gene)
                    & entropy_table['EntropyName'].eq(category)]
                assert len(record) == 1
                value = record['EntropyValue'].iloc[0]
                bar = axis.bar(index + shift, value, width=.32, color=color, linewidth=0,
                    hatch='////' if category == 'WeightedExpressionData' else None, edgecolor='white')[0]
                bar.set_gid(f'Figure5-entropy-{platform}-{gene}-{category}')
                assert bar.get_height() == value
                entropy_bars_checked += 1
        if platform == 'Long-reads':
            axis.set_ylabel('Entropy', labelpad=4)
        text_mm(figure5, left + width / 2, 162, platform, ha='center')
    # A corner-to-corner diagonal stays centered independently of hatch phase.
    from matplotlib.legend_handler import HandlerPatch
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as DrawingPath
    def entropy_legend_square(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
        side = min(width, height)
        x, y = -xdescent + (width-side)/2, -ydescent + (height-side)/2
        vertices = [(x,y), (x+side,y), (x+side,y+side), (x,y+side), (x,y)]
        codes = [DrawingPath.MOVETO, *[DrawingPath.LINETO]*3, DrawingPath.CLOSEPOLY]
        if orig_handle.get_label() == 'Observed':
            vertices += [(x,y), (x+side,y+side)]
            codes += [DrawingPath.MOVETO, DrawingPath.LINETO]
        return PathPatch(DrawingPath(vertices, codes))
    legend_items = [Patch(facecolor='white', edgecolor='#888888', linewidth=.6, label='Observed'),
                    Patch(facecolor='#888888', edgecolor='none', label='Hypothetical')]
    entropy_legend = axis.legend(handles=legend_items, loc='center right',
        bbox_to_anchor=(1, 32.5), bbox_transform=axis.get_yaxis_transform(),
        handler_map={Patch: HandlerPatch(patch_func=entropy_legend_square)},
        frameon=False, ncol=2, handlelength=.8, handleheight=.8, handletextpad=.5, columnspacing=1.5)
    for handle, name in zip(entropy_legend.legend_handles, ['Observed', 'Hypothetical']):
        handle.set_gid(f'Figure5-legend-{name}')
    text_mm(figure5, 94.5, 214, 'Gene', 'axis', ha='center')
    figure5.canvas.draw()
    legend_box = entropy_legend.get_window_extent(figure5.canvas.get_renderer())
    target_y_px = axis.transData.transform((0, 32.5))[1]
    np.testing.assert_allclose((legend_box.y0 + legend_box.y1)/2, target_y_px, atol=1e-8)
    figure5_legend_check = {'center_y_data': 32.5, 'bounds_mm': (legend_box.extents*25.4/figure5.dpi).tolist()}


# %% figure="Figure5"
if 'Figure5' in FIGURES_TO_BUILD:
    assert sum(record['finite_cells'] for record in figure5_cell_checks) == 20263  # BY Pearson + unchanged MI
    assert entropy_bars_checked == 42
    figure5.canvas.draw()
    square_checks = []
    for handle in entropy_legend.legend_handles:
        vertices = handle.get_transform().transform(handle.get_path().vertices)
        square = vertices[:4]
        np.testing.assert_allclose(np.ptp(square[:,0]), np.ptp(square[:,1]), atol=1e-10)
        if len(vertices) == 7:
            np.testing.assert_allclose(vertices[-2:].mean(axis=0), square.mean(axis=0), atol=1e-10)
        square_checks.append({'square': True, 'centered_diagonal': len(vertices) == 7,
                              'side_pt': float(np.ptp(square[:,0]) * 72 / figure5.dpi)})
    layout_checks = check_layout_spacing(figure5)
    # Preserve the official physical crop and write the three normal Figure5 files.
    from matplotlib.transforms import Bbox
    from matplotlib.text import Text
    from io import StringIO
    import base64
    import re
    import subprocess
    # Approved physical crop, in points and millimetres.
    figure5_layout = {'SVGHeader': '<svg xmlns:xlink="http://www.w3.org/1999/xlink" width="172.718661580mm" height="212.977337741mm" viewBox="1.9200000525 10.0799998090 489.5962194359 603.7152766458" xmlns="http://www.w3.org/2000/svg" version="1.1">', 'ViewBoxPoints': [1.9200000525, 10.079999809, 489.5962194359, 603.7152766458], 'WidthMM': 172.71866158, 'HeightMM': 212.977337741}
    left_pt, top_pt, width_pt, height_pt = figure5_layout['ViewBoxPoints']
    figure5_crop = Bbox.from_bounds(left_pt / 72, figure5.get_size_inches()[1] - (top_pt + height_pt) / 72,
        figure5_layout['WidthMM'] / 25.4, figure5_layout['HeightMM'] / 25.4)
    crop_bounds_mm = np.array(figure5_crop.extents) * 25.4
    for text_artist in figure5.findobj(Text):
        if text_artist.get_visible() and text_artist.get_text():
            bounds_mm = text_artist.get_window_extent(figure5.canvas.get_renderer()).extents * 25.4 / figure5.dpi
            assert (bounds_mm[:2] >= crop_bounds_mm[:2] - .02).all() and (bounds_mm[2:] <= crop_bounds_mm[2:] + .02).all(), text_artist.get_text()
    for figure5_axis in figure5.axes:
        for entropy_bar in figure5_axis.patches:
            if (entropy_bar.get_gid() or '').startswith('Figure5-entropy-'):
                assert 0 <= entropy_bar.get_height() <= figure5_axis.get_ylim()[1]
    for folder in ['svg', 'pdf', 'previews']:
        (figure5_output_dir / folder).mkdir(parents=True, exist_ok=True)
    figure5_svg_buffer = StringIO()
    figure5.savefig(figure5_svg_buffer, format='svg', metadata={'Creator': 'Combinatorics manuscript figures'})
    figure5_svg = figure5_svg_buffer.getvalue()
    figure5_font_css = ''
    for font_name, font_weight in [('OpenSans-Regular.ttf', 400), ('OpenSans-Bold.ttf', 700)]:
        font_data = base64.b64encode((FONT_DIR / font_name).read_bytes()).decode('ascii')
        figure5_font_css += "@font-face{font-family:'Open Sans';font-weight:" + str(font_weight) + ";src:url(data:font/ttf;base64," + font_data + ") format('truetype');}"
    figure5_svg = re.sub(r'<svg\b[^>]*>', lambda match: figure5_layout['SVGHeader'] + '<defs><style>' + figure5_font_css + '</style></defs>', figure5_svg, count=1)
    (figure5_output_dir / 'svg/Figure5.svg').write_text(figure5_svg)
    figure5.savefig(figure5_output_dir / 'pdf/Figure5.pdf', bbox_inches=figure5_crop, pad_inches=0,
        metadata={'Creator': 'Combinatorics manuscript figures'})
    subprocess.run(['pdftoppm', '-r', '600', '-singlefile', '-png',
        str(figure5_output_dir / 'pdf/Figure5.pdf'), str(figure5_output_dir / 'previews/Figure5')], check=True, timeout=120)

# %% [markdown] figure="Figure6"
# ## Figure 6 — octopus transcriptome and single-cell distributions
# Four eligible-gene cohorts retain their own denominators. The lower row
# uses the verified neuronal/non-neuronal FixedTMR50 caches. Original cumulative
# distributions, mean interpolation, reference bands and Mann–Whitney test are
# preserved. The full two-row assembly uses a 180 × 169 mm canvas.
#

# %% figure="Figure6"
if 'Figure6' in FIGURES_TO_BUILD:
    from matplotlib.patches import Rectangle
    figure6_table, figure6_metadata = load_table('figure6_distributions', {'styles': [{'cohort': 'WholeTranscriptome_TMR50', 'label': '>= 50 reads', 'color': 'purple', 'dash': 'dash'},
            {'cohort': 'WholeTranscriptome_TMR1000',
             'label': '>= 1000 reads',
             'color': 'green',
             'dash': 'dash'},
            {'cohort': 'SingleCell_Neuronal_FixedTMR50', 'label': 'Neural', 'color': 'red', 'dash': 'dash'},
            {'cohort': 'SingleCell_NonNeuronal_FixedTMR50',
             'label': 'Non-neural',
             'color': 'rgb(0,170,255)',
             'dash': 'dash'}]})
    figure6_cohorts, figure6_cohort_metadata = load_table('figure6_cohorts', {})
    figure6_metadata['denominators'] = figure6_cohorts.groupby('Cohort').size().to_dict()
    figure6_markers, _ = load_table('figure6_mean_markers', {})
    assert figure6_markers.Cohort.tolist() == ['WholeTranscriptome_TMR50', 'WholeTranscriptome_TMR1000']
    figure6_metadata['mean_markers'] = [
        {'cohort': r.Cohort, 'x': r.X, 'y': r.Y, 'color': color,
         'label': f'{r.X:.0f} distinct<br>proteins<br>(avg)' if i == 0 else f'{r.X:.0f} distinct proteins<br>(avg)'}
        for i, (r, color) in enumerate(zip(figure6_markers.itertuples(), ['purple', 'green']))]
    assert figure6_table.groupby('Cohort').DenominatorGenes.first().to_dict() == figure6_metadata['denominators']
    assert not figure6_cohorts.duplicated(['Cohort', 'Chrom']).any()
    figure6 = plt.figure(figsize=np.array([180, 174]) / 25.4)
    figure6_axes = []
    figure6_lower_y = figure6_table['CummulativeTranscripts'].min() ** 1.1
    figure6_points_checked = 0
    # Original SVG uses 9 px dashes on an 800 px canvas; retain that physical length.
    figure6_dash_pt=9*(180/800)*(72/25.4)


# %% figure="Figure6"
if 'Figure6' in FIGURES_TO_BUILD:
    for row, top in enumerate([18, 94]):
        axis = axes_mm(figure6, 20, top, 118, 63, grid='both')
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set(xlim=(0.8, 1200), ylim=(figure6_lower_y, 10**2.2))
        ticks = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
        axis.set_xticks(ticks, [str(v) for v in ticks])
        axis.set_yticks([.01, .1, 1, 10, 100], ['0.01', '0.1', '1', '10', '100'])
        axis.xaxis.set_minor_formatter(NullFormatter())
        minor_ticks = LogLocator(base=10, subs=[2, 5]).tick_values(*axis.get_ylim())
        axis.yaxis.set_minor_locator(mpl.ticker.FixedLocator([v for v in minor_ticks if axis.get_ylim()[0] <= v <= axis.get_ylim()[1]]))
        axis.yaxis.set_minor_formatter(FuncFormatter(lambda v,p: '2' if abs(np.log10(v/2)-round(np.log10(v/2)))<1e-8 else '5'))
        axis.grid(True, which='minor', axis='y', linewidth=.25, color='#dbe4ef')
        axis.tick_params(which='both', bottom=False, top=False, left=False, right=False)
        if row == 0:
            axis.tick_params(labelbottom=False)
            for left, width, color in [(1, 4, 'orange'), (5, 45, 'lightskyblue')]:
                axis.add_patch(Rectangle((left, figure6_lower_y), width, 100-figure6_lower_y,
                                        facecolor=color, alpha=.2, edgecolor='none', zorder=0))
        for style in figure6_metadata['styles'][row * 2:row * 2 + 2]:
            data = figure6_table.loc[figure6_table['Cohort'].eq(style['cohort'])]
            color = '#00aaff' if style['color'] == 'rgb(0,170,255)' else style['color']
            line, = axis.plot(data['NumOfProteins'], data['CummulativeTranscripts'],
                color=color, linestyle=(0, (figure6_dash_pt/.65, figure6_dash_pt/.65)) if row == 0 else (0, (figure6_dash_pt/.35, figure6_dash_pt/.35)), linewidth=.65 if row == 0 else .35,
                marker='o', markersize=1.7, markeredgewidth=0, clip_on=False)
            line.set_gid('Figure6-' + style['cohort'])
            check_xy(line, data['NumOfProteins'], data['CummulativeTranscripts'])
            figure6_points_checked += len(data)
        figure6_axes.append(axis)
        text_mm(figure6, 79, top - 2, ['Whole-transcriptome long-reads', 'Single-cell'][row], ha='center')
    text_mm(figure6, 20, 5, 'Octopus data')
    figure6_y_title = text_mm(figure6, 5, 87.5, '% of genes', 'axis', rotation=90, va='center', ha='center')
    figure6_x_title = text_mm(figure6, 79, 166, 'Distinct protein isoforms per gene', 'axis', ha='center')
    # Match visible title-to-tick-label gaps, rather than unequal anchor distances.
    figure6.canvas.draw()
    renderer = figure6.canvas.get_renderer()
    y_boxes = [t.get_window_extent(renderer) for a in figure6_axes
               for t in a.get_yticklabels(which='both') if t.get_visible() and t.get_text()]
    x_boxes = [t.get_window_extent(renderer) for t in figure6_axes[1].get_xticklabels() if t.get_visible()]
    y_clearance = min(b.x0 for b in y_boxes) - figure6_y_title.get_window_extent(renderer).x1
    x_top = min(b.y0 for b in x_boxes) - y_clearance
    shift = x_top - figure6_x_title.get_window_extent(renderer).y1
    figure6_x_title.set_y(figure6_x_title.get_position()[1] + shift / figure6.bbox.height)
    figure6_title_gap_mm = float(y_clearance * 25.4 / figure6.dpi)
    assert figure6_title_gap_mm > 0



# %% figure="Figure6"
if 'Figure6' in FIGURES_TO_BUILD:
    for index, marker in enumerate(figure6_metadata['mean_markers']):
        figure6_axes[0].plot(marker['x'], marker['y'], 'o', color=marker['color'], markersize=4, markeredgewidth=0)
        figure6_axes[0].annotate(marker['label'].replace('<br>', '\n'), (marker['x'], marker['y']),
            xytext=(-3, -4) if index == 0 else (3, 4), textcoords='offset points',
            fontsize=6, color=marker['color'], ha='right' if index == 0 else 'left',
            va='top' if index == 0 else 'bottom')
    for x, text in [(2.25, '~5 isoforms\nper gene\ndue to\nalternative\nsplicing'),
                    (17, '~50 distinct\npolypeptides\nper gene')]:
        figure6_axes[0].text(x, .3, text, fontsize=6, color='#555555', ha='center', va='center')
    figure6_statistics_text = (
        'Mann–Whitney U: P < 10$^{-100}$\n'
        'Wilcoxon signed-rank (2,981 shared genes): P < 10$^{-100}$'
    )
    figure6_axes[1].text(.97, .90, figure6_statistics_text, transform=figure6_axes[1].transAxes,
                        fontsize=6, ha='right', va='top', multialignment='left', linespacing=1.5,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85, 'pad': 1})
    for row, top in enumerate([18, 94]):
        text_mm(figure6, 146, top, ['Coverage', 'Cells'][row], 'legend_heading')
        for item, style in enumerate(figure6_metadata['styles'][row * 2:row * 2 + 2]):
            y = top + 5.2 + item * 4
            color = '#00aaff' if style['color'] == 'rgb(0,170,255)' else style['color']
            figure6.add_artist(Line2D([146/180, 151/180], [1-y/174]*2, transform=figure6.transFigure,
                color=color, linestyle=(0, (figure6_dash_pt/.65, figure6_dash_pt/.65)), linewidth=.65, marker='o', markersize=1.5, markeredgewidth=0))
            text_mm(figure6, 153, y, style['label'].replace('>=', '≥'), 'legend', va='center')


# %% figure="Figure6"
if 'Figure6' in FIGURES_TO_BUILD:
    assert figure6_points_checked == 192
    layout_checks = check_layout_spacing(figure6)
    export_figure(figure6, 'Figure6', {'equal_axis_title_clearance_mm':figure6_title_gap_mm, 'layout_clearance':layout_checks, 'four_curve_arrays_exact': True, 'points': figure6_points_checked,
        'denominators': figure6_metadata['denominators'], 'mean_markers_exact': 2,
        'fixed_TMR50_neuronal_cohorts': True, 'approved_P_value_annotation': figure6_statistics_text})


# %% [markdown] figure="Figure2_PCLO_top100_unique_reads_asset"
# ## Figure 2b — PCLO top-100 unique READ patterns
#
# Current 49 × 45.5 mm vector asset for BioRender, without a panel letter. Arial:
# axis titles 7 pt; rank numbers and legend 6 pt. Abundance uses NumOfReads; all
# 100 identities, cutoff ties, 144 site columns and NA states were verified against
# the full original table. Read the newly generated rank table under
# `MANUSCRIPT_FIGURE_SOURCES/Figure2`; write the standalone panel to the selected
# `MANUSCRIPT_FINAL_FIGURES` root. Supply the matching permitted font through
# `MANUSCRIPT_ARIAL`; assembly assets remain separate external inputs.
#

# %% figure="Figure2_PCLO_top100_unique_reads_asset"
if {'Figure2', 'Figure2b'} & set(FIGURES_TO_BUILD):
    # Load the approved matrix and the task-local Arial font.
    import shutil
    import subprocess
    from matplotlib.font_manager import FontProperties
    FIGURE2_DIR = Path(os.environ['MANUSCRIPT_FIGURE2_ASSETS']).resolve()
    FIGURE2_EVIDENCE = FINAL_FIGURES_DIR / 'render_work/figure2b/evidence'
    FIGURE2_EVIDENCE.mkdir(parents=True, exist_ok=True)
    figure2_input = figure_source('PCLO_top100_unique_reads')
    figure2_table = pd.read_table(figure2_input, dtype={'UniqueRead': str}, keep_default_na=False)
    assert not figure2_table.duplicated(['Gene', 'UniqueRead']).any()
    np.testing.assert_array_equal(figure2_table.UniqueReadRank, np.arange(1, 101))
    assert figure2_table.NumOfReads.is_monotonic_decreasing
    figure2_sites = [column for column in figure2_table if column.isdigit()]
    assert len(figure2_sites) == 144
    figure2_states = figure2_table[figure2_sites].to_numpy()
    assert figure2_states.shape == (100, 144) and (figure2_states == -1).sum() == 2255
    assert set(np.unique(figure2_states)) == {-1, 0, 1}
    figure2_ranking = figure2_table[['UniqueReadRank', 'UniqueRead', 'NumOfReads']].rename(
        columns={'UniqueReadRank': 'Rank', 'UniqueRead': 'PatternID'})
    figure2_ranking.to_csv(FIGURE2_EVIDENCE / 'Figure2b_top100_ranking.tsv', sep='\t', index=False)
    figure2_font = Path(os.environ['MANUSCRIPT_ARIAL']).resolve()
    assert figure2_font.is_file(), 'Supply the permitted matching Arial font'
    assert hashlib.sha256(figure2_font.read_bytes()).hexdigest() == '35c0f3559d8db569e36c31095b8a60d441643d95f59139de40e23fada819b833'


# %% figure="Figure2_PCLO_top100_unique_reads_asset"
if {'Figure2', 'Figure2b'} & set(FIGURES_TO_BUILD):
    # Local names and rc_context keep this Arial exception out of other figures.
    def draw_figure2b(states, font_path):
        from matplotlib.path import Path as VectorPath
        from matplotlib.patches import PathPatch
        import xml.etree.ElementTree as ET
        width_mm, height_mm = 49.0, 45.5
        np.testing.assert_array_equal(BATCH_LAYOUT['Figure2b']['canvas_mm'], [width_mm, height_mm])
        axis_font = FontProperties(family='Arial', fname=font_path, size=7)
        small_font = FontProperties(family='Arial', fname=font_path, size=6)
        assert axis_font.get_name() == 'Arial'
        with mpl.rc_context({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42,
                'svg.hashsalt': 'PCLO-Figure2b-BioRender-49x45p5',
                'savefig.facecolor': 'white', 'text.color': '#2a3f5f', 'axes.grid': False}):
            figure = plt.figure(figsize=(width_mm/25.4, height_mm/25.4))
            left, top, width, height = 7.6, 0.8, 27.1, 39.7
            axis = figure.add_axes([left/width_mm, 1-(top+height)/height_mm,
                                   width/width_mm, height/height_mm])
            colors = ['#c0c0c0', '#000000', '#4af626']
            # Compound vector paths retain every cell and avoid artificial white seams.
            axis.set_facecolor(colors[1])
            for state, color in zip([-1, 0, 1], colors):
                rows, columns = np.where(states == state)
                polygons = np.stack([np.column_stack([columns, rows]),
                    np.column_stack([columns+1, rows]), np.column_stack([columns+1, rows+1]),
                    np.column_stack([columns, rows+1])], axis=1)
                patch = PathPatch(VectorPath.make_compound_path_from_polys(polygons),
                    facecolor=color, edgecolor='none', antialiased=False, rasterized=False)
                patch.set_gid(f'Figure2b_state_{state}')
                axis.add_patch(patch)
            axis.set(xlim=(0, 144), ylim=(100, 0), xticks=[])
            axis.set_yticks(np.array([20, 40, 60, 80, 100])-.5, ['20', '40', '60', '80', '100'])
            axis.tick_params(axis='both', which='both', length=0, pad=1.5, colors='#2a3f5f')
            for spine in axis.spines.values(): spine.set_visible(False)
            for label in axis.get_yticklabels(): label.set_fontproperties(small_font)
            figure.text(1.6/width_mm, 1-(top+height/2)/height_mm,
                'Most expressed PCLO reads', rotation=90, ha='center', va='center', fontproperties=axis_font)
            figure.text((left+width/2)/width_mm, 1-43.4/height_mm,
                'Editing sites', ha='center', va='center', fontproperties=axis_font)
            legend_texts, legend_squares = [], []
            for center_y, label, color in [(17.0, 'Edited', colors[2]),
                    (20.65, 'Not-edited', colors[1]), (24.3, 'NA', colors[0])]:
                square = Rectangle((35.9/width_mm, 1-(center_y+1.1)/height_mm),
                    2.2/width_mm, 2.2/height_mm, transform=figure.transFigure,
                    facecolor=color, edgecolor='none')
                figure.add_artist(square)
                legend_squares.append(square)
                legend_texts.append(figure.text(38.9/width_mm, 1-center_y/height_mm, label,
                    ha='left', va='center', fontproperties=small_font))

            # Validate the complete physical text boxes before fixed-canvas export.
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            texts = [text for text in figure.findobj(Text) if text.get_visible() and text.get_text()]
            records = []
            for text in texts:
                box = text.get_window_extent(renderer)
                assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height, text.get_text()
                records.append({'text': text.get_text(), 'fontsize_pt': text.get_fontsize(),
                    'font': text.get_fontproperties().get_name(),
                    'bounds_mm': (np.array(box.extents)/figure.dpi*25.4).tolist()})
            for i, first in enumerate(texts):
                for second in texts[i+1:]:
                    assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)), (first.get_text(), second.get_text())
            for text, square in zip(legend_texts, legend_squares):
                assert text.get_window_extent(renderer).x0 > square.get_window_extent(renderer).x1
                assert text.get_window_extent(renderer).x0 > axis.get_window_extent(renderer).x1
            figure.savefig(FINAL_FIGURES_DIR / 'svg/Figure2b.svg', format='svg')
            figure.savefig(FINAL_FIGURES_DIR / 'pdf/Figure2b.pdf', format='pdf',
                metadata={'Title': 'PCLO top-100 unique-read patterns; Figure 2b'})
            figure.savefig(FINAL_FIGURES_DIR / 'previews/Figure2b.png', dpi=600)
            svg = ET.parse(FINAL_FIGURES_DIR / 'svg/Figure2b.svg')
            svg.getroot().set('width', '49mm'); svg.getroot().set('height', '45.5mm')
            ET.register_namespace('', 'http://www.w3.org/2000/svg')
            svg.write(FINAL_FIGURES_DIR / 'svg/Figure2b.svg', encoding='utf-8', xml_declaration=True)
            (FIGURE2_EVIDENCE / 'artists.json').write_text(json.dumps(records, indent=2))
        return figure

    figure2_rc_before = mpl.rcParams.copy()
    figure2_asset = draw_figure2b(figure2_states, figure2_font)
    assert dict(mpl.rcParams) == dict(figure2_rc_before)
    display(Image(filename=str(FINAL_FIGURES_DIR / 'previews/Figure2b.png')))


# %% [markdown] figure="FigureS1"
# ## Figure S1 — read coverage, quality and noise
# Exact read-level measurements and saved noise observations are persisted.
# Histogram bins are verified against every original bar. Original violin
# Bézier contours preserve the reference KDE; outliers use exact saved values.
# The full a–c assembly uses 180 × 190 mm, with space for 19 rotated gene labels.
#

# %% figure="FigureS1"
if 'FigureS1' in FIGURES_TO_BUILD:
    import hashlib
    from matplotlib.patches import PathPatch
    from matplotlib.transforms import Affine2D
    from manuscript_figures.violin_paths import original_violin_path
    s1_curves, s1_curve_meta = load_table('figure_s1_curves', {})
    s1_deletions, s1_deletion_meta = load_table('figure_s1_deletions', {})
    s1_noise, s1_noise_meta = load_table('figure_s1_noise', {'styles': [{'trace': 0, 'platform': 'PacBio', 'gene': 'GRIA2', 'color': '#00A08B'},
            {'trace': 1, 'platform': 'PacBio', 'gene': 'PCLO', 'color': '#AF0038'},
            {'trace': 2, 'platform': 'Illumina', 'gene': 'PCLO', 'color': '#AF0038'},
            {'trace': 3, 'platform': 'Illumina', 'gene': 'ABL', 'color': '#0D2A63'},
            {'trace': 4, 'platform': 'Illumina', 'gene': 'ACHA4', 'color': '#DA16FF'},
            {'trace': 5, 'platform': 'Illumina', 'gene': 'ANR17', 'color': '#FC0080'},
            {'trace': 6, 'platform': 'Illumina', 'gene': 'CA2D3', 'color': '#6C4516'},
            {'trace': 7, 'platform': 'Illumina', 'gene': 'CACB2', 'color': '#1616A7'},
            {'trace': 8, 'platform': 'Illumina', 'gene': 'CSKI1', 'color': '#222A2A'},
            {'trace': 9, 'platform': 'Illumina', 'gene': 'DGLA', 'color': '#620042'},
            {'trace': 10, 'platform': 'Illumina', 'gene': 'DOP1', 'color': '#A777F1'},
            {'trace': 11, 'platform': 'Illumina', 'gene': 'IQEC1', 'color': '#EB663B'},
            {'trace': 12, 'platform': 'Illumina', 'gene': 'K0513', 'color': '#750D86'},
            {'trace': 13, 'platform': 'Illumina', 'gene': 'KCNAS', 'color': '#FB0D0D'},
            {'trace': 14, 'platform': 'Illumina', 'gene': 'MTUS2', 'color': '#1CA71C'},
            {'trace': 15, 'platform': 'Illumina', 'gene': 'RIMS2', 'color': '#FB00D1'},
            {'trace': 16, 'platform': 'Illumina', 'gene': 'ROBO2', 'color': '#2E91E5'},
            {'trace': 17, 'platform': 'Illumina', 'gene': 'RUSC2', 'color': '#B68100'},
            {'trace': 18, 'platform': 'Illumina', 'gene': 'SCN1', 'color': '#6C7C32'},
            {'trace': 19, 'platform': 'Illumina', 'gene': 'TRIM2', 'color': '#B2828D'},
            {'trace': 20, 'platform': 'Illumina', 'gene': 'TWK7', 'color': '#862A16'}]})
    s1_octopus, s1_octopus_meta = load_table('figure_s1_octopus_histogram', {})
    s1_geometry_path = figure_source('figure_s1_violin_geometry', '.json')
    s1_geometry = json.loads(s1_geometry_path.read_text())
    s1_colors = {'GRIA2': '#00A08B', 'PCLO': '#AF0038'}
    figure_s1 = plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS1']['canvas_mm']) / 25.4)
    s1_points_checked = 0


# %% figure="FigureS1"
if 'FigureS1' in FIGURES_TO_BUILD:
    s1_a = axes_mm(figure_s1, 17, 24, 56, 49, grid='both')
    s1_b = axes_mm(figure_s1, 101, 24, 52, 49, grid=None)
    s1_b_percent = axes_mm(figure_s1, 101, 24, 52, 49, grid=None)
    s1_b_percent.patch.set_visible(False)
    s1_b_percent.yaxis.tick_right()
    s1_b_percent.yaxis.set_label_position('right')
    s1_b_percent.tick_params(which='both', left=False, right=False, labelleft=False, labelright=True,
                            bottom=False, top=False, labelbottom=False)
    for gene, color in s1_colors.items():
        for panel, axis, opacity, drawstyle in [('a', s1_a, .75, 'steps-post'), ('b', s1_b_percent, .5, 'default')]:
            data = s1_curves.loc[s1_curves.Panel.eq(panel) & s1_curves.Gene.eq(gene)]
            line, = axis.plot(data.X, data.Percent, color=color, alpha=opacity, linewidth=.65, drawstyle=drawstyle)
            check_xy(line, data.X, data.Percent)
            s1_points_checked += len(data)
        data = s1_deletions.loc[s1_deletions.Gene.eq(gene)]
        bars = s1_b.bar((data.Left+data.Right)/2, data.MeanDeletions, width=(data.Right-data.Left)*.9,
                       color=color, alpha=.5, linewidth=0)
        np.testing.assert_array_equal([bar.get_height() for bar in bars], data.MeanDeletions)
    s1_a.set(xlim=(0,6100), ylim=(0,105), xlabel='Read length', ylabel='Reads [%]')
    s1_a.set_xticks(np.arange(0,6001,1000))
    s1_a.set_yticks([0,25,50,75,100])
    s1_b.set(xlim=(.989975,1.000025), ylim=(0,s1_deletions.MeanDeletions.max()/.95),
             xlabel='Read quality', ylabel='Deletion events (avg)')
    s1_b_percent.set(xlim=s1_b.get_xlim(), ylim=(0,100), ylabel='Reads [%]')
    s1_b_percent.set_yticks([0,25,50,75,100])
    s1_b.set_xticks([.99,.992,.994,.996,.998,1], ['0.99','0.992','0.994','0.996','0.998','1'])
    s1_b.grid(axis='y', color='black', alpha=.4, linewidth=.3)
    s1_b_percent.grid(axis='y', color='lightpink', linestyle='--', linewidth=.4)
    s1_heading_gap_mm, s1_letter_gap_mm = 3.0, 4.0
    s1_heading_artists, s1_panel_letters = [], []
    for x, title, letter in [(17,'Long-reads coverage','a'),(101,'Long-reads quality','b')]:
        s1_heading_artists.append(text_mm(figure_s1,x,24-s1_heading_gap_mm,title,va='bottom'))
        s1_panel_letters.append(text_mm(figure_s1,x-s1_letter_gap_mm,24-s1_heading_gap_mm,
                                       letter,'panel',ha='right',va='bottom'))
    # Vertical right-side legends; each handle retains its trace's alpha and meaning.
    from matplotlib.patches import Patch
    for axis, left, filled in [(s1_a, 76, False), (s1_b, 168, True)]:
        handles = [Patch(facecolor=color, edgecolor='none', alpha=.5, label=gene) if filled else
                   Line2D([], [], color=color, alpha=.75, linewidth=.65, label=gene)
                   for gene, color in s1_colors.items()]
        axis.legend(handles=handles, title='Gene', loc='upper left',
            bbox_to_anchor=(left/180, 1-25/190), bbox_transform=figure_s1.transFigure,
            frameon=False, borderaxespad=0, handlelength=1, handletextpad=.5,
            labelspacing=.6, borderpad=0, alignment='left')


# %% figure="FigureS1"
if 'FigureS1' in FIGURES_TO_BUILD:
    s1_panel_letters.append(text_mm(figure_s1,17-s1_letter_gap_mm,92,'c','panel',ha='right'))
    text_mm(figure_s1,17,92,'Noise levels')
    for platform, x, width, title in [('PacBio',17,16,'Squid’s\nLong-reads'),
                                     ('Illumina',45,88,'Squid’s\nShort-reads')]:
        axis = axes_mm(figure_s1,x,110,width,46)
        styles = [s for s in s1_noise_meta['styles'] if s['platform'] == platform]
        for position, style in enumerate(styles):
            geometry = s1_geometry[style['trace']]
            path = original_violin_path(geometry['d'])
            transform = Affine2D().scale(len(styles)/geometry['plot_width_px'], -9/270).translate(-.5,9)
            axis.add_patch(PathPatch(path, transform=transform+axis.transData,
                facecolor=mpl.colors.to_rgba(style['color'],.5), edgecolor=style['color'], linewidth=.35))
            values = np.array(geometry['outlier_values'])
            data = s1_noise.loc[s1_noise.Trace.eq(style['trace']), 'PercentNoise']
            assert np.isin(values, data).all()
            points = axis.scatter(np.full(len(values),position),values,s=.65,color=style['color'],linewidths=0)
            check_xy(points,np.full(len(values),position),values,'scatter')
        axis.set(xlim=(-.5,len(styles)-.5), ylim=(0,9))
        axis.set_xticks(range(len(styles)), [s['gene'] for s in styles], rotation=60, ha='right')
        axis.set_yticks([0,2,4,6,8])
        text_mm(figure_s1,x+width/2,102,title,ha='center',va='center')
        text_mm(figure_s1,x+width/2,169.5,'Gene','axis',ha='center')
        if platform == 'PacBio':
            axis.set_ylabel('Per-site noise level [%]')
    axis = axes_mm(figure_s1,145,110,30,46,grid='x')
    axis.set_xscale('log')
    data = s1_octopus.loc[s1_octopus.Genes.gt(0)]
    bars = axis.barh((data.Left+data.Right)/2,data.Genes,height=data.Right-data.Left,color='black',linewidth=0)
    np.testing.assert_array_equal([b.get_width() for b in bars],data.Genes)
    axis.set(xlim=(.7,1100), ylim=(0,9), ylabel='Per-gene noise level [%]')
    axis.set_xticks([1,10,100,1000], ['1','10','100','1000'])
    axis.tick_params(axis='x',labelrotation=60)
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.tick_params(which='both', bottom=False, top=False, left=False, right=False)
    axis.set_yticks([0,2,4,6,8])
    text_mm(figure_s1,160,102,'Whole-transcriptome\noctopus data',ha='center',va='center')
    text_mm(figure_s1,160,169.5,'Genes','axis',ha='center')


# %% figure="FigureS1"
if 'FigureS1' in FIGURES_TO_BUILD:
    assert s1_points_checked == 348248 and s1_octopus.Genes.sum() == 1711
    c_bounds = np.array([a.get_position().bounds for a in figure_s1.axes[-3:]]) * [180,190,180,190]
    c_gutters = c_bounds[1:,0] - (c_bounds[:-1,0]+c_bounds[:-1,2])
    np.testing.assert_allclose(c_gutters, [12,12], atol=1e-10)
    # Measure visible text edges so the same gaps survive every fresh export.
    figure_s1.canvas.draw()
    renderer = figure_s1.canvas.get_renderer()
    mm_per_pixel = 25.4 / figure_s1.dpi
    title_gaps = [(t.get_window_extent(renderer).y0-a.get_window_extent(renderer).y1)*mm_per_pixel
                  for t,a in zip(s1_heading_artists,[s1_a,s1_b])]
    letter_gaps = [(a.get_window_extent(renderer).x0-t.get_window_extent(renderer).x1)*mm_per_pixel
                   for t,a in zip(s1_panel_letters,[s1_a,s1_b,figure_s1.axes[3]])]
    np.testing.assert_allclose(title_gaps,s1_heading_gap_mm,atol=1e-9)
    np.testing.assert_allclose(letter_gaps,s1_letter_gap_mm,atol=1e-9)
    b_letter = s1_panel_letters[1].get_window_extent(renderer)
    previous_panel_right = max(s1_a.get_window_extent(renderer).x1,s1_a.get_legend().get_window_extent(renderer).x1)
    assert (b_letter.x0-previous_panel_right)*mm_per_pixel > s1_letter_gap_mm
    layout_checks = check_layout_spacing(figure_s1)
    export_figure(figure_s1,'FigureS1',{'heading_gaps_mm':title_gaps, 'panel_letter_gaps_mm':letter_gaps, 'b_letter_closer_to_own_panel':True, 'layout_clearance':layout_checks, 'panel_c_gutters_mm':c_gutters.tolist(), 'exact_curve_points':s1_points_checked,'histogram_means':402,
        'noise_observations':len(s1_noise),'original_violin_contours':21,'exact_outliers':1118,
        'octopus_denominator':1711,'no_new_KDE_estimation':True})


# %% [markdown] figure="FigureS2"
# ## Figure S2 — original cohorts and corrected PCLO selection
#
# Panels a–f retain their 180 × 260 mm native canvas before the existing exterior crop.
# Panels a–d and f use the existing verified inputs. Panel e uses the fixed PCLO
# catalog/call-based selection in `various_plots.BQ30.ipynb`, retaining missing
# measurements and showing both current platforms where measurements are available.
#

# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    from matplotlib_venn import venn2
    from matplotlib.patches import Rectangle, PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D
    s2_venns,_=load_table('s2_venns', {})
    s2_pooled,s2_pooled_meta=load_table('s2_pooled', {})
    s2_distribution,_=load_table('s2_distribution', {})
    s2_boxes,_=load_table('s2_boxes', {})
    s2_pclo,_=load_table('s2_pclo', {})
    s2_motifs,s2_motif_meta=load_table('s2_motifs', {})
    s2_mm12=pd.read_table(figure_source('combined_s2_s15_mm12_panel_counts'),float_precision='round_trip')
    s2_mm12=s2_mm12.loc[s2_mm12.Figure.eq(5)]
    s2_mm12_meta={'dataset_order': ['SquidLongReads', 'SquidShortReads', 'Octopus'],
 'mismatch_order': ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G'],
 'mismatch_colors': {'A>C': '#2E91E5',
                     'A>G': '#E15F99',
                     'A>T': '#1CA71C',
                     'C>A': '#FB0D0D',
                     'C>G': '#DA16FF',
                     'C>T': '#222A2A',
                     'G>A': '#B68100',
                     'G>C': '#750D86',
                     'G>T': '#EB663B',
                     'T>A': '#511CFB',
                     'T>C': '#00A08B',
                     'T>G': '#FB00D1'}}
    figure_s2=plt.figure(figsize=np.array([180,260])/25.4)
    s2_counts={'a':0,'b':0,'c':0,'d':0,'e':0,'f':0}
    s2_extra_text=[]

    s2_statistics, _ = load_table('s2_statistics', {})
    assert len(s2_statistics) == 1
    s2_pooled_meta = s2_statistics.iloc[0].to_dict()


# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    text_mm(figure_s2,2,6,'a','panel')
    text_mm(figure_s2,15,6,'Mismatches above editing threshold, but below 100% mismatch frequency')
    text_mm(figure_s2,15,10,'Only genes with at most 3 suspected SNPs are included (editing detection enabled)','annotation')
    labels=['Squid’s Long-reads','Squid’s Short-reads','Whole-transcriptome octopus data']
    for column,(dataset,label) in enumerate(zip(s2_mm12_meta['dataset_order'],labels)):
        axis=axes_mm(figure_s2,16+54*column,17,46,29,hide_y=column!=0)
        data=s2_mm12.loc[s2_mm12.Dataset.eq(dataset)].set_index('Mismatch').loc[s2_mm12_meta['mismatch_order']]
        values=data.Percentage.to_numpy();axis.set(ylim=(0,103),xlim=(-.6,11.6))
        bars=axis.bar(range(12),values,color=[s2_mm12_meta['mismatch_colors'][m] for m in data.index],width=.7)
        axis.set_yticks(range(0,101,20));axis.set_xticks(range(12),data.index,rotation=-45,ha='left',rotation_mode='anchor')
        for bar,value in zip(bars,values):
            if value>0:
                inside=value>10
                axis.text(bar.get_x()+bar.get_width()/2,value-2 if inside else value+1,
                    f'{value:.1f}' if value<99.95 else '100',fontsize=6,
                    ha='center',va='top' if inside else 'bottom',rotation=90 if inside else 0)
        close_facet_title(axis,label)
        if column==0:axis.set_ylabel('Percentage',labelpad=2)
        if column==1:axis.set_xlabel('Mismatch',labelpad=0)
        s2_counts['a']+=len(data)

# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    text_mm(figure_s2,2,64,'b','panel')
    text_mm(figure_s2,22,64,'Squid’s Long-reads',ha='center')
    text_mm(figure_s2,115,64,'Squid’s Short-reads',ha='center')
    s2_octopus_title=text_mm(figure_s2,21,96,'Whole-transcriptome\n octopus data',ha='center')
    s2_venn_title_checks=[]
    s2_group_bottoms={}
    s2_group_extents={}
    for platform in ['PacBio','Illumina','Octopus']:
        data=s2_venns.loc[s2_venns.Platform.eq(platform)]
        for facet,record in enumerate(data.itertuples()):
            if platform=='PacBio':left,top,width,height=2+23*facet,73,20,12
            elif platform=='Illumina':left,top,width,height=54+17.5*(facet%7),73+27*(facet//7),17,12
            else:left,top,width,height=4,123,35,16
            s2_group_bottoms.setdefault(platform,[]).append(top+height+1.5)
            s2_group_extents.setdefault(platform,[]).append((left,left+width))
            axis=axes_mm(figure_s2,left,top,width,height,grid=None)
            diagram=venn2((record.NewOnly,record.KnownOnly,record.Both),set_labels=None,ax=axis)
            for artist in diagram.subset_labels:
                if artist is not None:artist.set(fontsize=6,color='black',fontfamily='Open Sans')
            axis.set_axis_off()
            if platform!='Octopus':
                title=text_mm(figure_s2,left+width/2,top-1.2,record.Gene,ha='center',va='bottom')
                s2_extra_text.append(title)
            else:
                title=s2_octopus_title;title.set_verticalalignment('bottom')
            s2_venn_title_checks.append((title,diagram,platform,facet))
            for center,label,count in [(left+width*.24,'De-novo',record.NewOnly+record.Both),
                                       (left+width*.76,'Known',record.KnownOnly+record.Both)]:
                s2_extra_text.append(text_mm(figure_s2,center,top+height+1.5,f'{label}\n({count})','legend',ha='center',va='top'))
            s2_counts['b']+=1

# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    text_mm(figure_s2,2,155,'c','panel')
    axis=axes_mm(figure_s2,16,157,42,35,grid='both')
    for label,color in [('Both','black'),('De-novo','green'),('Known','red')]:
        data=s2_pooled.loc[s2_pooled.Detected.eq(label)]
        points=axis.scatter(data.PercentEditing.fillna(0),data.PercentEditingKnown.fillna(0),
            s=1.8,color=color,label=label,linewidths=0)
        check_xy(points,data.PercentEditing.fillna(0),data.PercentEditingKnown.fillna(0),'scatter')
        s2_counts['c']+=len(data)
    axis.set(xlim=(0,100),ylim=(0,100),xticks=range(0,101,20),yticks=range(0,101,20))
    axis.set_xlabel('Editing level (targeted sequencing) [%]',labelpad=2)
    axis.set_ylabel('Previously reported\nediting level [%]',labelpad=2)
    close_facet_title(axis,'Pooled editing levels in squid')
    axis.legend(title='Detected',loc='upper left',bbox_to_anchor=(1.01,1),frameon=False,
                handlelength=.6,handletextpad=.3,borderaxespad=0,markerscale=1.5)
    r=s2_pooled_meta['pearson_r'];p=s2_pooled_meta['pearson_p']
    p_text=f'p-value = {p:.2g}' if p>0 else 'p-value < 1E−22'
    axis.text(.03,.97,f"Pearson’s r = {r:.2g}\n{p_text}",transform=axis.transAxes,
        fontsize=6,va='top',bbox={'facecolor':'white','edgecolor':'none','pad':1})


# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    text_mm(figure_s2,2,206,'d','panel')
    axis=axes_mm(figure_s2,16,208,42,22)
    for record,color in zip(s2_boxes.itertuples(),['black','green','red']):
        data=s2_distribution.loc[s2_distribution.Detected.eq(record.Detected)]
        axis.scatter(data.DisplayX,data.DisplayValue,s=.35,color=color,linewidths=0,alpha=.65)
        center=record.Category
        axis.add_patch(Rectangle((center-.245,record.Q1),.49,record.Q3-record.Q1,
            facecolor=color,edgecolor=color,alpha=.5,linewidth=.4))
        axis.plot([center-.245,center+.245],[record.Median]*2,color=color,linewidth=.5)
        for boundary,whisker in [(record.Q1,record.Low),(record.Q3,record.High)]:
            axis.plot([center,center],[boundary,whisker],color=color,linewidth=.5)
            axis.plot([center-.05,center+.05],[whisker]*2,color=color,linewidth=.5)
        s2_counts['d']+=len(data)
    axis.set(xlim=(-.5,2.5),ylim=(0,100),xticks=range(3),xticklabels=['Both','De-novo','Known'],yticks=[0,25,50,75,100])
    axis.set_xlabel('Detected',labelpad=1);axis.set_ylabel('Editing level [%]',labelpad=2)
    close_facet_title(axis,'Pooled editing levels in squid')


# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    text_mm(figure_s2,76,155,'e','panel')
    colors=['#269da5','#cf951c','#60a03b']
    columns=['PacBio','Illumina','EditingFrequencyKnown']
    labels=['Long-reads','Short-reads','Known (Short-reads)']
    for known,top,height in [(True,157,21),(False,183,21)]:
        axis=axes_mm(figure_s2,90,top,82,height,grid='both')
        data=s2_pclo.loc[s2_pclo.Edited_Known.eq(known)]
        if not known:data=data.loc[data.Edited_PacBio | data.Edited_Illumina]
        for column,color,label in zip(columns,colors,labels):
            if not known and column=='EditingFrequencyKnown':continue
            points=data.loc[data[column].notna()]  # Show measurements even when this platform did not call the site.
            artist=axis.scatter(points.Position,points[column],s=2,color=color,alpha=.95,linewidths=0,label=label)
            check_xy(artist,points.Position,points[column],'scatter');s2_counts['e']+=len(points)
        axis.set_xlim(s2_pclo.Position.min()-80,s2_pclo.Position.max()+80)
        axis.set_xticks([2500,3500,4500,5500])
        if known:
            axis.set(ylim=(0,100),yticks=[0,20,40,60,80,100]);axis.set_xticklabels([])
            handles=[Line2D([],[],marker='o',linestyle='none',color=color,
                markersize=3,markeredgewidth=0,alpha=1,label=label)
                for color,label in zip(colors,labels)]
            axis.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,1.12),ncol=3,frameon=False,
                columnspacing=.8,handlelength=.5,handletextpad=.3,borderaxespad=0)
        else:
            axis.set(yscale='linear',ylim=(0,5),yticks=[0,1,2,3,4,5])
            axis.yaxis.set_minor_locator(mpl.ticker.NullLocator())
            axis.set_xlabel('Position',labelpad=1)
        axis.tick_params(which='both',left=False,right=False,bottom=False,top=False)
        text_mm(figure_s2,175,top+height/2,'Known\npositions' if known else 'De-novo\npositions',
                'annotation',rotation=-90,ha='center',va='center')
    text_mm(figure_s2,82,180.5,'Editing level in squid’s PCLO [%]','axis',rotation=90,ha='center',va='center')


# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    text_mm(figure_s2,76,218,'f','panel')
    logo_colors={'A':'#008000','C':'#0000ff','G':'#ffa600','T':'#ff0000'}
    titles=['Squid\nLong-reads','Squid\nShort-reads','Whole-transcriptome\noctopus data']
    for col,(platform,title) in enumerate(zip(['PacBio','Illumina','Octopus'],titles)):
        axis=axes_mm(figure_s2,88+31*col,224,25,12,grid=None)
        data=s2_motifs.loc[s2_motifs.Platform.eq(platform)]
        for offset,group in data.groupby('Offset',sort=True):
            bottom=0
            for row in group.sort_values(['Probability','Base']).itertuples():
                probability=row.Probability
                if probability==0:continue
                glyph=TextPath((0,0),row.Base,size=1,prop=font_manager.FontProperties(family='Open Sans',weight='bold'))
                bounds=glyph.get_extents()
                transform=(Affine2D().translate(-bounds.x0,-bounds.y0)
                    .scale(.95/bounds.width,.98*probability/bounds.height)
                    .translate(offset-.475,bottom+.01*probability)+axis.transData)
                axis.add_patch(PathPatch(glyph,transform=transform,facecolor=logo_colors[row.Base],linewidth=0))
                bottom+=probability;s2_counts['f']+=1
            np.testing.assert_allclose(bottom,1,atol=1e-14)
        axis.set(xlim=(-1.5,1.5),ylim=(0,1),xticks=[-1,0,1],yticks=[0,.5,1])
        close_facet_title(axis,title)


# %%
if 'FigureS2' in FIGURES_TO_BUILD:
    np.testing.assert_allclose(max(s2_group_bottoms['Illumina']),max(s2_group_bottoms['Octopus']),atol=1e-10)
    assert min(x[0] for x in s2_group_extents['Illumina'])-max(x[1] for x in s2_group_extents['PacBio']) == 9
    assert not s2_pclo.duplicated(['Chrom','Position']).any()
    s2_expected_e=int(s2_pclo[['PacBio','Illumina']].notna().sum().sum()
        + s2_pclo.loc[s2_pclo.Edited_Known,'EditingFrequencyKnown'].notna().sum())
    assert s2_counts=={'a':36,'b':22,'c':1331,'d':1331,'e':s2_expected_e,'f':27},s2_counts
    figure_s2.canvas.draw()
    for title,diagram,platform,facet in s2_venn_title_checks:
        top=max(patch.get_window_extent().y1 for patch in diagram.patches if patch is not None)
        title.set_y((top+.6*figure_s2.dpi/25.4)/figure_s2.bbox.height)
    figure_s2.canvas.draw()
    s2_renderer=figure_s2.canvas.get_renderer()
    for title,diagram,platform,facet in s2_venn_title_checks:
        top=max(patch.get_window_extent().y1 for patch in diagram.patches if patch is not None)
        gap=(title.get_window_extent(s2_renderer).y0-top)*25.4/figure_s2.dpi
        np.testing.assert_allclose(gap,.6,atol=1e-8)
    boxes=[(artist,artist.get_window_extent(s2_renderer)) for artist in s2_extra_text]
    for index,(artist,box) in enumerate(boxes):
        for other,other_box in boxes[index+1:]:
            assert not box.overlaps(other_box),('S2 Venn labels overlap',artist.get_text(),other.get_text())
    s2_position_label=next(a.xaxis.label for a in figure_s2.axes if a.xaxis.label.get_text()=='Position')
    s2_ef_gap=(s2_position_label.get_window_extent(s2_renderer).y0
        - max(a.title.get_window_extent(s2_renderer).y1 for a in figure_s2.axes[-3:]))*25.4/figure_s2.dpi
    assert s2_ef_gap>=2,('S2 E/F clearance',s2_ef_gap)
    export_figure(figure_s2,'FigureS2',{'source':'Original A-D/F cohorts; E uses the fixed PCLO call-based cohort',
        'panels':s2_counts,'reference_validation':'finalization_2026-10-12/evidence/S2_recovery_validation.json',
        'venn_label_overlaps':0,'clipping_checked':True,'E_F_clearance_mm':s2_ef_gap,
        'pclo_sites':len(s2_pclo),'pclo_known_sites':int(s2_pclo.Edited_Known.sum()),
        'pclo_new_sites':int((~s2_pclo.Edited_Known).sum())})


# %% [markdown] figure="FigureS3"
# ## Figure S3 — non-synonymous substitutions per read
# The two PacBio genes appear above all 19 Illumina genes in original order.
# Integer-bin counts and original population mean/SD are preserved. Mean±SD
# labels are paired with Min/Max legend entries to avoid collisions over peaks.
# Canvas: 180 × 250 mm; the five supplementary facet rows retain 6 pt labels.
#

# %% figure="FigureS3"
if 'FigureS3' in FIGURES_TO_BUILD:
    s3_facets = []
    s3_counts,s3_count_meta=load_table('figure_s3_counts', {})
    s3_stats,s3_stat_meta=load_table('figure_s3_statistics', {})
    assert not s3_counts.duplicated(['Platform','Gene','Estimate','Substitutions']).any()
    figure_s3=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS3']['canvas_mm'])/25.4)
    s3_checked_bins=0


# %% figure="FigureS3"
if 'FigureS3' in FIGURES_TO_BUILD:
    for platform in ['PacBio','Illumina']:
        genes=s3_stats.loc[s3_stats.Platform.eq(platform),'Gene'].drop_duplicates().tolist()
        for index,gene in enumerate(genes):
            if platform=='PacBio':
                x,top,width,height=17+index*42,22,32,20
                first_column=index==0
            else:
                x,top,width,height=17+(index%4)*42,67+(index//4)*33,32,20
                first_column=index%4==0
            axis=axes_mm(figure_s3,x,top,width,height,hide_y=not first_column)
            styles=s3_stats.loc[s3_stats.Platform.eq(platform)&s3_stats.Gene.eq(gene)]
            for estimate_index,style in enumerate(styles.itertuples()):
                data=s3_counts.loc[s3_counts.Platform.eq(platform)&s3_counts.Gene.eq(gene)&s3_counts.Estimate.eq(style.Estimate)]
                bars=axis.bar(data.Substitutions,data.Reads,width=1,color=style.Color,alpha=.75,linewidth=0,clip_on=False)
                np.testing.assert_array_equal([b.get_height() for b in bars],data.Reads)
                assert data.Reads.sum()==style.DenominatorReads
                s3_checked_bins+=len(data)
                axis.plot([style.Mean]*2,[0,data.Reads.max()],color='white',linestyle=':',linewidth=.6)
                legend_y=(8000 if platform=='PacBio' else 300000)/style.YLimit-.11-estimate_index*.20
                axis.plot(.53,legend_y,'s',transform=axis.transAxes,color=style.Color,markersize=4,markeredgewidth=0)
                axis.text(.61,legend_y,f'{style.Estimate}: {style.Annotation}',transform=axis.transAxes,
                    fontsize=6,va='center',ha='left',bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':0})
            axis.set(xlim=(-.6,styles.XLimit.iloc[0]+.6),ylim=(0,styles.YLimit.iloc[0]))
            axis.set_xticks([0,20,40,60,80] if platform=='PacBio' else [0,10,20,30,40])
            axis.set_yticks([0,2000,4000,6000,8000] if platform=='PacBio' else [0,100000,200000,300000])
            axis.yaxis.set_major_formatter(FuncFormatter(lambda value,_: '0' if value==0 else f'{value/1000:g}k'))
            s3_facets.append(close_facet_title(axis, gene, grid_top=8000 if platform=='PacBio' else 300000))
        # Put each platform heading 1.5 mm above its rendered gene-title boxes.
        figure_s3.canvas.draw()
        renderer = figure_s3.canvas.get_renderer()
        group_facets = s3_facets[-len(genes):]
        gene_title_top = max(title.get_window_extent(renderer).y1 for _, title, _ in group_facets)
        heading_top = 250 - gene_title_top * 25.4 / figure_s3.dpi - 1.5
        text_mm(figure_s3,17,heading_top,'Squid’s Long-reads' if platform=='PacBio' else 'Squid’s Short-reads',va='bottom')
        text_mm(figure_s3,54 if platform=='PacBio' else 92,52 if platform=='PacBio' else 227,'Non-syn substitutions per read','axis',ha='center')
        text_mm(figure_s3,4,32 if platform=='PacBio' else 139,'Reads','axis',rotation=90,ha='center',va='center')


# %% figure="FigureS3"
if 'FigureS3' in FIGURES_TO_BUILD:
    assert s3_checked_bins==len(s3_counts)
    sizes_mm = np.array([a.get_position().bounds[2:] for a in figure_s3.axes]) * [180,250]
    np.testing.assert_allclose(sizes_mm, np.tile([32,20],(21,1)), atol=1e-10)
    title_attachment = check_facet_title_attachment(s3_facets)
    layout_checks = check_layout_spacing(figure_s3)
    export_figure(figure_s3,'FigureS3',{'layout_clearance':layout_checks, 'title_attachment':title_attachment, 'facets':21,'exact_count_bins':s3_checked_bins,
        'exact_mean_std_annotations':42,'denominators_verified':True,'original_cohort_gene_order':True,
        'boundary_bars_complete':True,'native_numerical_mean_std_order_preserved':True})


# %% [markdown] figure="FigureS4"
# ## Figure S4 — editing across PacBio3 samples
# The original 3 × 4 threshold grid retains all 144 percentages, gene order,
# replicate colors and the separate All-replicates intersection. The denominator
# definition is recorded in the input metadata; integer counts are not inferred
# from rounded ratios. Canvas: 180 × 129 mm.
#

# %% figure="FigureS4"
if 'FigureS4' in FIGURES_TO_BUILD:
    from matplotlib.patches import Rectangle
    s4_table,s4_metadata=load_table('figure_s4_plot', {'gene_order': ['GRIA2', 'ADAR1', 'IQEC1'], 'replicate_order': ['1', '2', '3', 'All']})
    s4_table['Color'] = s4_table.Color.map(lambda color: tuple(int(v)/255 for v in color.removeprefix('rgb(').removesuffix(')').split(',')))
    assert len(s4_table)==144 and s4_table.MinTotalCoverage.eq(0).all()
    figure_s4=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS4']['canvas_mm'])/25.4)
    s4_checked=0


# %% figure="FigureS4"
if 'FigureS4' in FIGURES_TO_BUILD:
    for row,top in enumerate([27,59,91]):
        for col,x in enumerate([17,52,87,122]):
            axis=axes_mm(figure_s4,x,top,30,27,hide_y=col>0)
            for rep_index,replicate in enumerate(s4_metadata['replicate_order']):
                data=s4_table.loc[s4_table.Row.eq(row)&s4_table.Column.eq(col)&s4_table.Replicate.eq(replicate)]
                assert data.Gene.tolist()==s4_metadata['gene_order']
                centers=np.arange(3)+(rep_index-1.5)*.2
                bars=axis.bar(centers,data.PercentEditedOfCoveredSites,width=.18,color=data.Color.iloc[0],linewidth=0)
                np.testing.assert_array_equal([b.get_height() for b in bars],data.PercentEditedOfCoveredSites)
                s4_checked+=len(data)
            axis.set(xlim=(-.5,2.5),ylim=(0,100/.95))
            axis.set_yticks(np.arange(0,101,10))
            axis.set_xticks(range(3),s4_metadata['gene_order'])
            axis.tick_params(labelbottom=row==2)
            if row==0:
                threshold=[.1,.3,.5,1][col]
                text_mm(figure_s4,x+15,26.4,f'Editing / sample ≥ {threshold:g}%',ha='center',va='bottom')
        text_mm(figure_s4,154.3,top+13.5,['All sites','Editing ≥ 1%','Editing ≥ 5%'][row],
                rotation=-90,ha='center',va='center')
    text_mm(figure_s4,4,72.5,'Sites showing editing [%]','axis',rotation=90,ha='center',va='center')
    text_mm(figure_s4,160,29,'Sample','legend_heading')
    for index,replicate in enumerate(s4_metadata['replicate_order']):
        y=35+index*4
        color=s4_table.loc[s4_table.Replicate.eq(replicate),'Color'].iloc[0]
        figure_s4.add_artist(Rectangle((160/180,1-(y+1)/129),2/180,2/129,
            transform=figure_s4.transFigure,facecolor=color,edgecolor='none'))
        text_mm(figure_s4,164,y,replicate,'legend',va='center')


# %% figure="FigureS4"
if 'FigureS4' in FIGURES_TO_BUILD:
    assert s4_checked==144
    layout_checks = check_layout_spacing(figure_s4)
    export_figure(figure_s4,'FigureS4',{'layout_clearance':layout_checks, 'exact_percentages':144,'original_threshold_grid':[3,4],
        'original_gene_replicate_order':True,'denominator_definition_preserved':True,
        'integer_counts_not_inferred':True})


# %% [markdown] figure="FigureS5"
# ## Figure S5 — editing by haplotype in six octopus genes
# The six original gene/site/haplotype cohorts retain 494 means, including six
# absent LIMD2 values. Only the 488 finite values become vector markers.
# Haplotype legends sit inside the upper-right of each wide, low panel, as in the reference.
# Canvas: 180 × 110 mm, original 2 × 3 facet order and haplotype colors.
#

# %% figure="FigureS5"
if 'FigureS5' in FIGURES_TO_BUILD:
    s5_table,s5_metadata=load_table('figure_s5_plot', {'gene_order': ['RBGP1', 'HEXB', 'TOLIP', 'LIMD2', 'RXR', 'SMUF2']})
    assert s5_table.MeanEditingPercent.isna().sum()==6
    assert len(s5_table[['Chrom','EditingSite']].drop_duplicates())==247
    figure_s5=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS5']['canvas_mm'])/25.4)
    s5_checked=0


# %% figure="FigureS5"
if 'FigureS5' in FIGURES_TO_BUILD:
    for facet,gene in enumerate(s5_metadata['gene_order']):
        x=18+(facet%3)*52
        top=13+(facet//3)*48
        axis=axes_mm(figure_s5,x,top,46,35,grid='both')
        cohort=s5_table.loc[s5_table.Facet.eq(facet)]
        xmin,xmax=cohort.EditingSite.min(),cohort.EditingSite.max()
        ymax=cohort.MeanEditingPercent.max()
        axis.set(xlim=(xmin-(xmax-xmin)*.06,xmax+(xmax-xmin)*.06),ylim=(-ymax*.07,ymax*1.07))
        for item,(trace,data) in enumerate(cohort.groupby('Trace',sort=False)):
            values=data.loc[data.MeanEditingPercent.notna()]
            points=axis.scatter(values.EditingSite,values.MeanEditingPercent,s=8,
                color=data.Color.iloc[0],alpha=.8,linewidths=0,clip_on=False)
            points.set_gid(f'FigureS5-trace-{trace}')
            check_xy(points,values.EditingSite,values.MeanEditingPercent,'scatter')
            s5_checked+=len(values)
        handles=[Line2D([],[],marker='o',linestyle='none',markersize=3.5,
            markeredgewidth=0,color=d.Color.iloc[0],label=d.Haplotype.iloc[0])
            for _,d in cohort.groupby('Trace',sort=False)]
        legend=axis.legend(handles=handles,title='Haplotype',loc='upper right',
            frameon=True,facecolor='white',edgecolor='#b0b0b0',framealpha=.9,
            borderpad=.35,handletextpad=.3,labelspacing=.25)
        legend.get_title().set_fontsize(7)
        legend.get_title().set_weight('bold')
        text_mm(figure_s5,x+23,top-3,gene,ha='center')
        # Use the observed maximum, before the unchanged axis headroom is added.
        y_tick_step = 10 if ymax <= 60 else 20
        axis.yaxis.set_major_locator(mpl.ticker.MultipleLocator(y_tick_step))
        axis.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4,integer=True,min_n_ticks=3))
        axis.set_xticks([v for v in axis.get_xticks() if axis.get_xlim()[0] <= v <= axis.get_xlim()[1]])
        axis.set_yticks([v for v in axis.get_yticks() if axis.get_ylim()[0] <= v <= axis.get_ylim()[1]])
        axis.yaxis.set_major_formatter(FuncFormatter(lambda value,_:f'{value:g}'))
    text_mm(figure_s5,93,105.5,'Editing site','axis',ha='center')
    text_mm(figure_s5,4,54.5,'Mean editing [%]','axis',rotation=90,ha='center',va='center')


# %% figure="FigureS5"
if 'FigureS5' in FIGURES_TO_BUILD:
    assert s5_checked==488
    s5_y_ticks = []
    for facet, (gene, axis) in enumerate(zip(s5_metadata['gene_order'], figure_s5.axes)):
        maximum = s5_table.loc[s5_table.Facet.eq(facet), 'MeanEditingPercent'].max()
        expected_step = 10 if maximum <= 60 else 20
        np.testing.assert_allclose(np.diff(axis.get_yticks()), expected_step)
        s5_y_ticks.append({'gene': gene, 'maximum_data_y': float(maximum), 'dtick': expected_step})
    # Candidate boxes include the rendered frame and the full marker radius.
    s5_legend_checks = place_clear_legends(figure_s5, s5_metadata['gene_order'])
    hexb = s5_table.loc[s5_table.Gene.eq('HEXB') & s5_table.EditingSite.eq(2365), 'MeanEditingPercent']
    np.testing.assert_allclose(sorted(hexb), [52.49221183800623,58.70206489675516], rtol=0, atol=1e-12)
    layout_checks = check_layout_spacing(figure_s5)
    export_figure(figure_s5,'FigureS5',{'legend_clearance':s5_legend_checks, 'layout_clearance':layout_checks, 'y_tick_rule': s5_y_ticks, 'points':488,'original_site_haplotype_rows':494,
        'missing_means_preserved':6,'identified_sites':247,'exact_plotted_arrays':True,
        'reference_in_panel_legends':True,'cohort_and_denominator_definition_preserved':True})


# %% [markdown] figure="FigureS6"
# ## Figure S6 — dispersion, PCLO recoding and codon/isoform relationship
# The original a–c order and cohorts are preserved. Panel b uses all 131
# unfiltered codons; it continues on a second axis to retain every label at
# 6 pt. Its original category order (including three final A codons) is retained.
# Panel c displays 21 maxima, while its original Pearson calculation uses all
# 336 solution repetitions. Canvas: 180 × 260 mm.
#

# %% figure="FigureS6"
if 'FigureS6' in FIGURES_TO_BUILD:
    from matplotlib.patches import Rectangle
    s6_dispersion,s6_dispersion_meta=load_table('figure_s6a_squid', {})
    s6_histogram,s6_histogram_meta=load_table('figure_s6a_histogram', {})
    s6_recoding,s6_recoding_meta=load_table('figure_s6b_recoding', {})
    s6_maxima,s6_maxima_meta=load_table('figure_s6c_maxima', {'original_x_limits': [0, 144.10000000000002], 'original_y_limits': [0, 20958.300000000003]})
    s6_solutions,s6_solutions_meta=load_table('figure_s6c_solutions', {})
    assert len(s6_solutions)==336 and len(s6_maxima)==21
    assert s6_recoding.Count.sum()==1725 and s6_histogram.Genes.sum()==1262
    figure_s6=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS6']['canvas_mm'])/25.4)

    s6_recoding_meta['codon_order'] = s6_recoding.sort_values('CodonOrder').Codon.drop_duplicates().tolist()

    s6_statistics, _ = load_table('figure_s6c_statistics', {})
    assert len(s6_statistics) == 1 and s6_statistics.N.iloc[0] == len(s6_maxima)
    s6_maxima_meta['annotation'] = s6_statistics.Annotation.iloc[0]


# %%
s6_histogram

# %% figure="FigureS6"
if 'FigureS6' in FIGURES_TO_BUILD:
    text_mm(figure_s6,2,8,'a','panel')
    text_mm(figure_s6,16,8,'Dispersion of estimates for the number of distinct protein isoforms')
    for platform,x,width,title in [('PacBio',16,16,'Squid’s\nLong-reads'),('Illumina',44,89,'Squid’s\nShort-reads')]:
        axis=axes_mm(figure_s6,x,26,width,28)
        data=s6_dispersion.loc[s6_dispersion.Platform.eq(platform)]
        bars=axis.bar(np.arange(len(data)),data.PercentDispersion,color=data.Color,width=.8,linewidth=0)
        np.testing.assert_array_equal([b.get_height() for b in bars],data.PercentDispersion)
        axis.set(xlim=(-.5,len(data)-.5),ylim=(0,2))
        axis.set_xticks(range(len(data)),data.Gene,rotation=60,ha='right')
        axis.set_yticks([0,.5,1,1.5,2],['0','0.5','1','1.5','2'])
        text_mm(figure_s6,x+width/2,19,title,ha='center',va='center')
        text_mm(figure_s6,x+width/2,66.3,'Gene','axis',ha='center')
    axis=axes_mm(figure_s6,145,26,30,28,grid='x')
    axis.set_xscale('log')
    data=s6_histogram.loc[s6_histogram.Genes.gt(0)]
    bars=axis.barh((data.Left+data.Right)/2,data.Genes,height=(data.Right-data.Left)*.8,color='black',linewidth=0)
    np.testing.assert_array_equal([b.get_width() for b in bars],data.Genes)
    axis.set(xlim=(.7,1600),ylim=(-.25,50.25))
    axis.set_xticks([1,10,100,1000],['1','10','100','1000'])
    axis.tick_params(axis='x',labelrotation=60)
    axis.set_yticks([0,10,20,30,40,50])
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.tick_params(which='both',left=False,right=False,bottom=False,top=False)
    text_mm(figure_s6,159,19,'Whole-transcriptome\noctopus data',ha='center',va='center')
    text_mm(figure_s6,160,66.3,'Genes','axis',ha='center')
    text_mm(figure_s6,4,40,'Dispersion [%]','axis',rotation=90,ha='center',va='center')


# %% figure="FigureS6"
if 'FigureS6' in FIGURES_TO_BUILD:
    s6_upper_artists = [*figure_s6.axes, *figure_s6.texts, *figure_s6.artists]
    text_mm(figure_s6,2,76,'b','panel')
    text_mm(figure_s6,18,76,'Recoding events in PCLO 100 most expressed isoforms')
    text_mm(figure_s6,18,83,'Recoding substitution','legend_heading',va='center')
    legend=s6_recoding.loc[s6_recoding.ShowLegend].drop_duplicates('Trace')
    for index,style in enumerate(legend.itertuples()):
        x=59+index*10.3
        figure_s6.add_artist(Rectangle((x/180,1-84/260),2/180,2/260,
            transform=figure_s6.transFigure,facecolor=style.Color,edgecolor='none'))
        text_mm(figure_s6,x+2.7,83,style.NewAA,'legend',va='center')
    s6_bars_checked=0
    for start,stop,top in [(0,66,88),(66,131,128)]:
        axis=axes_mm(figure_s6,18,top,157,20)
        stacked=np.zeros(stop-start,dtype=int)
        for trace,data in s6_recoding.groupby('Trace',sort=False):
            data=data.loc[data.CodonOrder.ge(start)&data.CodonOrder.lt(stop)]
            if data.empty:continue
            positions=data.CodonOrder.to_numpy()-start
            bars=axis.bar(positions,data.Count,bottom=stacked[positions],width=.8,
                          color=data.Color.iloc[0],linewidth=0)
            np.testing.assert_array_equal([b.get_height() for b in bars],data.Count)
            np.testing.assert_array_equal([b.get_y() for b in bars],stacked[positions])
            stacked[positions]+=data.Count.to_numpy()
            s6_bars_checked+=len(data)
        axis.set(xlim=(-.5,stop-start-.5),ylim=(0,90))
        axis.set_yticks([0,20,40,60,80])
        axis.set_xticks(range(stop-start),s6_recoding_meta['codon_order'][start:stop],rotation=90)
    text_mm(figure_s6,4,118,'Count','axis',rotation=90,ha='center',va='center')
    text_mm(figure_s6,96.5,167,'Codon position (genomically encoded AA)','axis',ha='center')
    assert s6_bars_checked==len(s6_recoding)


# %% figure="FigureS6"
if 'FigureS6' in FIGURES_TO_BUILD:
    s6_middle_artists = [a for a in [*figure_s6.axes, *figure_s6.texts, *figure_s6.artists]
                         if a not in s6_upper_artists]
    text_mm(figure_s6,2,178,'c','panel')
    axis=axes_mm(figure_s6,19,182,93,51,grid='both')
    for style in s6_maxima.itertuples():
        marker='*' if style.Marker=='star' else 'o'
        point=axis.scatter([style.EditableCodons],[style.MaximumIsoforms],marker=marker,
            s=36 if marker=='*' else 12,color=style.Color,linewidths=0,clip_on=False)
        check_xy(point,[style.EditableCodons],[style.MaximumIsoforms],'scatter')
    axis.set(xlim=s6_maxima_meta['original_x_limits'],ylim=s6_maxima_meta['original_y_limits'],
             xlabel='Editable codons',ylabel='Distinct protein isoforms')
    axis.set_xticks(np.arange(0,141,20))
    axis.set_yticks([0,5000,10000,15000,20000],['0','5k','10k','15k','20k'])
    axis.text(110,2500,s6_maxima_meta['annotation'],fontsize=6,ha='center',va='center',
              bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':1})
    for platform,heading_top in [('Long-reads',182),('Short-reads',199)]:
        text_mm(figure_s6,123,heading_top,platform,'legend_heading')
        data=s6_maxima.loc[s6_maxima.Platform.eq(platform)]
        for index,style in enumerate(data.itertuples()):
            column=index//10 if platform=='Short-reads' else 0
            row=index%10 if platform=='Short-reads' else index
            x=123+column*27
            y=heading_top+5+row*3.6
            figure_s6.add_artist(Line2D([x/180],[1-y/260],transform=figure_s6.transFigure,
                marker='*' if style.Marker=='star' else 'o',color=style.Color,
                markersize=5 if style.Marker=='star' else 3.5,markeredgewidth=0,linestyle='none'))
            text_mm(figure_s6,x+3,y,style.Gene,'legend',va='center')


# %% figure="FigureS6"
if 'FigureS6' in FIGURES_TO_BUILD:
    bounds_mm = np.array([a.get_position().bounds for a in figure_s6.axes[:3]]) * [180,260,180,260]
    np.testing.assert_allclose(bounds_mm[1:,0]-(bounds_mm[:-1,0]+bounds_mm[:-1,2]), [12,12], atol=1e-10)
    s6_lower_artists = [a for a in [*figure_s6.axes, *figure_s6.texts, *figure_s6.artists]
                        if a not in s6_upper_artists+s6_middle_artists]
    s6_panel_spacing = balance_middle_panel(figure_s6, s6_upper_artists, s6_middle_artists, s6_lower_artists)
    layout_checks = check_layout_spacing(figure_s6)
    export_figure(figure_s6,'FigureS6',{'panel_gap_balance':s6_panel_spacing, 'layout_clearance':layout_checks, 'dispersion_genes':1283,'histogram_denominator':1262,
        'unfiltered_codons':131,'isoforms_in_recoding_panel':100,'exact_recoding_counts':1725,
        'recoding_table_rows_checked':s6_bars_checked,'original_codon_category_order':True,
        'maxima_plotted':21,'pearson_solution_rows':336,'original_pearson_annotation_verified':True})


# %% [markdown] figure="FigureS8"
# ## Figure S8 — relative expression and original fitted curves
# Full 2-gene long-read section above the original 19-gene short-read grid.
# Canvas: 180 × 260 mm; all original data, colors and facet order retained.
# The selected reference is the assembled manuscript PDF, pages 14/15;
# loose SVG assets contain different historical curves. Shared S8/S9 data
# differ from the separately saved Figure 4 merged-cache plot.
# S9 Top-N markers now use the exact 1-based Rank=N value; the previous offset is retained in evidence.
#
# This publication section loads the fixed selected-expression handoff documented in README.md; comparison variants remain separate.
#

# %% figure="FigureS8"
if 'FigureS8' in FIGURES_TO_BUILD:
    s8_facets = []
    s8_expression,s8_meta=load_table('figure_s8_s9_expression', {'styles': [{'Platform': 'Long-reads', 'Gene': 'GRIA2', 'Facet': 0, 'Color': '#00A08B'},
            {'Platform': 'Long-reads', 'Gene': 'PCLO', 'Facet': 1, 'Color': '#AF0038'},
            {'Platform': 'Short-reads', 'Gene': 'RUSC2', 'Facet': 0, 'Color': '#B68100'},
            {'Platform': 'Short-reads', 'Gene': 'TRIM2', 'Facet': 1, 'Color': '#B2828D'},
            {'Platform': 'Short-reads', 'Gene': 'CA2D3', 'Facet': 2, 'Color': '#6C4516'},
            {'Platform': 'Short-reads', 'Gene': 'ABL', 'Facet': 3, 'Color': '#0D2A63'},
            {'Platform': 'Short-reads', 'Gene': 'DGLA', 'Facet': 4, 'Color': '#620042'},
            {'Platform': 'Short-reads', 'Gene': 'K0513', 'Facet': 5, 'Color': '#750D86'},
            {'Platform': 'Short-reads', 'Gene': 'KCNAS', 'Facet': 6, 'Color': '#FB0D0D'},
            {'Platform': 'Short-reads', 'Gene': 'ACHA4', 'Facet': 7, 'Color': '#DA16FF'},
            {'Platform': 'Short-reads', 'Gene': 'ANR17', 'Facet': 8, 'Color': '#FC0080'},
            {'Platform': 'Short-reads', 'Gene': 'TWK7', 'Facet': 9, 'Color': '#862A16'},
            {'Platform': 'Short-reads', 'Gene': 'SCN1', 'Facet': 10, 'Color': '#6C7C32'},
            {'Platform': 'Short-reads', 'Gene': 'CACB2', 'Facet': 11, 'Color': '#1616A7'},
            {'Platform': 'Short-reads', 'Gene': 'RIMS2', 'Facet': 12, 'Color': '#FB00D1'},
            {'Platform': 'Short-reads', 'Gene': 'PCLO', 'Facet': 13, 'Color': '#AF0038'},
            {'Platform': 'Short-reads', 'Gene': 'DOP1', 'Facet': 14, 'Color': '#A777F1'},
            {'Platform': 'Short-reads', 'Gene': 'IQEC1', 'Facet': 15, 'Color': '#EB663B'},
            {'Platform': 'Short-reads', 'Gene': 'CSKI1', 'Facet': 16, 'Color': '#222A2A'},
            {'Platform': 'Short-reads', 'Gene': 'MTUS2', 'Facet': 17, 'Color': '#1CA71C'},
            {'Platform': 'Short-reads', 'Gene': 'ROBO2', 'Facet': 18, 'Color': '#2E91E5'}]})
    s8_details,s8_detail_meta=load_table('figure_s8_fits', {})
    figure_s8=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS8']['canvas_mm'])/25.4)
    s8_points=0
    s8_detail_points=0
    # One shared range contains every original point and persisted fit, on both platforms.
    s8_all_y = np.r_[s8_expression.RelativeExpression, s8_details.FittedRelativeExpression]
    assert np.isfinite(s8_all_y).all() and (s8_all_y > 0).all()
    s8_y_limits = (10**-4.144520345527512, 10**1.5315880109715003)
    if s8_all_y.min() < s8_y_limits[0] or s8_all_y.max() > s8_y_limits[1]:
        s8_y_limits = (min(s8_y_limits[0], s8_all_y.min()/1.05),
                       max(s8_y_limits[1], s8_all_y.max()*1.05))
    s8_equation_checks = []

    s8_equations, _ = load_table('figure_s8_equations', {})
    assert not s8_equations.duplicated(['Platform', 'Gene']).any()
    s8_meta['styles'] = pd.DataFrame(s8_meta['styles']).merge(s8_equations[['Platform', 'Gene', 'Equation']], on=['Platform', 'Gene'], validate='one_to_one').to_dict('records')


# %% figure="FigureS8"
if 'FigureS8' in FIGURES_TO_BUILD:
    for style in s8_meta['styles']:
        platform,gene,facet=style['Platform'],style['Gene'],style['Facet']
        long_reads=platform=='Long-reads'
        column=facet if long_reads else facet%4
        row=0 if long_reads else facet//4
        left=18+column*39
        top=22 if long_reads else 68+row*36.2
        width,height=(34,24)
        axis=axes_mm(figure_s8,left,top,width,height,grid='both',hide_y=column!=0)
        data=s8_expression.loc[s8_expression.Platform.eq(platform)&s8_expression.Gene.eq(gene)]
        details=s8_details.loc[s8_details.Platform.eq(platform)&s8_details.Gene.eq(gene)]

        points=axis.scatter(data.Rank,data.RelativeExpression,s=1.5 if long_reads else .7,
            color=style['Color'],alpha=.2,linewidths=0)
        check_xy(points,data.Rank,data.RelativeExpression,'scatter')
        line,=axis.plot(details.Rank,details.FittedRelativeExpression,color='gray',linestyle='--',linewidth=1)
        check_xy(line,details.Rank,details.FittedRelativeExpression)
        axis.set(xscale='log',yscale='log',xlim=(.8,80000),
            ylim=s8_y_limits)
        axis.set_xticks([1,100,10000], ['1','100','10k'])
        axis.set_yticks([.0001,.001,.01,.1,1,10], ['0.0001','0.001','0.01','0.1','1','10'])
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.yaxis.set_minor_locator(mpl.ticker.NullLocator())
        axis.yaxis.set_minor_formatter(NullFormatter())
        equation=style['Equation'].replace(' * ',' × ')
        equation_artist=axis.text(.97,np.sqrt(10),equation,transform=axis.get_yaxis_transform(),
            fontsize=6,ha='right',va='center',color='gray')
        figure_s8.canvas.draw()
        box=equation_artist.get_window_extent(figure_s8.canvas.get_renderer())
        target=axis.transData.transform((1,np.sqrt(10)))[1]
        np.testing.assert_allclose((box.y0+box.y1)/2,target,atol=1e-8)
        margin_box=box.padded(.12*figure_s8.dpi/25.4)
        point_xy=axis.transData.transform(np.c_[data.Rank,data.RelativeExpression])
        assert not any(margin_box.contains(*point) for point in point_xy), (platform,gene,'equation overlaps data')
        fit_xy=axis.transData.transform(np.c_[details.Rank,details.FittedRelativeExpression])
        assert not mpl.path.Path(fit_xy).intersects_bbox(margin_box,filled=False), (platform,gene,'equation overlaps fit')
        assert box.x0>=axis.bbox.x0 and box.x1<=axis.bbox.x1
        s8_equation_checks.append({'platform':platform,'gene':gene,'center_y':float(np.sqrt(10)),
                                  'bounds_mm':(box.extents*25.4/figure_s8.dpi).tolist()})

        axis.tick_params(which='both',left=False,right=False,bottom=False,top=False)
        s8_facets.append(close_facet_title(axis, gene))
        s8_points+=len(data);s8_detail_points+=len(details)
    text_mm(figure_s8,18,14,'Squid’s Long-reads')
    text_mm(figure_s8,18,60,'Squid’s Short-reads')
    text_mm(figure_s8,54.5,55,'Isoform rank','axis',ha='center')
    text_mm(figure_s8,93.5,248.8,'Distinct protein rank','axis',ha='center')
    for center in [34,156.4]:
        text_mm(figure_s8,3.5,center,'Relative expression [%]',
                'axis',rotation=90,ha='center',va='center')
    assert s8_points==len(s8_expression)==400407
    assert s8_detail_points==len(s8_details)


# %% figure="FigureS8"
if 'FigureS8' in FIGURES_TO_BUILD:
    sizes_mm = np.array([a.get_position().bounds[2:] for a in figure_s8.axes]) * [180,260]
    np.testing.assert_allclose(sizes_mm, np.tile([34,24],(21,1)), atol=1e-10)
    assert s8_expression.Rank.between(.8, 80000).all()
    assert s8_details.Rank.between(.8, 80000).all()
    x_tick_offsets_mm = []
    for axis in figure_s8.axes:
        np.testing.assert_allclose(axis.get_xlim(), [.8, 80000])
        np.testing.assert_allclose(axis.get_xticks(), [1, 100, 10000])
        # Verify actual physical tick positions, not only tick values or labels.
        ticks_px = axis.transData.transform([[rank, 1] for rank in axis.get_xticks()])[:,0]
        x_tick_offsets_mm.append(((ticks_px-axis.bbox.x0)*25.4/figure_s8.dpi).tolist())
    np.testing.assert_allclose(x_tick_offsets_mm, np.tile(x_tick_offsets_mm[0], (21,1)), atol=1e-10)
    for axis in figure_s8.axes:
        np.testing.assert_allclose(axis.get_ylim(), s8_y_limits)
        np.testing.assert_allclose(axis.get_yticks(), [.0001,.001,.01,.1,1,10])
    assert s8_all_y.min() >= s8_y_limits[0] and s8_all_y.max() <= s8_y_limits[1]
    title_attachment = check_facet_title_attachment(s8_facets)
    layout_checks = check_layout_spacing(figure_s8)
    export_figure(figure_s8,'FigureS8',{'common_x_limits': [.8, 80000], 'shared_y_limits': s8_y_limits,
        'equation_centers': s8_equation_checks, 'no_refitting_this_revision': True, 'x_tick_offsets_mm': x_tick_offsets_mm,
        'layout_clearance':layout_checks, 'title_attachment':title_attachment, 'points':s8_points,
        'exact_detail_points':s8_detail_points,'facets':21,'original_reference_layout':True, 'fits_reused_from_completed_correction':True,
        'source_indexing_preserved':True,'expression_version_identical_to_Figure4':True})


# %% [markdown] figure="FigureS9"
# ## Figure S9 — cumulative isoform expression
# Full 2-gene long-read section above the original 19-gene short-read grid.
# Canvas: 180 × 260 mm; all original data, colors and facet order retained.
# The selected reference is the assembled manuscript PDF, pages 14/15;
# loose SVG assets contain different historical curves. Shared S8/S9 data
# differ from the separately saved Figure 4 merged-cache plot.
# S9 Top-N markers now use the exact 1-based Rank=N value; the previous offset is retained in evidence.
#
# This publication section loads the fixed selected-expression handoff documented in README.md; comparison variants remain separate.
#

# %% figure="FigureS9"
if 'FigureS9' in FIGURES_TO_BUILD:
    s9_facets = []
    s9_expression,s9_meta=load_table('figure_s8_s9_expression', {'styles': [{'Platform': 'Long-reads', 'Gene': 'GRIA2', 'Facet': 0, 'Color': '#00A08B'},
            {'Platform': 'Long-reads', 'Gene': 'PCLO', 'Facet': 1, 'Color': '#AF0038'},
            {'Platform': 'Short-reads', 'Gene': 'RUSC2', 'Facet': 0, 'Color': '#B68100'},
            {'Platform': 'Short-reads', 'Gene': 'TRIM2', 'Facet': 1, 'Color': '#B2828D'},
            {'Platform': 'Short-reads', 'Gene': 'CA2D3', 'Facet': 2, 'Color': '#6C4516'},
            {'Platform': 'Short-reads', 'Gene': 'ABL', 'Facet': 3, 'Color': '#0D2A63'},
            {'Platform': 'Short-reads', 'Gene': 'DGLA', 'Facet': 4, 'Color': '#620042'},
            {'Platform': 'Short-reads', 'Gene': 'K0513', 'Facet': 5, 'Color': '#750D86'},
            {'Platform': 'Short-reads', 'Gene': 'KCNAS', 'Facet': 6, 'Color': '#FB0D0D'},
            {'Platform': 'Short-reads', 'Gene': 'ACHA4', 'Facet': 7, 'Color': '#DA16FF'},
            {'Platform': 'Short-reads', 'Gene': 'ANR17', 'Facet': 8, 'Color': '#FC0080'},
            {'Platform': 'Short-reads', 'Gene': 'TWK7', 'Facet': 9, 'Color': '#862A16'},
            {'Platform': 'Short-reads', 'Gene': 'SCN1', 'Facet': 10, 'Color': '#6C7C32'},
            {'Platform': 'Short-reads', 'Gene': 'CACB2', 'Facet': 11, 'Color': '#1616A7'},
            {'Platform': 'Short-reads', 'Gene': 'RIMS2', 'Facet': 12, 'Color': '#FB00D1'},
            {'Platform': 'Short-reads', 'Gene': 'PCLO', 'Facet': 13, 'Color': '#AF0038'},
            {'Platform': 'Short-reads', 'Gene': 'DOP1', 'Facet': 14, 'Color': '#A777F1'},
            {'Platform': 'Short-reads', 'Gene': 'IQEC1', 'Facet': 15, 'Color': '#EB663B'},
            {'Platform': 'Short-reads', 'Gene': 'CSKI1', 'Facet': 16, 'Color': '#222A2A'},
            {'Platform': 'Short-reads', 'Gene': 'MTUS2', 'Facet': 17, 'Color': '#1CA71C'},
            {'Platform': 'Short-reads', 'Gene': 'ROBO2', 'Facet': 18, 'Color': '#2E91E5'}]})
    s9_details,s9_detail_meta=load_table('figure_s9_markers', {})
    figure_s9=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS9']['canvas_mm'])/25.4)
    from manuscript_figures.batch_helpers import annotate_rank_markers
    s9_annotation_checks=[]
    s9_points=0
    s9_detail_points=0


# %% figure="FigureS9"
if 'FigureS9' in FIGURES_TO_BUILD:
    for style in s9_meta['styles']:
        platform,gene,facet=style['Platform'],style['Gene'],style['Facet']
        long_reads=platform=='Long-reads'
        column=facet if long_reads else facet%4
        row=0 if long_reads else facet//4
        left=18+column*39
        top=22 if long_reads else 68+row*36.2
        width,height=(34,24)
        axis=axes_mm(figure_s9,left,top,width,height,grid='both',hide_y=column!=0)
        data=s9_expression.loc[s9_expression.Platform.eq(platform)&s9_expression.Gene.eq(gene)]
        details=s9_details.loc[s9_details.Platform.eq(platform)&s9_details.Gene.eq(gene)]

        line,=axis.plot(data.Rank,data.CumulativeExpression,color=style['Color'],linewidth=.6)
        check_xy(line,data.Rank,data.CumulativeExpression)
        points=axis.scatter(details.RankLabel,details.CumulativeExpression,s=14,marker='s',
            color='black',alpha=.5,linewidths=0,zorder=3)
        check_xy(points,details.RankLabel,details.CumulativeExpression,'scatter')
        assert (details.SourceArrayIndex==details.RankLabel-1).all()
        np.testing.assert_array_equal(details.CumulativeExpression,data.set_index('Rank').loc[details.RankLabel,'CumulativeExpression'])
        axis.set(xlim=(-.025*len(data),1.025*len(data)),ylim=(-2,105))
        axis.set_yticks([0,20,40,60,80,100])
        step=500 if len(data)<2000 else 2000 if len(data)<=10005 else 5000 if len(data)<=21711 else 10000 if len(data)<50000 else 20000
        axis.set_xticks(np.arange(0,len(data)+1,step))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x/1000:g}k' if x>=1000 else f'{x:g}'))
        s9_annotation_checks.extend(annotate_rank_markers(axis,details,line,platform,gene))

        axis.tick_params(which='both',left=False,right=False,bottom=False,top=False)
        s9_facets.append(close_facet_title(axis, gene))
        s9_points+=len(data);s9_detail_points+=len(details)
    text_mm(figure_s9,18,14,'Squid’s Long-reads')
    text_mm(figure_s9,18,60,'Squid’s Short-reads')
    text_mm(figure_s9,54.5,55,'Isoform rank','axis',ha='center')
    text_mm(figure_s9,93.5,248.8,'Distinct protein rank','axis',ha='center')
    for center in [34,156.4]:
        text_mm(figure_s9,3.5,center,'Cumulative relative expression [%]',
                'axis',rotation=90,ha='center',va='center')
    assert s9_points==len(s9_expression)==400407
    assert s9_detail_points==len(s9_details)


# %% figure="FigureS9"
if 'FigureS9' in FIGURES_TO_BUILD:
    sizes_mm = np.array([a.get_position().bounds[2:] for a in figure_s9.axes]) * [180,260]
    np.testing.assert_allclose(sizes_mm, np.tile([34,24],(21,1)), atol=1e-10)
    title_attachment = check_facet_title_attachment(s9_facets)
    layout_checks = check_layout_spacing(figure_s9)
    export_figure(figure_s9,'FigureS9',{'layout_clearance':layout_checks, 'title_attachment':title_attachment, 'points':s9_points,
        'exact_detail_points':s9_detail_points,'facets':21,'original_reference_PDF':True,
        'exact_one_based_top_n':True, 'annotations':s9_annotation_checks,'expression_version_distinct_from_Figure4':True})


# %% [markdown] figure="FigureS10"
# ## Figure S10 — subsampling and UMI deduplication
# Original a–c assembly, 180 × 250 mm. Panels a/b compare the original PacBio
# full and 10% data; c compares PacBio3 reads with recognizable barcodes with
# their UMI-deduplicated counterpart. These are distinct cohorts.
#

# %% figure="FigureS10"
if 'FigureS10' in FIGURES_TO_BUILD:
    s10_paired,s10_meta=load_table('figure_s10_paired', {'styles': [{'Panel': 'a', 'Gene': 'GRIA2', 'Facet': 0, 'Color': '#00A08B'},
            {'Panel': 'a', 'Gene': 'PCLO', 'Facet': 1, 'Color': '#AF0038'},
            {'Panel': 'b',
             'Gene': 'GRIA2',
             'Facet': 0,
             'Group': 'Subsampled data (10%)',
             'Color': 'rgb(4, 180, 161)',
             'Alpha': 0.7},
            {'Panel': 'b',
             'Gene': 'PCLO',
             'Facet': 1,
             'Group': 'Subsampled data (10%)',
             'Color': 'rgb(179, 20, 78)',
             'Alpha': 0.7},
            {'Panel': 'b',
             'Gene': 'GRIA2',
             'Facet': 0,
             'Group': 'Full data',
             'Color': 'rgb(16, 240, 227)',
             'Alpha': 0.7},
            {'Panel': 'b',
             'Gene': 'PCLO',
             'Facet': 1,
             'Group': 'Full data',
             'Color': 'rgb(191, 80, 144)',
             'Alpha': 0.7},
            {'Panel': 'c',
             'Gene': 'GRIA2',
             'Facet': 0,
             'Group': 'Full data',
             'Color': 'rgb(4, 180, 161)',
             'Alpha': 0.5},
            {'Panel': 'c',
             'Gene': 'ADAR1',
             'Facet': 1,
             'Group': 'Full data',
             'Color': 'rgb(77, 8, 229)',
             'Alpha': 0.5},
            {'Panel': 'c',
             'Gene': 'IQEC1',
             'Facet': 2,
             'Group': 'Full data',
             'Color': 'rgb(239, 122, 81)',
             'Alpha': 0.5},
            {'Panel': 'c',
             'Gene': 'GRIA2',
             'Facet': 0,
             'Group': 'Deduplicated data (UMIs)',
             'Color': 'rgb(16, 240, 227)',
             'Alpha': 0.5},
            {'Panel': 'c',
             'Gene': 'ADAR1',
             'Facet': 1,
             'Group': 'Deduplicated data (UMIs)',
             'Color': 'rgb(97, 108, 255)',
             'Alpha': 0.5},
            {'Panel': 'c',
             'Gene': 'IQEC1',
             'Facet': 2,
             'Group': 'Deduplicated data (UMIs)',
             'Color': 'rgb(251, 182, 147)',
             'Alpha': 0.5}],
 'a_limits': [0, 0.506815955603998]})
    s10_fits,s10_fits_meta=load_table('figure_s10_fits', {})
    s10_cumulative,s10_cumulative_meta=load_table('figure_s10_cumulative', {})
    figure_s10=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS10']['canvas_mm'])/25.4)
    for style in s10_meta['styles']:
        if style['Color'].startswith('rgb('):
            style['Color']=tuple(float(v)/255 for v in style['Color'][4:-1].split(','))
    s10_points=0


# %% figure="FigureS10"
if 'FigureS10' in FIGURES_TO_BUILD:
    text_mm(figure_s10,2,12,'a','panel')
    s10_equations = []
    s10_equation_artists = []
    s10_a_height = 69 * .4 / np.ptp(s10_meta['a_limits'])
    assert s10_paired['%RelativeExpression1.0'].max() < .4 and s10_fits.Y.max() < .4
    for style in [s for s in s10_meta['styles'] if s['Panel']=='a']:
        x=19+style['Facet']*83
        axis=axes_mm(figure_s10,x,16,69,s10_a_height,grid='both',hide_y=style['Facet']!=0)
        data=s10_paired.loc[s10_paired.Gene.eq(style['Gene'])]
        fits=s10_fits.loc[s10_fits.Gene.eq(style['Gene'])]
        points=axis.scatter(data['%RelativeExpression0.1'],data['%RelativeExpression1.0'],
            s=2.5,color=style['Color'],linewidths=0,clip_on=False)
        check_xy(points,data['%RelativeExpression0.1'],data['%RelativeExpression1.0'],'scatter')
        line,=axis.plot(fits.X,fits.Y,color='gray',alpha=.5,linestyle='--',linewidth=1)
        check_xy(line,fits.X,fits.Y)
        axis.set(xlim=s10_meta['a_limits'],ylim=(0,.4),aspect='equal')
        axis.set_xticks(np.arange(0,.51,.1));axis.set_yticks(np.arange(0,.41,.1))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
        axis.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
        text_mm(figure_s10,x+34.5,12,style['Gene'],ha='center')
        # Read the equation from the saved straight line; do not refit the data.
        first, last = fits.loc[fits.X.idxmin()], fits.loc[fits.X.idxmax()]
        slope = (last.Y-first.Y)/(last.X-first.X)
        intercept = first.Y-slope*first.X
        np.testing.assert_allclose(fits.Y, slope*fits.X+intercept, rtol=0, atol=1e-14)
        equation = f'y = {slope:.2f}x'
        assert intercept < 0
        if round(intercept, 2) != 0:
            equation += f' − {abs(intercept):.2f}'
        equation_text = axis.text(.05,.25,equation,transform=axis.get_yaxis_transform(),
                                  fontsize=6,ha='left',va='center')
        s10_equation_artists.append((axis,equation_text))
        s10_equations.append({'gene':style['Gene'],'slope':slope,'intercept':intercept,
            'equation':equation,'annotation_x_axes':.05,'annotation_y_data':.25,'max_observed_y':float(data['%RelativeExpression1.0'].max()),
            'max_saved_fit_y':float(fits.Y.max()),'source':'persisted fit coordinates; no refit'})
        s10_points+=len(data)+len(fits)
    text_mm(figure_s10,95,16+s10_a_height+9,'Subsampled data (10%) relative expression [%]','axis',ha='center')
    text_mm(figure_s10,7,16+s10_a_height/2,'Full data relative expression [%]','axis',ha='center',va='center',rotation=90)
    s10_upper_artists = [*figure_s10.axes, *figure_s10.texts, *figure_s10.artists]


# %% figure="FigureS10"
if 'FigureS10' in FIGURES_TO_BUILD:
    for panel,letter_top,top,height in [('b',100,107,39),('c',166,181,32)]:
        text_mm(figure_s10,2,letter_top,panel,'panel')
        genes=['GRIA2','PCLO'] if panel=='b' else ['GRIA2','ADAR1','IQEC1']
        for facet,gene in enumerate(genes):
            left=19+facet*(83 if panel=='b' else 55)
            width=69 if panel=='b' else 45
            axis=axes_mm(figure_s10,left,top,width,height,grid='both',hide_y=facet!=0)
            for style in [s for s in s10_meta['styles'] if s['Panel']==panel and s['Gene']==gene]:
                data=s10_cumulative.loc[s10_cumulative.Panel.eq(panel)&s10_cumulative.Gene.eq(gene)&s10_cumulative.Group.eq(style['Group'])]
                points=axis.scatter(data.Rank,data.CumulativeExpression,s=.9,color=style['Color'],
                    alpha=style['Alpha'],linewidths=0,clip_on=False)
                check_xy(points,data.Rank,data.CumulativeExpression,'scatter')
                axis.scatter([],[],s=6,color=style['Color'],label=style['Group'])
                s10_points+=len(data)
            axis.set(xscale='log',xlim=(.65,22000),ylim=(-2,105))
            axis.set_xticks([1,10,100,1000,10000],['1','10','100','1k','10k'])
            axis.set_yticks([0,20,40,60,80,100])
            if panel=='b':
                # Rank-1 coordinates stay exact; move only the boundary grid to the padded edge.
                axis.xaxis.get_major_ticks()[0].gridline.set_visible(False)
                axis.axvline(.65,color='#dbe4ef',linewidth=.35,zorder=0)
            axis.xaxis.set_minor_formatter(NullFormatter())
            axis.tick_params(which='both',left=False,right=False,bottom=False,top=False)
            axis.legend(loc='upper left',
                frameon=False,handletextpad=.3,borderaxespad=.2,labelspacing=.4)
            text_mm(figure_s10,left+width/2,top-3,gene,ha='center')
        text_mm(figure_s10,95,top+height+9,'Isoform rank','axis',ha='center')
        text_mm(figure_s10,7,top+height/2,'Cumulative relative expression [%]','axis',ha='center',va='center',rotation=90)
        if panel=='b':
            s10_middle_artists = [a for a in [*figure_s10.axes, *figure_s10.texts, *figure_s10.artists]
                                 if a not in s10_upper_artists]
    text_mm(figure_s10,19,166,'Cumulative expression vs. distinct unique proteins')
    text_mm(figure_s10,19,172,'Data are based on reads with recognizable barcodes','subtitle')


# %% figure="FigureS10"
if 'FigureS10' in FIGURES_TO_BUILD:
    assert s10_points==len(s10_paired)+len(s10_fits)+len(s10_cumulative)
    sizes_mm = np.array([a.get_position().bounds[2:] for a in figure_s10.axes]) * [180,250]
    np.testing.assert_allclose(sizes_mm[:2], [[69,s10_a_height],[69,s10_a_height]], atol=1e-10)
    np.testing.assert_allclose(sizes_mm[:4,0], 69, atol=1e-10)
    assert all(a.get_aspect()==1 and a.get_ylim()==(0,.4) for a in figure_s10.axes[:2])
    s10_lower_artists = [a for a in [*figure_s10.axes, *figure_s10.texts, *figure_s10.artists]
                         if a not in s10_upper_artists+s10_middle_artists]
    s10_panel_spacing = balance_middle_panel(figure_s10, s10_upper_artists, s10_middle_artists, s10_lower_artists)
    renderer = figure_s10.canvas.get_renderer()
    for axis, label in s10_equation_artists:
        bounds = label.get_window_extent(renderer)
        np.testing.assert_allclose(axis.transData.inverted().transform(bounds.get_points().mean(axis=0))[1], .25, atol=1e-12)
    layout_checks = check_layout_spacing(figure_s10)
    export_figure(figure_s10,'FigureS10',{'panel_gap_balance':s10_panel_spacing, 'y_title_x_mm':7, 'saved_regression_equations':s10_equations, 'equal_xy_unit_scale':True, 'layout_clearance':layout_checks, 'points':s10_points,'paired_protein_IDs_verified':True,
        'cumulative_ranks_and_values_exact':True,'regressions_preserved_without_refitting':True,
        'distinct_recognizable_barcode_cohort_in_c':True})


# %% [markdown] figure="FigureS11"
# ## Figure S11 — GRIA2 editing-site interdependencies
# Complete a–b assembly, 180 × 255 mm. Both native vector matrices retain
# 109 ordered sites, revised BY/lower-triangle display masks and exact colors.
#

# %% figure="FigureS11"
if 'FigureS11' in FIGURES_TO_BUILD:
    from manuscript_figures.batch_helpers import continuous_colorbar
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LinearSegmentedColormap,Normalize
    from matplotlib.cm import ScalarMappable
    s11_table,s11_meta=load_table('figure_s11_matrices.BY', {'color_scales': {'Pearson': [[0, '#67001f'],
                              [0.1, '#b2182b'],
                              [0.2, '#d6604d'],
                              [0.3, '#f4a582'],
                              [0.4, '#fddbc7'],
                              [0.5, '#f7f7f7'],
                              [0.6, '#d1e5f0'],
                              [0.7, '#92c5de'],
                              [0.8, '#4393c3'],
                              [0.9, '#2166ac'],
                              [1, '#053061']],
                  'MI': [[0, '#f7f7f7'],
                         [0.25, '#d1e5f0'],
                         [0.5, '#67a9cf'],
                         [0.75, '#2166ac'],
                         [1, '#053061']]},
 'checks': [{'reference_image_size': [374, 374]}, {'reference_image_size': [374, 374]}]})
    # Same correction as Figure 5a; the MI rows remain unchanged.
    figure_s11=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS11']['canvas_mm'])/25.4)
    s11_finite_cells=0

    s11_meta['site_order'] = s11_table.loc[s11_table.ColumnIndex.eq(0), 'RowSite'].drop_duplicates().tolist()


# %% figure="FigureS11"
if 'FigureS11' in FIGURES_TO_BUILD:
    text_mm(figure_s11,90,8,'Editing-site interdependencies in GRIA2 (Long-reads)',ha='center')
    for row,measure in enumerate(['Pearson','MI']):
        top=[18,136][row]
        data=s11_table.loc[s11_table.Measure.eq(measure)]
        n=len(s11_meta['site_order'])
        values=data.Value.to_numpy().reshape(n,n)
        mask=data.Masked.to_numpy().reshape(n,n)
        np.testing.assert_array_equal(mask,np.isnan(values))
        axis=axes_mm(figure_s11,30,top,105,105,grid=None)
        axis.set(xlim=(0,n),ylim=(n,0),xticks=[],yticks=[])
        scale=s11_meta['color_scales'][measure]
        stops=np.array([s[0] for s in scale])
        channels=np.array([[int(s[1][j:j+2],16) for j in [1,3,5]] for s in scale])
        normalized=(values[~mask]+1)/2 if measure=='Pearson' else values[~mask]
        colors=np.round(np.column_stack([np.interp(normalized,stops,channels[:,j]) for j in range(3)]))/255
        reference_height,reference_width=s11_meta['checks'][row]['reference_image_size']
        gx,gy=n/reference_width,n/reference_height
        rows,columns=np.where(~mask)
        cells=[Rectangle((x+gx/2,y+gy/2),1-gx,1-gy) for y,x in zip(rows,columns)]
        collection=PatchCollection(cells,facecolors=colors,edgecolors='none',antialiased=False)
        collection.set_gid(f'FigureS11-{measure}')
        axis.add_collection(collection)
        np.testing.assert_array_equal(collection.get_facecolors()[:,:3],colors)
        assert len(collection.get_paths())==len(cells)
        s11_finite_cells+=len(cells)
        text_mm(figure_s11,7,top-3,['a','b'][row],'panel')
        text_mm(figure_s11,85,top-.6,['Site-site correlations','Normalized mutual information'][row],ha='center',va='bottom')
        color_axis=axes_mm(figure_s11,143,top+15,4,75,grid=None)
        continuous_colorbar(color_axis,scale,-1 if measure=='Pearson' else 0,1,
            np.arange(-1,1.01,.25) if measure=='Pearson' else np.arange(0,1.01,.2))


# %% figure="FigureS11"
if 'FigureS11' in FIGURES_TO_BUILD:
    assert s11_finite_cells==s11_table.Value.notna().sum()
    layout_checks = check_layout_spacing(figure_s11)
    export_figure(figure_s11,'FigureS11',{'layout_clearance':layout_checks, 'native_heatmap_cells':s11_finite_cells,'ordered_sites':109,
        'matrix_values_masks_exact':True,'original_reference_colors_exact':True})


# %% [markdown] figure="FigureS12"
# ## Figure S12 — significant site correlations by distance
# Complete two-gene long-read section above the original 19-gene short-read
# grid. Canvas: 180 × 260 mm. Reuse persisted FDR flags, 100 bp bins, exact
# site-pair correlations, Plotly quartiles and original outlier definitions.
# Reference QA accounts separately for Plotly's one-pixel median-ink clamp;
# the median values themselves remain exact scientific medians.
#

# %% figure="FigureS12"
if 'FigureS12' in FIGURES_TO_BUILD:
    from matplotlib.patches import Rectangle
    s12_boxes,s12_meta=load_table('figure_s12_boxes', {'sources': [{'Platform': 'Long-reads',
              'gene_order': ['GRIA2', 'PCLO'],
              'affine_y_to_svg': [-220.41458468878855, 235.00565554910213]},
             {'Platform': 'Short-reads',
              'gene_order': ['ABL',
                             'ACHA4',
                             'ANR17',
                             'CA2D3',
                             'CACB2',
                             'CSKI1',
                             'DGLA',
                             'DOP1',
                             'IQEC1',
                             'K0513',
                             'KCNAS',
                             'MTUS2',
                             'PCLO',
                             'RIMS2',
                             'ROBO2',
                             'RUSC2',
                             'SCN1',
                             'TRIM2',
                             'TWK7'],
              'affine_y_to_svg': [-84.3669061143725, 88.83159481616087]}],
 'styles': [{'Platform': 'Long-reads', 'Label': '0 - 100', 'Color': 'rgb(68, 1, 84)'},
            {'Platform': 'Long-reads', 'Label': '100 - 200', 'Color': 'rgb(67, 56, 128)'},
            {'Platform': 'Long-reads', 'Label': '200 - 300', 'Color': 'rgb(49, 104, 142)'},
            {'Platform': 'Long-reads', 'Label': '300 - 400', 'Color': 'rgb(34, 144, 140)'},
            {'Platform': 'Long-reads', 'Label': '400 - 500', 'Color': 'rgb(53, 183, 121)'},
            {'Platform': 'Long-reads', 'Label': '500 - 600', 'Color': 'rgb(146, 214, 65)'},
            {'Platform': 'Long-reads', 'Label': '600 - 700', 'Color': 'rgb(253, 231, 37)'},
            {'Platform': 'Short-reads', 'Label': '0 - 100', 'Color': 'rgb(68, 1, 84)'},
            {'Platform': 'Short-reads', 'Label': '100 - 200', 'Color': 'rgb(72, 40, 120)'},
            {'Platform': 'Short-reads', 'Label': '200 - 300', 'Color': 'rgb(62, 73, 137)'},
            {'Platform': 'Short-reads', 'Label': '300 - 400', 'Color': 'rgb(49, 104, 142)'},
            {'Platform': 'Short-reads', 'Label': '400 - 500', 'Color': 'rgb(38, 130, 142)'},
            {'Platform': 'Short-reads', 'Label': '500 - 600', 'Color': 'rgb(31, 158, 137)'},
            {'Platform': 'Short-reads', 'Label': '600 - 700', 'Color': 'rgb(53, 183, 121)'},
            {'Platform': 'Short-reads', 'Label': '700 - 800', 'Color': 'rgb(110, 206, 88)'},
            {'Platform': 'Short-reads', 'Label': '800 - 900', 'Color': 'rgb(181, 222, 43)'},
            {'Platform': 'Short-reads', 'Label': '900 - 1000', 'Color': 'rgb(253, 231, 37)'}]})
    s12_outliers,s12_outlier_meta=load_table('figure_s12_outliers', {})
    s12_pairs,s12_pair_meta=load_table('figure_s12_pairs', {})
    figure_s12=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS12']['canvas_mm'])/25.4)
    s12_box_count=0
    s12_point_count=0


# %% figure="FigureS12"
if 'FigureS12' in FIGURES_TO_BUILD:
    for source in s12_meta['sources']:
        platform=source['Platform']
        long_reads=platform=='Long-reads'
        for facet,gene in enumerate(source['gene_order']):
            column=facet if long_reads else facet%5
            row=0 if long_reads else facet//5
            left=18+column*(72 if long_reads else 27)
            top=22 if long_reads else 88+row*40
            width,height=(60,43) if long_reads else (23,23)
            axis=axes_mm(figure_s12,left,top,width,height,hide_y=column!=0)
            data=s12_boxes.loc[s12_boxes.Platform.eq(platform)&s12_boxes.Gene.eq(gene)]
            ncolors=sum(s['Platform']==platform for s in s12_meta['styles'])
            for box in data.itertuples():
                color=tuple(float(v)/255 for v in box.Color[4:-1].split(','))
                stats={'q1':box.q1,'med':box.median,'q3':box.q3,'whislo':box.whislo,'whishi':box.whishi}
                artists=axis.bxp([stats],positions=[box.DistanceBin],widths=49/ncolors,
                    showfliers=False,manage_ticks=False,patch_artist=True,
                    boxprops={'facecolor':color,'edgecolor':color,'linewidth':.3,'alpha':.5},
                    whiskerprops={'color':color,'linewidth':.3},capprops={'color':color,'linewidth':.3},
                    medianprops={'color':color,'linewidth':.3})
                np.testing.assert_array_equal(artists['medians'][0].get_ydata(),[box.median,box.median])
                np.testing.assert_array_equal(artists['whiskers'][0].get_ydata(),[box.q1,box.whislo])
                np.testing.assert_array_equal(artists['whiskers'][1].get_ydata(),[box.q3,box.whishi])
                s12_box_count+=1
            for color,data in s12_outliers.loc[s12_outliers.Platform.eq(platform)&s12_outliers.Gene.eq(gene)].groupby('Color',sort=False):
                rgb=tuple(float(v)/255 for v in color[4:-1].split(','))
                points=axis.scatter(data.DistanceBin,data.r,s=1.7,color=rgb,linewidths=0,clip_on=False)
                check_xy(points,data.DistanceBin,data.r,'scatter');s12_point_count+=len(data)
            axis.set(xlim=(20,3050 if long_reads else 1050),ylim=(-.52,1.04) if long_reads else (-.44,1.04))
            axis.set_xticks([500,1000,1500,2000,2500,3000] if long_reads else [200,600,1000])
            axis.set_yticks(np.arange(-.5,1.01,.1) if long_reads else np.arange(-.4,1.01,.2))
            axis.yaxis.set_major_formatter(FuncFormatter(lambda v,p:'0' if abs(v)<1e-12 else f'{v:.1f}'.removesuffix('.0')))
            text_mm(figure_s12,left+width/2,top-3,gene,ha='center')
        heading_top=14 if long_reads else 80
        text_mm(figure_s12,18,heading_top,'Squid’s '+platform)
        text_mm(figure_s12,84,75 if long_reads else 243,'Absolute distance between site pairs [bp]','axis',ha='center')
        text_mm(figure_s12,4,43.5 if long_reads else 159.5,'Pearson’s r','axis',rotation=90,ha='center',va='center')
        legend_top=25 if long_reads else 91
        text_mm(figure_s12,157,legend_top,'Pairs','legend_heading')
        for index,style in enumerate(s for s in s12_meta['styles'] if s['Platform']==platform):
            color=tuple(float(v)/255 for v in style['Color'][4:-1].split(','))
            y=legend_top+4.5+index*3.6
            figure_s12.add_artist(Rectangle((157/180,1-(y+1)/260),1.8/180,1.8/260,
                transform=figure_s12.transFigure,facecolor=color,edgecolor='none'))
            text_mm(figure_s12,160,y,style['Label'],'legend',va='center')


# %% figure="FigureS12"
if 'FigureS12' in FIGURES_TO_BUILD:
    assert s12_box_count==len(s12_boxes) and s12_point_count==len(s12_outliers)
    assert s12_boxes.Pairs.sum()==len(s12_pairs)==34597
    export_figure(figure_s12,'FigureS12',{'exact_site_pairs':len(s12_pairs),'box_summaries':s12_box_count,
        'exact_outlier_points':s12_point_count,'original_FDR_flags_and_100bp_bins':True,
        'original_quartile_definition':True,'legacy_reference_median_ink_clamp_explained':True})


# %% [markdown] figure="FigureS13"
# ## Figure S13 — nine strongly diversified octopus transcripts
# Original 3 × 3 order and tissue legend, 180 × 150 mm. Each selected protein
# contributes once per sample. Coverage counts unique reads at eligible CDS
# positions. All 63 reference points and their colors match the cache recovery.
#

# %% figure="FigureS13"
if 'FigureS13' in FIGURES_TO_BUILD:
    s13_data,s13_meta=load_table('figure_s13_plot', {'legend': [{'Tissue': 'Axial nerve cord', 'Color': 'rgb(51, 102, 204)'},
            {'Tissue': 'Frontal & vertical lobe', 'Color': 'rgb(255, 153, 0)'},
            {'Tissue': 'Pedunculate & olfactory lobe', 'Color': 'rgb(16, 150, 24)'},
            {'Tissue': 'Stellate g. & visceral g.', 'Color': 'rgb(153, 0, 153)'},
            {'Tissue': 'Sucker', 'Color': 'rgb(0, 153, 198)'},
            {'Tissue': 'Retina & optic lobe', 'Color': 'rgb(221, 68, 119)'},
            {'Tissue': 'Non-neuronal tissues mix', 'Color': 'rgb(220, 57, 18)'}]})
    # The requested <=200 check failed: 12 coverage values exceed 200 (maximum 254).
    assert s13_data.NumOfReads.max()==254 and s13_data.NumOfProteinsPerGenePerSample.max()==171
    assert s13_data.NumOfReads.le(260).all() and s13_data.NumOfProteinsPerGenePerSample.le(200).all()
    figure_s13=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS13']['canvas_mm'])/25.4)
    s13_points=0

    s13_meta['gene_order'] = s13_data.sort_values('Facet')[['Chrom', 'Transcript', 'PooledIsoforms']].drop_duplicates().to_dict('records')


# %% figure="FigureS13"
if 'FigureS13' in FIGURES_TO_BUILD:
    for facet,gene in enumerate(s13_meta['gene_order']):
        row,column=divmod(facet,3)
        left,top=18+column*38,22+row*43
        axis=axes_mm(figure_s13,left,top,31,30,grid='both',hide_y=column!=0)
        axis.tick_params(axis='both', pad=3)  # Approved S13 tick-label spacing.
        data=s13_data.loc[s13_data.Facet.eq(facet)]
        for point in data.itertuples():
            color=tuple(float(v)/255 for v in point.Color[4:-1].split(','))
            artist=axis.scatter([point.NumOfReads],[point.NumOfProteinsPerGenePerSample],
                s=11,color=color,alpha=.7,linewidths=0,clip_on=False)
            check_xy(artist,[point.NumOfReads],[point.NumOfProteinsPerGenePerSample],'scatter')
            s13_points+=1
        axis.set(xlim=(-5,260),ylim=(-5,205))
        axis.set_xticks([0,50,100,150,200,250]);axis.set_yticks([0,50,100,150,200])
        label=gene['Transcript'].split('_')[0]+f" ({gene['PooledIsoforms']})"
        text_mm(figure_s13,left+15.5,top-2,label,ha='center')
    text_mm(figure_s13,71.5,147,'Coverage','axis',ha='center')
    text_mm(figure_s13,4,80,'Distinct protein isoforms','axis',ha='center',va='center',rotation=90)
    text_mm(figure_s13,133,23,'Tissue','legend_heading')
    for index,style in enumerate(s13_meta['legend']):
        y=28+index*4.5
        color=tuple(float(v)/255 for v in style['Color'][4:-1].split(','))
        figure_s13.add_artist(Line2D([134/180],[1-y/150],transform=figure_s13.transFigure,
            marker='o',markersize=3,color=color,alpha=.7,markeredgewidth=0,linestyle='none'))
        text_mm(figure_s13,137,y,style['Tissue'],'legend',va='center')


# %% figure="FigureS13"
if 'FigureS13' in FIGURES_TO_BUILD:
    assert s13_points==len(s13_data)==63
    export_figure(figure_s13,'FigureS13',{'exact_points':s13_points,'transcripts':9,'samples':7,
        'selected_protein_membership_and_coverage_read_membership_verified':True,
        'all_reference_coordinates_and_tissue_colors_verified':True})


# %% [markdown] figure="FigureS14"
# ## Figure S14 — neural versus non-neural pooled octopus genes
# Corrected classification on the unchanged 1,711-gene whole-transcriptome cohort.
# The two curves include 178 neural and 1,460 non-neural genes; 73 unclassified genes
# remain in the cohort table and are excluded from the curves and test. The pooled analysis notebook produces the corrected annotations, fixed-cohort
# curves and Mann–Whitney result from native processed inputs. The original 180 × 112 mm canvas and styling are retained.
#

# %% figure="FigureS14"
if 'FigureS14' in FIGURES_TO_BUILD:
    s14_input_dir = FIGURE_SOURCES / 'FigureS14'
    s14_data = pd.read_table(s14_input_dir / 'figure_s14_curves.tsv', float_precision='round_trip')
    s14_cohort = pd.read_table(s14_input_dir / 'figure_s14_cohorts.tsv', float_precision='round_trip')
    s14_test = pd.read_table(s14_input_dir / 'figure_s14_test.tsv', float_precision='round_trip')
    assert len(s14_cohort) == 1711 and s14_cohort.Chrom.is_unique
    assert s14_cohort.IsNeural.value_counts().to_dict() == {'Yes': 178, 'No': 1460, 'Missing': 73}
    assert s14_data.groupby('IsNeural').Denominator.unique().map(list).to_dict() == {'No': [1460], 'Yes': [178]}
    assert len(s14_test) == 1 and s14_test[['NeuralGenes', 'NonNeuralGenes']].iloc[0].tolist() == [178, 1460]
    s14_data['Label'] = s14_data.IsNeural.map({'Yes': 'Neural', 'No': 'Non-neural'})
    s14_data['Color'] = s14_data.IsNeural.map({'Yes': 'red', 'No': '#00aaff'})
    assert s14_data[['Label', 'Color']].notna().all().all()
    s14_meta = {'denominators': {'Yes': 178, 'No': 1460, 'Missing': 73},
        'y_limits': [0.05130229287703567, 158.48931924611142],
        'mannwhitneyu': {'statistic': float(s14_test.Statistic.iloc[0]),
            'pvalue': float(s14_test.PValue.iloc[0]), 'display': f'p-val = {s14_test.PValue.iloc[0]:.2e}'}}
    figure_s14 = plt.figure(figsize=np.array([180.0, 112.0]) / 25.4)


# %% figure="FigureS14"
if 'FigureS14' in FIGURES_TO_BUILD:
    axis=axes_mm(figure_s14,18,22,120,75,grid='both')
    for label,data in s14_data.groupby('Label',sort=False):
        color=data.Color.iloc[0]
        line,=axis.plot(data.NumOfProteins,data.CummulativeTranscripts,color=color,
            marker='o',markersize=1.6,markeredgewidth=0,linestyle='--',linewidth=.35,label=label)
        check_xy(line,data.NumOfProteins,data.CummulativeTranscripts)
    axis.set(xscale='log',yscale='log',xlim=(.7,1250),ylim=s14_meta['y_limits'],
             xlabel='Distinct protein isoforms per gene',ylabel='% of genes')
    axis.set_xticks([1,2,5,10,20,50,100,200,500,1000],['1','2','5','10','20','50','100','200','500','1000'])
    axis.set_yticks([.1,.2,.5,1,2,5,10,20,50,100],['0.1','0.2','0.5','1','2','5','10','20','50','100'])
    axis.xaxis.set_minor_formatter(NullFormatter());axis.yaxis.set_minor_formatter(NullFormatter())
    axis.tick_params(which='both',left=False,right=False,bottom=False,top=False)
    axis.text(10,1.3,s14_meta['mannwhitneyu']['display'],ha='center',fontsize=6,weight='bold')
    axis.text(10,.9,'(Mann-Whitney U)',ha='center',fontsize=6)
    axis.legend(title='Genes',loc='upper left',bbox_to_anchor=(1.06,1),frameon=False,
                borderaxespad=0,handlelength=1.5,handletextpad=.5)
    text_mm(figure_s14,18,11,'Whole-transcriptome octopus data')


# %% figure="FigureS14"
if 'FigureS14' in FIGURES_TO_BUILD:
    # Each point is the percentage of its own group with at least x isoforms.
    assert set(s14_data.IsNeural) == {'Yes', 'No'}
    for condition, curve in s14_data.groupby('IsNeural', sort=False):
        values = s14_cohort.loc[s14_cohort.IsNeural.eq(condition), 'NumOfProteins']
        assert curve.Denominator.eq(len(values)).all()
        np.testing.assert_array_equal(curve.NumOfProteins, np.sort(values.unique()))
        expected = [100 * values.ge(count).sum() / len(values) for count in curve.NumOfProteins]
        np.testing.assert_allclose(curve.CummulativeTranscripts, expected, rtol=0, atol=1e-12)
    export_figure(figure_s14, 'FigureS14', {
        'exact_curve_points': len(s14_data), 'denominators': s14_meta['denominators'],
        'corrected_curve_definition_verified': True, 'input_source': 'native pooled cohort and corrected neural annotations',
        'MannWhitney_result_source': 'Pooled analysis notebook: figure_s14_test.tsv; renderer does not test',
        'MannWhitney_U': s14_meta['mannwhitneyu']['statistic'],
        'MannWhitney_P': s14_meta['mannwhitneyu']['pvalue']})


# %% [markdown] figure="FigureS16"
# ## Figure S16 — empirical graph assessment
# Original a–b assembly, 180 × 235 mm. Recover saved mean/SD summaries for
# three mock genes, three NA probabilities and five coverage fractions.
# No simulation or graph search is rerun. Both error-bar endpoints were
# checked against the original reference vector coordinates.
#

# %% figure="FigureS16"
if 'FigureS16' in FIGURES_TO_BUILD:
    s16_data,s16_meta=load_table('figure_s16_summaries', {'unknown_probabilities': [0.09, 0.13, 0.18]})
    figure_s16=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS16']['canvas_mm'])/25.4)
    s16_checked_points=0
    s16_dashes={'solid':'-','dash':'--','dot':':'}

    def s16_errorbar(axis,data):
        style=data.iloc[0]
        artists=axis.errorbar(data.NumOfReads,data.Mean,yerr=data.SD,
            color=style.Color,alpha=style.Opacity,linestyle=s16_dashes[style.Dash],
            marker='o',markersize=1.8,markeredgewidth=0,linewidth=1,
            elinewidth=.8,capsize=1.2,capthick=.8)
        check_xy(artists.lines[0],data.NumOfReads,data.Mean)
        expected=np.array([[[x,y-sd],[x,y+sd]] for x,y,sd in zip(data.NumOfReads,data.Mean,data.SD)])
        np.testing.assert_array_equal(artists.lines[2][0].get_segments(),expected)
        return len(data)


# %% figure="FigureS16"
if 'FigureS16' in FIGURES_TO_BUILD:
    text_mm(figure_s16,2,12,'a','panel')
    for row,probability in enumerate(s16_meta['unknown_probabilities']):
        for column,gene in enumerate(['CBPC1','MY18A','VPS8']):
            left,top=18+column*38.4,25+row*25.6
            axis=axes_mm(figure_s16,left,top,34.4,20,grid='both',hide_y=column!=0)
            selected=s16_data.loc[s16_data.Panel.eq('a')&s16_data.Gene.eq(gene)&s16_data.UnknownProbability.eq(probability)]
            for trace,data in selected.groupby('Trace',sort=False):
                s16_checked_points+=s16_errorbar(axis,data)
            axis.set(xlim=(15000,105000),ylim=(-3000,85000))
            axis.set_xticks([20000,40000,60000,80000,100000],['20k','40k','60k','80k','100k'])
            axis.set_yticks([0,20000,40000,60000,80000],['0','20k','40k','60k','80k'])
            axis.tick_params(labelbottom=row==2)
            if row==0:text_mm(figure_s16,left+17.2,21,gene,ha='center')
        text_mm(figure_s16,136,top+10,f'P[NA]={probability:g}',rotation=270,ha='center',va='center')
    text_mm(figure_s16,73.6,106,'Simulated reads','axis',ha='center')
    text_mm(figure_s16,4,60.6,'Distinct proteins','axis',ha='center',va='center',rotation=90)
    handles=[Line2D([],[],color=color,linewidth=1,alpha=.7,label=label) for color,label in
             [('black','Complete'),('grey','Errored + partially NA')]]
    figure_s16.legend(handles=handles,loc='center right',bbox_to_anchor=(140/180,1-10/235),
                      frameon=False,ncol=2,handlelength=1.5,columnspacing=1.5)


# %% figure="FigureS16"
if 'FigureS16' in FIGURES_TO_BUILD:
    text_mm(figure_s16,10,121,'b','panel')
    for row,metric,title,top,height,limits in [
        (0,'CoveredPercent','Mean % of covered isoforms',122,44,(46,86)),
        (1,'FalsePositivePercent','Mean % of false-positive isoforms',171,44,(3.8,10.4))]:
        axis=axes_mm(figure_s16,22,top,66.7,height,grid='both')
        for trace,data in s16_data.loc[s16_data.Panel.eq('b')&s16_data.Metric.eq(metric)].groupby('Trace',sort=False):
            s16_checked_points+=s16_errorbar(axis,data)
        axis.set(xlim=(15000,105000),ylim=limits,ylabel=title)
        axis.set_xticks([20000,40000,60000,80000,100000],['20k','40k','60k','80k','100k'])
        axis.set_yticks(np.arange(50,86,5) if row==0 else np.arange(4,11))
        axis.tick_params(labelbottom=row==1)
    text_mm(figure_s16,55.35,224,'Simulated reads','axis',ha='center')
    text_mm(figure_s16,99,122,'NA probability, mock gene','legend_heading',va='top')
    for group,probability in enumerate(s16_meta['unknown_probabilities']):
        top=132
        left=99+group*25
        text_mm(figure_s16,left,top,f'{probability:g}','legend_heading')
        for index,(gene,color) in enumerate([('CBPC1','red'),('MY18A','blue'),('VPS8','green')]):
            y=top+4.5+index*4
            figure_s16.add_artist(Line2D([(left)/180,(left+5)/180],[1-y/235]*2,
                transform=figure_s16.transFigure,color=color,alpha=.7,linewidth=1,
                linestyle=['-','--',':'][group],marker='o',markevery=[1],markersize=1.8,markeredgewidth=0))
            text_mm(figure_s16,left+7,y,gene,'legend',va='center')


# %% figure="FigureS16"
if 'FigureS16' in FIGURES_TO_BUILD:
    assert s16_checked_points==len(s16_data)==180
    layout_checks = check_layout_spacing(figure_s16)
    export_figure(figure_s16,'FigureS16',{'layout_clearance':layout_checks, 'exact_means':180,'exact_SD_intervals':180,
        'original_genes_probabilities_fractions_colors_and_dashes':True,'simulations_rerun':False})


# %% [markdown] figure="FigureS17"
# ## Figure S17 — editing in common and rare isoforms
# Complete 180 × 260 mm assembly: two long-read facets above 19 short-read
# facets. Exact saved observations and native statistical summaries, with
# original vector violin contours; no KDE or biological analysis rerun.
#

# %% figure="FigureS17"
if 'FigureS17' in FIGURES_TO_BUILD:
    import hashlib
    from manuscript_figures.violin_paths import original_violin_path
    from matplotlib.patches import PathPatch
    from matplotlib.transforms import Affine2D
    s17_observations,s17_meta=load_table('figure_s17_observations', {})
    s17_stats,s17_stats_meta=load_table('figure_s17_statistics', {})
    s17_geometry_path=figure_source('figure_s17_geometry', '.json')
    s17_geometry=json.loads(s17_geometry_path.read_text())
    figure_s17=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS17']['canvas_mm'])/25.4)
    s17_checked=0


# %% figure="FigureS17"
if 'FigureS17' in FIGURES_TO_BUILD:
    for platform,top,width,height,pitch,columns in [
        ('Long-reads',23,26.8,25,0,2),('Short-reads',72,26.8,25,43,5)]:
        genes=s17_stats.loc[s17_stats.Platform.eq(platform),'Gene'].drop_duplicates()
        text_mm(figure_s17,18,top-9,platform)
        for facet,gene in enumerate(genes):
            row,column=divmod(facet,columns)
            left=22+column*31
            axis=axes_mm(figure_s17,left,top+row*pitch,width,height)
            subset=s17_stats.loc[s17_stats.Platform.eq(platform)&s17_stats.Gene.eq(gene)]
            geometries=[g for g in s17_geometry if g['Platform']==platform and g['Gene']==gene]
            for summary,geometry in zip(subset.itertuples(),geometries):
                contour=original_violin_path(geometry['path'])
                transform=Affine2D().scale(2/geometry['width_px'],1/geometry['y_scale']).translate(-.5,-geometry['y_offset']/geometry['y_scale'])
                axis.add_patch(PathPatch(transform.transform_path(contour),facecolor='white',edgecolor=summary.Color,linewidth=.6,zorder=2))
                stats={'q1':summary.Q1,'med':summary.Median,'q3':summary.Q3,
                    'whislo':summary.WhiskerLow,'whishi':summary.WhiskerHigh,
                    'mean':summary.Mean,'fliers':geometry['outliers']}
                artists=axis.bxp([stats],positions=[summary.Category],widths=.1224,
                    showmeans=True,meanline=True,showcaps=False,manage_ticks=False,
                    boxprops={'color':summary.Color,'linewidth':.5},
                    medianprops={'color':summary.Color,'linewidth':.5},
                    whiskerprops={'color':summary.Color,'linewidth':.5},
                    meanprops={'color':summary.Color,'linewidth':.5,'linestyle':'--'},
                    flierprops={'marker':'.','markersize':1.5,'markeredgewidth':0,'markerfacecolor':summary.Color},zorder=3)
                np.testing.assert_array_equal(artists['medians'][0].get_ydata(),[summary.Median]*2)
                np.testing.assert_array_equal(artists['means'][0].get_ydata(),[summary.Mean]*2)
                np.testing.assert_array_equal(artists['fliers'][0].get_ydata(),geometry['outliers'])
                values=s17_observations.loc[s17_observations.Platform.eq(platform)&s17_observations.Gene.eq(gene)&s17_observations.Abundancy.eq(summary.Abundancy),'MeanEditedPositions']
                np.testing.assert_array_equal(np.quantile(values,[.25,.5,.75],method='hazen'),[summary.Q1,summary.Median,summary.Q3])
                s17_checked+=len(values)
            g=geometries[0]
            axis.set(xlim=(-.5,1.5),ylim=((g['height_px']-g['y_offset'])/g['y_scale'],-g['y_offset']/g['y_scale']),
                xticks=[0,1],xticklabels=['1%\nmost\nexpressed','30%\nleast\nexpressed'])
            axis.set_yticks([0,20,40,60] if platform=='Long-reads' else [0,10,20,30,40])
            if column: axis.tick_params(labelleft=False)
            axis.set_title(gene,fontsize=7,pad=4)
    text_mm(figure_s17,5,35.5,'Mean edited positions in\noriginally-supporting reads',rotation=90,ha='center',va='center')
    text_mm(figure_s17,5,149,'Mean edited positions in originally-supporting reads',rotation=90,ha='center',va='center')


# %% figure="FigureS17"
if 'FigureS17' in FIGURES_TO_BUILD:
    assert s17_checked==len(s17_observations)
    sizes_mm = np.array([a.get_position().bounds[2:] for a in figure_s17.axes]) * [180,260]
    np.testing.assert_allclose(sizes_mm, np.tile([26.8,25],(21,1)), atol=1e-10)
    layout_checks = check_layout_spacing(figure_s17)
    export_figure(figure_s17,'FigureS17',{'layout_clearance':layout_checks, 'observations':s17_checked,'violin_groups':42,'native_contours':True,'exact_statistics':True})


# %% [markdown] figure="FigureS18"
# ## Figure S18 — editing at each site in common and rare isoforms
#
# PacBio and Illumina bars now show original-read G/(A+G) percentages.
# Two-sided Fisher tests and within-gene BH correction determine significance;
# the previous isoform SEM bars have been removed (2026-11-09). The author approved smaller facet text for S18 only, to display every site without
# label overlap on one 180 × 260 mm page. The source arrangement is retained (two full-width PacBio rows, then three columns of Illumina genes).
#

# %% figure="FigureS18"
if 'FigureS18' in FIGURES_TO_BUILD:
    # Read only revised read-level means and Fisher/BH results.
    s18_table = pd.read_table(figure_source('figure_s18_fisher'), float_precision='round_trip')
    assert s18_table.groupby(['Platform', 'Gene']).ngroups == 21
    assert not s18_table.duplicated(['Platform', 'Gene', 'EditingSite']).any()
    s18_width, s18_height = 180, 260  # Author-approved single-page maximum; facet text may shrink.
    figure_s18 = plt.figure(figsize=np.array([s18_width, s18_height]) / 25.4)
    s18_colors = {'Common': ['#238bca', '#0f275c'], 'Rare': ['#ff9324', '#ffcf66']}
    s18_axes, s18_checks = [], []
    # Move the complete Short-reads group up 6 mm; keep all facet sizes/row pitches.
    s18_group_top_mm = {'Long-reads': 18, 'Short-reads': 83}
    s18_headers = {}
    for letter, platform in zip('ab', s18_group_top_mm, strict=True):
        top = s18_group_top_mm[platform] - 7.5
        s18_headers[platform] = (
            text_mm(figure_s18, 3, top, letter, 'panel'),
            text_mm(figure_s18, 16, top, platform))


# %% figure="FigureS18"
if 'FigureS18' in FIGURES_TO_BUILD:
    from matplotlib.patches import Patch
    for platform in ['Long-reads', 'Short-reads']:
        genes = s18_table.loc[s18_table.Platform.eq(platform), 'Gene'].drop_duplicates()
        for index, gene in enumerate(genes):
            row, column = (index, 0) if platform == 'Long-reads' else divmod(index, 3)
            left = 16 if platform == 'Long-reads' else 16 + column * 54
            top = s18_group_top_mm[platform] + row * (29 if platform == 'Long-reads' else 24)
            width = 156 if platform == 'Long-reads' else 48
            panel_height_mm = 18 if platform == 'Long-reads' else 13
            axis = axes_mm(figure_s18, left, top, width, panel_height_mm, hide_y=column > 0)
            data = s18_table.loc[s18_table.Platform.eq(platform) & s18_table.Gene.eq(gene)].sort_values('SiteOrder')
            both = data.Testable
            x = np.arange(len(data))
            for shift, group, value in [(-.2125, 'Common', 'CommonEditingPercent'),
                                       (.2125, 'Rare', 'RareEditingPercent')]:
                colors = [s18_colors[group][0 if covered else 1] for covered in both]
                bars = axis.bar(x + shift, data[value], width=.40375, color=colors, edgecolor='none', zorder=3)
                np.testing.assert_array_equal([bar.get_height() for bar in bars], data[value])
            upper = data[['CommonEditingPercent', 'RareEditingPercent']].max(axis=1)
            stars = data.Significant.fillna(False)
            star_height = upper + max(1, .03 * upper.max())
            for position, height in zip(x[stars], star_height[stars]):
                axis.text(position, height, '*', color='black', fontsize=6, ha='center', va='center')
            axis.set(xlim=(-.5, len(data)-.5), ylim=(0, 107))
            axis.set_xticks(x, data.EditingSite.astype(str), rotation=90)
            axis.set_yticks(np.arange(0, 101, 10))
            if column == 0: axis.set_ylabel('Editing [%]', labelpad=3)
            axis.tick_params(axis='y', labelsize=4 if platform == 'Long-reads' else 3)
            axis.set_title(gene, fontsize=7, y=1 + 2.1/panel_height_mm, pad=0, va='bottom')
            details = f'{len(data)} sites, {stars.sum()} ({100 * stars.sum() / len(data):.1f}%) significant'
            annotation = axis.text(.5, 1 + .25/panel_height_mm, details, fontsize=3.5,
                                   transform=axis.transAxes, ha='center', va='bottom', color='#2a3f5f')
            annotation.set_gid('s18-author-approved-compact-text')
            assert star_height.max() < 107
            s18_axes.append(axis)
            s18_checks.append({'platform': platform, 'gene': gene, 'sites': len(data), 'significant': int(stars.sum())})
    assert s18_table.Testable.all()  # Both original-read groups cover every plotted site.
    legend_handles = [Patch(color=s18_colors[group][0], label=group) for group in ['Common', 'Rare']]
    s18_legend = figure_s18.legend(handles=legend_handles, loc='lower right', bbox_to_anchor=(.97, .025), frameon=False,
                     title='Both groups covered', fontsize=6, title_fontsize=7)


# %% figure="FigureS18"
if 'FigureS18' in FIGURES_TO_BUILD:
    figure_s18.canvas.draw()
    renderer = figure_s18.canvas.get_renderer()
    for axis, check in zip(s18_axes, s18_checks, strict=True):
        labels = axis.get_xticklabels()
        slot = axis.get_window_extent(renderer).width / len(labels)
        widest = max(label.get_window_extent(renderer).width for label in labels)
        size = min(5, np.floor(10 * 6 * .82 * slot / widest) / 10)
        assert size > 0
        for artist in labels + axis.get_yticklabels() + axis.texts:
            artist.set_gid('s18-author-approved-compact-text')
        for label in labels: label.set_fontsize(size)
        for star in axis.texts:
            if star.get_text() == '*': star.set_fontsize(size)
        check['site_label_pt'] = size
        check['y_tick_pt'] = axis.yaxis.get_major_ticks()[0].label1.get_fontsize()
    figure_s18.canvas.draw()
    for axis in s18_axes:
        bounds = [label.get_window_extent(renderer) for label in axis.get_xticklabels() if label.get_visible()]
        assert all(right.x0 >= left.x1 for left, right in zip(bounds, bounds[1:])), 'S18 site labels overlap'
    # Both platform headings sit 1.5 mm above their first gene title; letters share Y.
    heading_checks = []
    for platform, (letter, heading) in s18_headers.items():
        first_axis = next(axis for axis, check in zip(s18_axes, s18_checks, strict=True)
                          if check['platform'] == platform)
        title_top = first_axis.title.get_window_extent(renderer).y1
        target_bottom = title_top + 1.5 * figure_s18.dpi / 25.4
        shift_y = (target_bottom - heading.get_window_extent(renderer).y0) / figure_s18.bbox.height
        new_y = heading.get_position()[1] + shift_y
        for artist in [letter, heading]:
            artist.set_position((artist.get_position()[0], new_y))
        gap_mm = (heading.get_window_extent(renderer).y0 - title_top) * 25.4 / figure_s18.dpi
        np.testing.assert_allclose(gap_mm, 1.5, atol=1e-8)
        assert letter.get_position()[1] == heading.get_position()[1]
        heading_checks.append({'platform': platform, 'letter_and_heading_y': new_y,
                               'heading_to_first_gene_title_mm': gap_mm})
    figure_s18.canvas.draw()
    # Align the full rendered legend box with the last Short-reads plotting row.
    last_row_bottom = min(axis.get_window_extent(renderer).y0
                          for axis, check in zip(s18_axes, s18_checks, strict=True)
                          if check['platform'] == 'Short-reads')
    old_legend_bottom = s18_legend.get_window_extent(renderer).y0
    shift = (last_row_bottom - old_legend_bottom) / figure_s18.bbox.height
    s18_legend.set_bbox_to_anchor((.97, .025 + shift), transform=figure_s18.transFigure)
    figure_s18.canvas.draw()
    legend_bottom = s18_legend.get_window_extent(renderer).y0
    np.testing.assert_allclose(legend_bottom, last_row_bottom, atol=1e-8)
    legend_alignment = {'legend_bottom_mm_from_page_bottom': legend_bottom * 25.4 / figure_s18.dpi,
        'last_row_bottom_mm_from_page_bottom': last_row_bottom * 25.4 / figure_s18.dpi,
        'legend_raised_mm': (legend_bottom - old_legend_bottom) * 25.4 / figure_s18.dpi}
    spacing_checks = check_layout_spacing(figure_s18)
    export_figure(figure_s18, 'FigureS18', {'canvas_mm': [s18_width, s18_height], 'facets': s18_checks,
        'text_exception': 'Author explicitly permits smaller facet text to show every site on one 180 x 260 mm page.',
        'spacing_checks': spacing_checks, 'legend_alignment': legend_alignment,
        'platform_heading_checks': heading_checks, 'short_reads_shift_up_mm': 6,
        'all_sites_and_gene_platforms_present': True, 'site_tick_labels_do_not_overlap': True,
        'FDR': 'Two-sided Fisher; BH over testable sites within each gene',
        'rare_unit': 'original supporting read', 'common_mean': '100*G/(A+G), original supporting reads',
        'error_bars': 'none; historical isoform SEM removed'})


# %% [markdown] figure="FigureS19"
# Figure S19 — completed from the original Julia simulation, rerun with documented
# seed 1892, N=10,000, K=10 and the original nine probability pairs. This new run uses
# one Julia thread; the historical run's thread count/RNG state were not recorded,
# so numerical identity to the old artwork is not claimed. All read matrices,
# correlation matrices, RNG checkpoints and actual 15-bin edges/counts are saved.
# The following Plotly cells read those saved tables only.

# %%
if 'FigureS19' in FIGURES_TO_BUILD:
    import plotly.graph_objects as go
    import hashlib
    from manuscript_figures.plotly_export import export_plotly
    s19_input = FIGURE_SOURCES / 'FigureS19'
    s19_histograms = pd.read_table(s19_input / 's19_histograms.tsv', float_precision='round_trip')
    s19_correlations = pd.read_table(s19_input / 's19_correlations.tsv', float_precision='round_trip')
    assert len(s19_histograms) == 135 and len(s19_correlations) == 900
    assert s19_histograms.groupby(['EditingProbability', 'ReadRemovalProbability']).Count.sum().eq(45).all()
    s19_values = s19_correlations.loc[s19_correlations.Displayed.eq(1), 'Correlation']
    assert len(s19_values) == 405
    s19_xrange = [0.15, 0.45]  # Author-approved shared range; original bins are unchanged.
    assert s19_histograms.Left.min() >= s19_xrange[0]
    assert s19_histograms.Right.max() <= s19_xrange[1]
    assert s19_histograms.ExpectedCorrelation.between(*s19_xrange).all()
    s19_yrange = [0, float(2*np.ceil(s19_histograms.Count.max()/2))]
    assert s19_histograms.Count.max() < s19_yrange[1]
    s19_pt = 4 * 25.4 / 72  # Four SVG coordinates per mm, including font sizes.
    s19_canvas_mm = 160
    figure_s19 = go.Figure(layout=dict(width=4*s19_canvas_mm, height=4*s19_canvas_mm, template='none',
        margin=dict(l=0,r=0,t=0,b=0), paper_bgcolor='white',plot_bgcolor='white',
        font=dict(family='Open Sans',size=6*s19_pt,color='#2a3f5f'), showlegend=False))

# %%
if 'FigureS19' in FIGURES_TO_BUILD:
    s19_scale = s19_canvas_mm/180
    s19_gap_mm = 9*s19_scale  # Equal 8 mm gaps horizontally and vertically.
    s19_panel_width_mm = 48*s19_scale
    s19_panel_height_mm = (137*s19_scale - 2*s19_gap_mm)/3
    s19_left_mm, s19_top_mm = 11*s19_scale, 19*s19_scale
    for row, editing in enumerate([.001,.002,.003]):
        for column, removal in enumerate([.001,.002,.003]):
            number = row*3 + column + 1
            suffix = '' if number == 1 else str(number)
            data = s19_histograms.loc[s19_histograms.EditingProbability.eq(editing)
                                     & s19_histograms.ReadRemovalProbability.eq(removal)]
            width,height = s19_panel_width_mm,s19_panel_height_mm
            left = s19_left_mm + (width+s19_gap_mm)*column
            top = s19_top_mm + (height+s19_gap_mm)*row
            xdomain = [left/s19_canvas_mm,(left+width)/s19_canvas_mm]
            ydomain = [1-(top+height)/s19_canvas_mm,1-top/s19_canvas_mm]
            common = dict(showline=False,zeroline=False,ticks='',showgrid=True,
                          gridcolor='#dbe4ef',gridwidth=.5,tickfont=dict(size=6*s19_pt))
            figure_s19.update_layout({f'xaxis{suffix}':dict(common,domain=xdomain,anchor='y'+suffix,
                range=s19_xrange,tick0=.15,dtick=.05,tickformat='.2f',ticklabelstandoff=8,showticklabels=row==2),
                f'yaxis{suffix}':dict(common,domain=ydomain,anchor='x'+suffix,
                range=s19_yrange,dtick=2,showticklabels=column==0)})
            # Complete the grid at both Y boundaries using the ordinary grid style.
            for boundary in s19_yrange:
                figure_s19.add_shape(type='line',xref='x'+suffix+' domain',yref='y'+suffix,
                    x0=0,x1=1,y0=boundary,y1=boundary,layer='below',
                    line=dict(color=common['gridcolor'],width=common['gridwidth']))
            # Complete the grid at both X boundaries; labels sit 2 mm below Y=0.
            for boundary in s19_xrange:
                figure_s19.add_shape(type='line',xref='x'+suffix,yref='y'+suffix+' domain',
                    x0=boundary,x1=boundary,y0=0,y1=1,layer='below',
                    line=dict(color=common['gridcolor'],width=common['gridwidth']))
            figure_s19.add_trace(go.Bar(x=(data.Left+data.Right)/2,y=data.Count,
                width=data.Right-data.Left,marker=dict(color='#0072b2',line_width=0),
                xaxis='x'+suffix,yaxis='y'+suffix,hoverinfo='skip'))
            expected = float(data.ExpectedCorrelation.iloc[0])
            figure_s19.add_shape(type='line',xref='x'+suffix,yref='y'+suffix,
                x0=expected,x1=expected,y0=0,y1=s19_yrange[1],line=dict(color='red',width=1.5,dash='dot'))
            occupied = data.loc[data.Count.gt(0)]
            free_intervals = [(s19_xrange[0], min(expected, occupied.Left.min())),
                              (max(expected, occupied.Right.max()), s19_xrange[1])]
            text_interval = max(free_intervals, key=lambda interval: interval[1]-interval[0])
            annotation_x = sum(text_interval)/2
            figure_s19.add_annotation(xref='x'+suffix,yref='y'+suffix,x=annotation_x,y=s19_yrange[1]/2,
                text=f'expected =<br>{expected:.3f}',xanchor='center',align='center',
                yanchor='middle',showarrow=False,font=dict(color='red',size=6*s19_pt))
            assert int(data.Count.sum()) == 45
    # Shared column and row headings replace the nine repeated panel titles.
    for column, removal in enumerate([.001,.002,.003]):
        center_x = s19_left_mm + (s19_panel_width_mm+s19_gap_mm)*column + s19_panel_width_mm/2
        figure_s19.add_annotation(xref='paper',yref='paper',x=center_x/s19_canvas_mm,y=1-(s19_top_mm-2)/s19_canvas_mm,
            text=f'Remove read = {removal*100:g}%',xanchor='center',yanchor='bottom',
            showarrow=False,font=dict(size=7*s19_pt))
    for row, editing in enumerate([.001,.002,.003]):
        center_y = s19_top_mm + (s19_panel_height_mm+s19_gap_mm)*row + s19_panel_height_mm/2
        figure_s19.add_annotation(xref='paper',yref='paper',x=1-3/s19_canvas_mm,y=1-center_y/s19_canvas_mm,
            text=f'Edit site = {editing*100:g}%',textangle=90,xanchor='center',yanchor='middle',
            showarrow=False,font=dict(size=7*s19_pt))
    figure_s19.add_annotation(xref='paper',yref='paper',x=92/180,y=14/180,text="Pearson's r",
        showarrow=False,yanchor='middle',font=dict(size=7*s19_pt))
    figure_s19.add_annotation(xref='paper',yref='paper',x=3/s19_canvas_mm,y=1-87.5/180,text='Sites',textangle=-90,
        showarrow=False,xanchor='center',font=dict(size=7*s19_pt))

# %%
if 'FigureS19' in FIGURES_TO_BUILD:
    assert len(figure_s19.data) == 9 and len(figure_s19.layout.shapes) == 45
    for trace, (_, data) in zip(figure_s19.data, s19_histograms.groupby(
            ['EditingProbability','ReadRemovalProbability'], sort=True)):
        np.testing.assert_array_equal(trace.y, data.Count)
        np.testing.assert_array_equal(trace.width, data.Right-data.Left)
    expected_lines = [line for line in figure_s19.layout.shapes if line.line.color == 'red']
    assert len(expected_lines) == 9
    for line in expected_lines:
        np.testing.assert_allclose([line.y0,line.y1],s19_yrange)
        assert line.line.width == 1.5
    export_plotly(figure_s19, 'FigureS19', s19_canvas_mm, s19_canvas_mm,
                  output_dir=FINAL_FIGURES_DIR / 'render_work/raw', font_dir=FONT_DIR)
    display(Image(filename=str(FINAL_FIGURES_DIR / 'render_work/raw/previews/FigureS19.png')))

# %% [markdown] figure="FigureS20"
# ## Figure S20 — common expression and rare chimeras
# Complete 180 × 205 mm a–b assembly, using native saved values and cached
# cohort denominators. Panel a retains the original three rows per gene
# for the three rare-fraction settings; these are not biological replicates.
# Panel b preserves the reference's displayed multi-category gene order.
#

# %% figure="FigureS20"
if 'FigureS20' in FIGURES_TO_BUILD:
    s20_expression,s20_meta=load_table('figure_s20_expression', {})
    s20_boxes,_=load_table('figure_s20_boxes', {})
    s20_chimeras,_=load_table('figure_s20_chimeras', {})
    figure_s20=plt.figure(figsize=np.array(BATCH_LAYOUT['FigureS20']['canvas_mm'])/25.4)
    s20_labels={'Illumina1':'IL','PacBio1':'PB'}


# %% figure="FigureS20"
if 'FigureS20' in FIGURES_TO_BUILD:
    text_mm(figure_s20,10,17,'a',weight='bold')
    axis=axes_mm(figure_s20,25,15,65,53.6)
    for summary in s20_boxes.itertuples():
        data=s20_expression.loc[s20_expression.Platform.eq(summary.Platform)]
        stats={'q1':summary.Q1,'med':summary.Median,'q3':summary.Q3,
            'whislo':summary.WhiskerLow,'whishi':summary.WhiskerHigh,'fliers':[]}
        artists=axis.bxp([stats],positions=[summary.Order],widths=.49,
            patch_artist=True,manage_ticks=False,showfliers=False,
            boxprops={'facecolor':mpl.colors.to_rgba(summary.Color,.5),'edgecolor':summary.Color,'linewidth':.7},
            medianprops={'color':summary.Color,'linewidth':.7},
            whiskerprops={'color':summary.Color,'linewidth':.7},capprops={'color':summary.Color,'linewidth':.7})
        points=axis.scatter(summary.Order+data.JitterOffset,data.CumulativeExpression,
            s=5,color=summary.Color,edgecolors='none',zorder=3)
        np.testing.assert_array_equal(points.get_offsets()[:,1],data.CumulativeExpression)
        np.testing.assert_array_equal(artists['medians'][0].get_ydata(),[summary.Median]*2)
    axis.set(xlim=(-.8,1.65),ylim=(18,78),xticks=[0,1],xticklabels=['Illumina','PacBio'],
        yticks=np.arange(20,80,10),ylabel='Cumulative expression of\n1% most-common proteins [%]')


# %% figure="FigureS20"
if 'FigureS20' in FIGURES_TO_BUILD:
    text_mm(figure_s20,10,93.6,'b',weight='bold')
    for column,comparison in enumerate(['Strict comparison','Soft comparison']):
        data=s20_chimeras.loc[s20_chimeras.Comparison.eq(comparison)].sort_values('DisplayOrder')
        axis=axes_mm(figure_s20,23+column*79,95.6,68,45.6)
        bars=axis.bar(data.DisplayOrder,data.ChimericPercent,width=.8,color=data.Color,linewidth=0,zorder=3)
        np.testing.assert_array_equal([bar.get_height() for bar in bars],data.ChimericPercent)
        np.testing.assert_array_equal(100*data.RecoveredChimericCount/data.ActualYRareProteins,data.ChimericPercent)
        axis.set(xlim=(-.6,20.6),ylim=(0,37.02),yticks=np.arange(0,36,5),
            xticks=data.DisplayOrder,xticklabels=data.Gene)
        axis.tick_params(axis='x',labelrotation=90)
        if column: axis.tick_params(labelleft=False)
        axis.set_title(comparison,fontsize=7,pad=5)
        for boundary in [-.5,18.5,20.5]:
            axis.plot([boundary]*2,[0,-.39],transform=axis.get_xaxis_transform(),
                color='#555555',linewidth=.45,clip_on=False)
        for platform,group in data.groupby('Platform',sort=False):
            axis.text(group.DisplayOrder.mean(),-.32,s20_labels[platform],
                transform=axis.get_xaxis_transform(),fontsize=7,ha='center',va='top')
    text_mm(figure_s20,12.5,118.4,'Chimeric rare proteins [%]',rotation=90,ha='center',va='center')


# %% figure="FigureS20"
if 'FigureS20' in FIGURES_TO_BUILD:
    layout_checks = check_layout_spacing(figure_s20)
    export_figure(figure_s20,'FigureS20',{'layout_clearance':layout_checks, 'box_observations':len(s20_expression),
        'bars':len(s20_chimeras),'original_denominators':True,'native_jitter':True})


# %% [markdown] figure="FigureS6c_gene_maxima_review"
# ## Figure S6c — approved statistical unit
# The author accepted Pearson on the 21 displayed gene-platform maxima (12 October 2026):
# r = 0.7435724549605074, p = 0.00011191897358298873. The sampled Illumina cohort is unchanged.
# The 336 algorithm rows remain evidence; the 21 points are not 21 biological replicates.
# The following legacy draft is disabled and excluded from the final collection.
#

# %% figure="FigureS6c_gene_maxima_review"
if 'FigureS6c_gene_maxima_review' in FIGURES_TO_BUILD:
    from scipy.stats import pearsonr
    s6_review_maxima,_=load_table('figure_s6c_maxima')
    s6_review_solutions,_=load_table('figure_s6c_solutions', {})
    assert not s6_review_maxima.duplicated(['Platform','Gene']).any()
    s6_review_result=pearsonr(s6_review_maxima.EditableCodons,s6_review_maxima.MaximumIsoforms)
    s6_legacy_result=pearsonr(s6_review_solutions.EditableAAs,s6_review_solutions.NumOfProteins)
    np.testing.assert_allclose(s6_review_result,[.7435724549605073,.00011191897358298873],rtol=0,atol=1e-15)
    np.testing.assert_allclose(s6_legacy_result,[.742006575020387,5.784191507143405e-60],rtol=1e-12,atol=0)
    s6_review_figure=plt.figure(figsize=np.array([180,87])/25.4)


# %% figure="FigureS6c_gene_maxima_review"
if 'FigureS6c_gene_maxima_review' in FIGURES_TO_BUILD:
    _,s6_review_meta=load_table('figure_s6c_maxima')
    text_mm(s6_review_figure,19,7,'REVIEW DRAFT — 21 gene maxima')
    axis=axes_mm(s6_review_figure,19,19,93,56,grid='both')
    for style in s6_review_maxima.itertuples():
        marker='*' if style.Marker=='star' else 'o'
        point=axis.scatter([style.EditableCodons],[style.MaximumIsoforms],marker=marker,
            s=36 if marker=='*' else 12,color=style.Color,linewidths=0,clip_on=False)
        check_xy(point,[style.EditableCodons],[style.MaximumIsoforms],'scatter')
    axis.set(xlim=s6_review_meta['original_x_limits'],ylim=s6_review_meta['original_y_limits'],
             xlabel='Editable codons',ylabel='Distinct protein isoforms')
    axis.set_xticks(np.arange(0,141,20))
    axis.set_yticks([0,5000,10000,15000,20000],['0','5k','10k','15k','20k'])
    axis.text(110,2500,f'Pearson r = {s6_review_result.statistic:.3f}\np = {s6_review_result.pvalue:.3g}\nn = 21 gene maxima',fontsize=6,ha='center',va='center',
              bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':1})
    for platform,heading_top in [('Long-reads',18),('Short-reads',35)]:
        text_mm(s6_review_figure,123,heading_top,platform,'legend_heading')
        data=s6_review_maxima.loc[s6_review_maxima.Platform.eq(platform)]
        for index,style in enumerate(data.itertuples()):
            column=index//10 if platform=='Short-reads' else 0
            row=index%10 if platform=='Short-reads' else index
            x=123+column*27
            y=heading_top+5+row*3.6
            s6_review_figure.add_artist(Line2D([x/180],[1-y/87],transform=s6_review_figure.transFigure,
                marker='*' if style.Marker=='star' else 'o',color=style.Color,
                markersize=5 if style.Marker=='star' else 3.5,markeredgewidth=0,linestyle='none'))
            text_mm(s6_review_figure,x+3,y,style.Gene,'legend',va='center')


# %% figure="FigureS6c_gene_maxima_review"
if 'FigureS6c_gene_maxima_review' in FIGURES_TO_BUILD:
    review_evidence={'review_only':True,'statistical_unit':'one maximum per gene/platform',
        'n':21,'r':s6_review_result.statistic,'p':s6_review_result.pvalue,
        'legacy_n':336,'legacy_r':s6_legacy_result.statistic,'legacy_p':s6_legacy_result.pvalue,
        'legacy_unit':'dependent algorithm/solution rows repeated within 21 genes'}
    (BATCH_DIR/'revision_2026-10-11/evidence/s6c_statistical_units.json').write_text(json.dumps(review_evidence,indent=2))
    export_figure(s6_review_figure,'FigureS6c_gene_maxima_review',review_evidence)


# %% [markdown]
# ## Final publication layout and 600 DPI exports
#
# The 14 October layout is applied at native physical scale. Internal placement changes are specified above in millimetres. The final step crops only measured exterior white space; it never scales an assembled page. Crop rectangles use PDF points (72 per inch), origin at the original page top left. Unchanged artwork is bound to the approved backup by SHA-256. Figure 1 awaits the author’s BioRender revision.
#

# %% [markdown]
# ## Figure 2 — approved PDF assembly
#
# Use the updated `Fig 2 wo panel b.pdf` as the clean base, and
# `Fig 2 partially vectoric.pdf` as the placement reference. Frozen, hashed copies
# of the two manual assets are preserved in
# `figure_rebuild_work/figure2_assembly_2026-10-13/inputs/`.
#
# First run `figure2-asset-*` to render this run’s native PCLO rank table. Assembly
# uses that new `pdf/Figure2b.pdf`, with its source hash recorded.
#
# Place the whole 49 × 45.5 mm panel page at x=367.08 pt, y=10.065 pt, measured
# from the top-left of the base page (x right, y down). Keep the original 510 ×
# 629.04 pt MediaBox/CropBox: 179.916667 × 221.911333 mm. The existing b is retained.
# Docker performs composition and SVG export; no scientific calculations or font
# installation on the host. The final PNG is an uncropped PDF render at 600 DPI.
#

# %%
if 'Figure2' in FIGURES_TO_BUILD:
    from manuscript_figures.docker_export import prepare_figure2, run_figure2_stage, copy_new_figure
    FIGURE2_ASSEMBLY = FINAL_FIGURES_DIR / 'render_work/Figure2'
    figure2_assets = Path(os.environ['MANUSCRIPT_FIGURE2_ASSETS']).resolve()
    prepare_figure2(figure2_assets, FIGURE2_ASSEMBLY, FINAL_FIGURES_DIR / 'pdf/Figure2b.pdf')


# %%
if 'Figure2' in FIGURES_TO_BUILD:
    # A temporary, network-disabled Docker container inserts the vector PDF page.
    # Resource limits: two CPUs, 4 GiB RAM, 180 seconds per stage.
    run_figure2_stage(FIGURE2_ASSEMBLY, 'assemble')


# %%
if 'Figure2' in FIGURES_TO_BUILD:
    # Export SVG with live Unicode text and embedded original glyphs, plus 600 DPI PNG.
    run_figure2_stage(FIGURE2_ASSEMBLY, 'export')


# %%
if 'Figure2' in FIGURES_TO_BUILD:
    import shutil
    from manuscript_figures.docker_export import export_cropped_page
    run_figure2_stage(FIGURE2_ASSEMBLY, 'validate')
    figure2_raw = FIGURE2_ASSEMBLY / 'raw'
    for ext in ['svg', 'pdf']:
        (figure2_raw / ext).mkdir(parents=True, exist_ok=False)
        shutil.copy2(FIGURE2_ASSEMBLY / 'exports' / f'Figure2.{ext}', figure2_raw / ext / f'Figure2.{ext}')
    figure2_finished = export_cropped_page(figure2_raw, FIGURE2_ASSEMBLY / 'finished', 'Figure2',
        [0.0, 5.279999999999999, 510.0, 613.6799780273437], FONT_DIR)
    for ext, folder in [('svg', 'svg'), ('pdf', 'pdf'), ('png', 'previews')]:
        target = FINAL_FIGURES_DIR / folder / f'Figure2.{ext}'
        assert not target.exists(), target
        shutil.copy2(figure2_finished / folder / f'Figure2.{ext}', target)
    display(Image(filename=str(figure2_finished / 'previews/Figure2.png')))


# %%
# Crop only newly rendered, explicitly selected full canvases at physical scale one.
# Figure5 and S7 have dedicated exporters; Figure2 is assembled/cropped below.
from manuscript_figures.docker_export import export_cropped_page
import shutil
# Exact PDF-point rectangles from the approved page exports (not rounded SVG units).
publication_crops = {'Figure3': [0.0, 5.279999999999999, 510.2362060546875, 510.2362060546875],
 'Figure4': [1.92, 11.52, 491.0362060546875, 627.7152880859375],
 'Figure6': [4.32, 1.92, 476.1562060546875, 493.22833251953125],
 'FigureS1': [15.84, 45.12, 510.2362060546875, 487.70270263671875],
 'FigureS2': [0.0, 4.32, 510.2362060546875, 686.1278735351563],
 'FigureS3': [0.0, 41.76, 503.0362060546875, 652.5014379882813],
 'FigureS4': [0.0, 60.48, 484.3162060546875, 352.7092810058594],
 'FigureS5': [0.0, 16.32, 483.3562060546875, 307.49103515625],
 'FigureS6': [0.0, 10.08, 503.0362060546875, 690.9278735351562],
 'FigureS8': [0.0, 27.36, 485.7562060546875, 713.9678735351563],
 'FigureS9': [0.0, 27.36, 485.7562060546875, 713.9678735351563],
 'FigureS10': [0.0, 23.04, 500.1562060546875, 636.6614379882812],
 'FigureS11': [12.959999999999999, 10.08, 442.07620605468753, 689.7146557617187],
 'FigureS12': [1.92, 27.36, 489.1162060546875, 697.6478735351562],
 'FigureS13': [0.0, 44.64, 475.6762060546875, 425.19683837890625],
 'FigureS14': [17.28, 18.72, 464.1562060546875, 308.36031616210937],
 'FigureS15': [1.92, 7.199999999999999, 502.0762060546875, 674.1014379882813],
 'FigureS16': [0.0, 19.2, 468.4762060546875, 642.1417236328125],
 'FigureS17': [0.0, 27.36, 496.7962060546875, 675.0878735351563],
 'FigureS18': [0.0, 15.84, 510.2362060546875, 694.2878735351562],
 'FigureS19': [0.0, 26.4, 453.5433044433594, 428.1033044433594],
 'FigureS20': [21.599999999999998, 36.96, 488.6362060546875, 457.74235595703124]}
for figure_name in list(FIGURES_TO_BUILD) + (['Figure3'] if RENDER_FIGURE3 and 'Figure3' not in FIGURES_TO_BUILD else []):
    if figure_name not in publication_crops:
        continue
    raw_dir = FINAL_FIGURES_DIR / 'render_work/raw'
    finished = FINAL_FIGURES_DIR / 'render_work/cropped' / figure_name
    destinations = [(folder, ext) for folder, ext in [('svg', 'svg'), ('pdf', 'pdf'), ('previews', 'png')]]
    assert all(not (FINAL_FIGURES_DIR / folder / f'{figure_name}.{ext}').exists() for folder, ext in destinations), 'Choose a new FINAL_FIGURES_DIR'
    export_cropped_page(raw_dir, finished, figure_name, publication_crops[figure_name], FONT_DIR)
    if figure_name in {'FigureS11', 'FigureS18'}:
        # These approved pages use Matplotlib's physical crop and Poppler at 600 dpi.
        from matplotlib.transforms import Bbox
        import re
        import subprocess
        figure = figure_s11 if figure_name == 'FigureS11' else figure_s18
        header = re.search(r'<svg\b[^>]*>', (finished / 'svg' / f'{figure_name}.svg').read_text()).group()
        view = [float(value) for value in re.search(r'viewBox="([^"]+)"', header)[1].split()]
        width_mm = float(re.search(r'\bwidth="([\d.]+)mm"', header)[1])
        height_mm = float(re.search(r'\bheight="([\d.]+)mm"', header)[1])
        crop = Bbox.from_bounds(view[0] / 72, figure.get_size_inches()[1] - (view[1] + view[3]) / 72,
                               width_mm / 25.4, height_mm / 25.4)
        figure.savefig(finished / 'pdf' / f'{figure_name}.pdf', bbox_inches=crop, pad_inches=0,
                       metadata={'Creator': 'Combinatorics manuscript figures'})
        subprocess.run(['pdftoppm', '-r', '600', '-singlefile', '-png',
            str(finished / 'pdf' / f'{figure_name}.pdf'), str(finished / 'previews' / figure_name)],
            check=True, timeout=120)
    for folder, ext in destinations:
        shutil.copy2(finished / folder / f'{figure_name}.{ext}', FINAL_FIGURES_DIR / folder / f'{figure_name}.{ext}')

# %% [markdown]
# ## Fisher/BY result handoff
#
# Run the S18 and BY subsections in both scientific notebooks, with a new output directory. Assemble their tables and the S4 workbook using the documented producer functions before rendering. Run the individual Figure 5, S11 and S18 sections above; the selected-page export preserves their final physical crop. The renderer does not run scientific inference or refresh the historical collection.
#

# %%
# Scientific assembly belongs to the producer step, before rendering:
# manuscript_plot_tables.assemble_fisher_tables(producer_dir, new_output_dir)
# manuscript_plot_tables.apply_by_matrices(matrix_table, producer_dir)
# manuscript_workbooks.export_fisher_workbook(template, new_workbook, summary)
# See README.md for native inputs and the explicit order.

# %%
# Render Figure5, FigureS11 and FigureS18 through their individual sections above.
# No nested notebook execution or historical collection refresh is performed.

# %% [markdown]
# ## Entropy comparison — 2026-11-16
#
# The Figure 5c producer compares Observed, Old hypothetical, All reads and
# Originally-supporting reads. Run the comparison setup and producer/export cells
# in combinatorics2. Figures use SVG, PDF and 600 DPI PNG; repeated filenames receive
# a new version. S7 has one dedicated canonical producer in the following section.
#

# %%
COMPARISON_FIGURES_TO_BUILD = []


# %%
# Independent plotting setup using the established Open Sans figure typography.
import sys
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

assert Path(sys.prefix).name == "combinatorics2"
comparison_work = Path("/private6/projects/Combinatorics/Code/Notebooks/entropy_audit_2026-11-12/updated_entropy_and_s7_comparison_2026-11-16")
comparison_notebooks = Path("/private6/projects/Combinatorics/Code/Notebooks")
comparison_final = comparison_notebooks / "final_figures"
comparison_inputs = comparison_notebooks / "figure_rebuild_work/correction_2026-10-12/inputs"
for name in ["OpenSans-Regular.ttf", "OpenSans-Bold.ttf"]:
    font_manager.fontManager.addfont(comparison_inputs.parent / "fonts" / name)
mpl.rcParams.update({"font.family": "Open Sans", "font.size": 6, "axes.labelsize": 7,
    "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": .35, "grid.linewidth": .35, "grid.color": "#dbe4ef",
    "xtick.major.width": .35, "ytick.major.width": .35, "savefig.facecolor": "white"})
comparison_outputs = []


# %%
def export_comparison_figure(figure, stem, output_directory=None):
    manifest_directory = comparison_work if output_directory is None else output_directory
    destinations = [("svg", "svg"), ("pdf", "pdf"), ("previews", "png")]
    output_stem, version = stem, 1
    while True:
        paths = [comparison_final / folder / f"{output_stem}.{extension}" for folder, extension in destinations]
        manifest = manifest_directory / f"{output_stem}.manifest.json"
        if not any(path.exists() for path in [*paths, manifest]):
            break
        version += 1
        output_stem = f"{stem}_v{version}"
    for path in paths:
        figure.savefig(path, dpi=600, facecolor="white")
        comparison_outputs.append({"Figure": output_stem, "Path": str(path), "Format": path.suffix[1:],
            "Bytes": path.stat().st_size, "SHA256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "WidthMM": float(figure.get_figwidth()*25.4), "HeightMM": float(figure.get_figheight()*25.4),
            "PNG_DPI": 600 if path.suffix == ".png" else "vector"})
    with manifest.open("x") as handle:
        json.dump([r for r in comparison_outputs if r["Figure"] == output_stem], handle, indent=2)
    print("Written:", *[str(path) for path in paths], sep="\n")



# %% [markdown]
# ### Figure 5c — Observed, Old hypothetical, All reads and Originally-supporting reads
#

# %%
if "Figure5c_entropy_comparison" in COMPARISON_FIGURES_TO_BUILD:
    support_figure_work = Path("/private6/projects/Combinatorics/Code/Notebooks/entropy_audit_2026-11-12/originally_supporting_reads_2026-11-16")
    comparison_entropy = pd.read_table(support_figure_work / "entropy_population_comparison.tsv", float_precision="round_trip")
    comparison_entropy_previous = pd.read_table(comparison_work / "entropy_updated_per_gene.tsv", float_precision="round_trip")
    comparison_entropy_style = json.loads((comparison_inputs / "squid_original_figure5c_entropy_plot.json").read_text())
    assert len(comparison_entropy) == 21 and not comparison_entropy.duplicated(["Platform", "Gene"]).any()
    for platform, label in [("PacBio", "Long-reads"), ("Illumina", "Short-reads")]:
        assert comparison_entropy.loc[comparison_entropy.Platform.eq(platform), "Gene"].tolist() == comparison_entropy_style["gene_order"][label]
    for label, previous in [("Observed", "CurrentObservedEntropyBits"), ("Old hypothetical", "CurrentOldHypotheticalEntropyBits"),
                            ("All reads", "UpdatedHypotheticalEntropyBits")]:
        np.testing.assert_array_equal(comparison_entropy[label], comparison_entropy_previous[previous])


# %%
if "Figure5c_entropy_comparison" in COMPARISON_FIGURES_TO_BUILD:
    entropy_comparison_figure = plt.figure(figsize=(180/25.4, 88/25.4))
    entropy_comparison_bar_count = 0
    entropy_series = ["Observed", "Old hypothetical", "All reads", "Originally-supporting reads"]
    entropy_upper_limit = max(40, np.ceil(comparison_entropy[entropy_series].max().max()/10)*10)
    for platform, title, left, width in [("PacBio", "Long-reads (PacBio1)", 14, 27),
                                        ("Illumina", "Short-reads (Illumina)", 48, 129)]:
        axis = entropy_comparison_figure.add_axes([left/180, 24/88, width/180, 51/88])
        rows = comparison_entropy.loc[comparison_entropy.Platform.eq(platform)]
        colors = [comparison_entropy_style["gene_colors"][gene] for gene in rows.Gene]
        x = np.arange(len(rows))
        for offset, column, alpha, hatch in [(-.30, "Observed", 1, "////"), (-.10, "Old hypothetical", .35, None),
                (.10, "All reads", 1, None), (.30, "Originally-supporting reads", 1, "...")]:
            outlined = column == "Originally-supporting reads"
            bars = axis.bar(x + offset, rows[column], width=.18, color="white" if outlined else colors,
                alpha=alpha, hatch=hatch, edgecolor=colors if outlined else "white", linewidth=.45, zorder=3)
            np.testing.assert_allclose([bar.get_height() for bar in bars], rows[column], rtol=0, atol=0, equal_nan=True)
            for j, value in enumerate(rows[column]):
                if pd.isna(value): axis.text(x[j]+offset, .8, "N/A", rotation=90, ha="center", fontsize=5)
            entropy_comparison_bar_count += len(bars)
        axis.set(ylim=(0, entropy_upper_limit), xlim=(-.65, len(rows)-.35), title=title)
        axis.set_xticks(x, rows.Gene, rotation=55, ha="right", rotation_mode="anchor")
        axis.set_yticks(np.arange(0, entropy_upper_limit+1, 10))
        axis.spines[:].set_visible(False)
        axis.grid(axis="y", zorder=0)
        axis.tick_params(axis="both", length=0, pad=2, labelleft=platform == "PacBio")
        if platform == "PacBio": axis.set_ylabel("Shannon entropy [bits]")
    entropy_comparison_figure.legend(handles=[
        Patch(facecolor="#555555", edgecolor="white", hatch="////", label="Observed"),
        Patch(facecolor="#bbbbbb", label="Old hypothetical"), Patch(facecolor="#555555", label="All reads"),
        Patch(facecolor="white", edgecolor="#555555", hatch="...", linewidth=.45, label="Originally-supporting reads")],
        loc="lower center", bbox_to_anchor=(.5, .015), ncol=4, frameon=False, columnspacing=1.5)
    assert entropy_comparison_bar_count == 84


# %%
if "Figure5c_entropy_comparison" in COMPARISON_FIGURES_TO_BUILD:
    export_comparison_figure(entropy_comparison_figure, "Figure5c_entropy_comparison", output_directory=support_figure_work)
    plt.close(entropy_comparison_figure)


# %% [markdown]
# ## Figure S7 — corrected AA-to-binary mapping
#
# This is the sole S7 producer. Both panels use the validated 16 November corrected
# scores; panel b connects consecutive cluster sizes with lines and circle markers.
# Run this section in combinatorics2 after selecting FigureS7 and setting a new `FINAL_FIGURES_DIR` in the common path cell. It loads the generated score table and exports at the fixed publication size (600 DPI PNG). It does not run clustering, use a live scientific kernel or replace existing artwork.
#

# %%
# Independent plotting setup for the sole canonical S7 producer.
from pathlib import Path
import os
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import Bbox
BUILD_FIGURE_S7 = 'FigureS7' in FIGURES_TO_BUILD
s7_notebooks = Path.cwd()
from manuscript_paths import figure_source
s7_final_dir = FINAL_FIGURES_DIR
s7_work = s7_final_dir / 'render_work/S7'


# %%
# Preserve the original publication font roles and drawing style.
S7_TEXT_PT = {'panel': 7, 'axis': 7, 'heading': 7, 'legend_heading': 7,
           'tick': 6, 'legend': 6, 'annotation': 6}
for font in ['OpenSans-Regular.ttf', 'OpenSans-Bold.ttf']:
    font_manager.fontManager.addfont(FONT_DIR / font)
assert Path(font_manager.findfont('Open Sans', fallback_to_default=False)).name == 'OpenSans-Regular.ttf'
mpl.rcParams.update({'font.family': 'Open Sans', 'font.size': 6,
    'axes.labelsize': 7, 'axes.titlesize': 7, 'xtick.labelsize': 6,
    'ytick.labelsize': 6, 'legend.fontsize': 6, 'legend.title_fontsize': 7,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.linewidth': 0.35, 'grid.linewidth': 0.35, 'grid.color': '#dbe4ef',
    'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'svg.hashsalt': 'combinatorics-figure3', 'savefig.facecolor': 'white'})

S7_INK = "#2a3f5f"


# %%
# Preserve the original millimetre coordinates and unobtrusive axes.
def s7_axes(figure, x, top, width, height):
    canvas_width, canvas_height = figure.get_size_inches() * 25.4
    return figure.add_axes([x / canvas_width, 1 - (top + height) / canvas_height,
                            width / canvas_width, height / canvas_height])

def s7_text(figure, x, top, label, role, **kwargs):
    canvas_width, canvas_height = figure.get_size_inches() * 25.4
    artist = figure.text(x / canvas_width, 1 - top / canvas_height, label,
        fontsize=S7_TEXT_PT[role], color=S7_INK, **kwargs)
    artist.set_gid(f"role-{role}-{len(figure.texts)}")
    return artist

def s7_style(axis, hide_y=False):
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)
    axis.grid(True, axis="y", zorder=0)
    axis.tick_params(pad=1.5, labelleft=not hide_y)
    axis.tick_params(axis="y", which="both", left=False, right=False)
    axis.tick_params(axis="x", which="both", bottom=False, top=False)
    axis.set_axisbelow(True)
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value / 1000:g}k" if value >= 1000 else f"{value:g}"))
    axis.tick_params(colors=S7_INK)
    for text in [axis.title, axis.xaxis.label, axis.yaxis.label]:
        text.set_color(S7_INK)



# %%
# The scientific producer writes all corrected scores; style is fixed here.
s7_scores = pd.read_table(figure_source('s7_scores'), float_precision='round_trip')
assert s7_scores.columns.tolist() == ['Panel', 'Platform', 'Gene', 'Dissimilarity', 'ClusterSize', 'MeanSilhouetteScore']
assert len(s7_scores) == 1202 and not s7_scores.duplicated(s7_scores.columns[:-1]).any()
assert np.isfinite(s7_scores.MeanSilhouetteScore).all()
s7_colors = {'GRIA2': '#00A08B',
 'PCLO': '#AF0038',
 'ABL': '#0D2A63',
 'ACHA4': '#DA16FF',
 'ANR17': '#FC0080',
 'CA2D3': '#6C4516',
 'CACB2': '#1616A7',
 'CSKI1': '#222A2A',
 'DGLA': '#620042',
 'DOP1': '#A777F1',
 'IQEC1': '#EB663B',
 'K0513': '#750D86',
 'KCNAS': '#FB0D0D',
 'MTUS2': '#1CA71C',
 'RIMS2': '#FB00D1',
 'ROBO2': '#2E91E5',
 'RUSC2': '#B68100',
 'SCN1': '#6C7C32',
 'TRIM2': '#B2828D',
 'TWK7': '#862A16'}
s7_reference_layout = {'s7': {'a': {'groups': ['Long-reads', 'Short-reads'],
       'legend_order': ['ABL',
                        'ACHA4',
                        'ANR17',
                        'CA2D3',
                        'CACB2',
                        'CSKI1',
                        'DGLA',
                        'DOP1',
                        'GRIA2',
                        'IQEC1',
                        'K0513',
                        'KCNAS',
                        'MTUS2',
                        'PCLO',
                        'RIMS2',
                        'ROBO2',
                        'RUSC2',
                        'SCN1',
                        'TRIM2',
                        'TWK7'],
       'limits': {'x': [0, 500.0], 'y': [0, 0.24371506529974304]},
       'reference_asset': 'Fig S7/Mean silhouette score of MiniBatchKMeans vs. '
                          'K size - PacBio1 vs. Illumina - dense.svg',
       'reference_sha256': '9159b391d8183f3324e5efc9c4222d85843a876c720fddde22a6e3af1c67e7d2'},
 'b': {'groups': ['Miyata', 'Grantham 100'],
       'legend_order': ['GRIA2', 'PCLO'],
       'limits': {'x': [0, 42.53156864069147], 'y': [0, 0.10582010582010584]},
       'reference_asset': 'Fig S7/Mean silhouette score of MiniBatchKMeans vs. '
                          'K size - PacBio1 - dissimilar.svg',
       'reference_sha256': '26bbba2fc6a39577165d6e72d15b0fea489748d5601dea9dc03cc3e25b643c92'}}}
s7_tables = {}
for panel in ['a', 'b']:
    data = s7_scores.loc[s7_scores.Panel.eq(panel)].copy()
    data['Platform'] = data.Platform.map({'PacBio': 'Long-reads', 'Illumina': 'Short-reads'})
    data['Dissimilarity'] = data.Dissimilarity.replace({'Grantham100': 'Grantham 100'})
    data['Color'] = data.Gene.map(s7_colors)
    assert data.Platform.notna().all() and data.Color.notna().all()
    assert len(data) == (1050 if panel == 'a' else 152)
    s7_tables[panel] = data


# %%
# Fixed publication crop, independent of the output being replaced.
# Values are copied exactly from the verified pre-edit FigureS7.svg.
s7_crop = [1.9200000525, 8.6400000487, 445.9162182421, 444.9033069513]
s7_width_mm, s7_height_mm = 157.309328247, 156.951999068
s7_bbox = Bbox.from_bounds(s7_crop[0] / 72,
    160 / 25.4 - s7_crop[1] / 72 - s7_height_mm / 25.4,
    s7_width_mm / 25.4, s7_height_mm / 25.4)
if BUILD_FIGURE_S7:
    s7_figure = plt.figure(figsize=(180 / 25.4, 160 / 25.4), dpi=150)
    s7_artist_checks = []


# %%
if BUILD_FIGURE_S7:
    # Panel a: 50 stored points per gene/platform, joined in their original order.
    # Both panels use the same 55 mm columns, 8 mm gutter and legend column.
    S7_COLUMN_X_MM, S7_COLUMN_WIDTH_MM = [18, 81], 55
    S7_X_TITLE_CENTER_MM, S7_LEGEND_X_MM = 77, 144
    s7a_layout = s7_reference_layout['s7']['a']
    for platform, x in zip(s7a_layout['groups'], S7_COLUMN_X_MM):
        axis = s7_axes(s7_figure, x, 12, S7_COLUMN_WIDTH_MM, 60)
        platform_table = s7_tables['a'].loc[s7_tables['a']['Platform'].eq(platform)]
        for gene in s7a_layout['legend_order']:
            curve = platform_table.loc[platform_table['Gene'].eq(gene)]
            if curve.empty:  # A shared gene legend covers both platform cohorts.
                continue
            assert curve['ClusterSize'].tolist() == list(range(10, 501, 10))
            line, = axis.plot(curve['ClusterSize'], curve['MeanSilhouetteScore'],
                              color=s7_colors[gene], linewidth=0.6)
            line.set_gid(f'S7a-{platform}-{gene}')
            np.testing.assert_array_equal(line.get_xydata(), curve[['ClusterSize', 'MeanSilhouetteScore']])
            s7_artist_checks.append({'panel': 'a', 'group': platform, 'gene': gene, 'points': len(curve)})
        axis.set(xlim=s7a_layout['limits']['x'], ylim=s7a_layout['limits']['y'])
        axis.set_xticks(range(0, 501, 100))
        axis.set_yticks(np.arange(0, 0.201, 0.05))
        axis.set_title(platform, pad=4)
        if platform == 'Long-reads':
            axis.set_ylabel('Mean silhouette score of MiniBatchKMeans', labelpad=3)
        s7_style(axis, hide_y=platform == 'Short-reads')
        axis.grid(True, axis='x')
    s7_text(s7_figure, S7_X_TITLE_CENTER_MM, 81, 'Number of clusters (k)', 'axis', ha='center')
    s7_text(s7_figure, 3, 7, 'a', 'panel', weight='bold')

# %%
if BUILD_FIGURE_S7:
    # One alphabetical legend, matching the reference's trace and color ordering.
    s7_text(s7_figure, S7_LEGEND_X_MM, 12, 'Gene', 'legend_heading')
    for row, gene in enumerate(s7a_layout['legend_order']):
        top = 17.2 + row * 2.7
        s7_figure.add_artist(Line2D(np.array([S7_LEGEND_X_MM, S7_LEGEND_X_MM + 3]) / 180,
            [1 - top / 160] * 2, transform=s7_figure.transFigure,
            color=s7_colors[gene], linewidth=0.6))
        s7_text(s7_figure, S7_LEGEND_X_MM + 4.5, top, gene, 'legend', va='center')

# %%
if BUILD_FIGURE_S7:
    # Panel b: connect the stored scores for 3–40 dissimilar clusters.
    s7b_layout = s7_reference_layout['s7']['b']
    for category, x in zip(s7b_layout['groups'], S7_COLUMN_X_MM):
        axis = s7_axes(s7_figure, x, 100, S7_COLUMN_WIDTH_MM, 48)
        for gene in s7b_layout['legend_order']:
            points = s7_tables['b'].loc[s7_tables['b']['Dissimilarity'].eq(category)
                                       & s7_tables['b']['Gene'].eq(gene)]
            assert points['ClusterSize'].tolist() == list(range(3, 41))
            line, = axis.plot(points['ClusterSize'], points['MeanSilhouetteScore'],
                color=s7_colors[gene], linewidth=0.6, linestyle='-',
                marker='o', markersize=1.8, markeredgewidth=0.25,
                clip_on=False)  # Keep complete circles at the k=40 boundary.
            line.set_gid(f'S7b-{category.replace(" ", "")}-{gene}')
            np.testing.assert_array_equal(line.get_xydata(), points[['ClusterSize', 'MeanSilhouetteScore']])
            s7_artist_checks.append({'panel': 'b', 'group': category, 'gene': gene, 'points': len(points)})
        # Remove the reference's automatic X padding so 40 aligns with panel a's 500.
        # PCLO / Grantham 100 / k=5 is 0.109803, above the old 0.105820 limit.
        axis.set(xlim=(0, 40), ylim=(0, 0.12))
        axis.set_xticks(range(0, 41, 10))
        axis.set_yticks(np.arange(0, 0.121, 0.02))
        axis.set_title(category, pad=4)
        if category == 'Miyata':
            axis.set_ylabel('Mean silhouette score of MiniBatchKMeans', labelpad=3)
        s7_style(axis, hide_y=category == 'Grantham 100')
        axis.grid(True, axis='x')
    s7_text(s7_figure, 3, 87, 'b', 'panel', weight='bold')
    s7_text(s7_figure, 18, 90, 'Long-reads', 'heading')
    s7_text(s7_figure, S7_X_TITLE_CENTER_MM, 157, 'Number of dissimilar clusters (k)', 'axis', ha='center')
    s7_text(s7_figure, S7_LEGEND_X_MM, 100, 'Gene', 'legend_heading')
    for row, gene in enumerate(s7b_layout['legend_order']):
        top = 105.2 + row * 2.7  # Same heading offset and item spacing as panel a.
        s7_figure.add_artist(Line2D(np.array([S7_LEGEND_X_MM, S7_LEGEND_X_MM + 3]) / 180,
            [1 - top / 160] * 2, transform=s7_figure.transFigure, linewidth=0.6,
            marker='o', markersize=1.8, markeredgewidth=0.25,
            color=s7_colors[gene], linestyle='-'))
        s7_text(s7_figure, S7_LEGEND_X_MM + 4.5, top, gene, 'legend', va='center')

# %%
# Verify every plotted value, original axes, and text inside the fixed crop.
if BUILD_FIGURE_S7:
    assert sum(item["points"] for item in s7_artist_checks if item["panel"] == "a") == 1050
    assert sum(item["points"] for item in s7_artist_checks if item["panel"] == "b") == 152
    s7_figure.canvas.draw()
    s7_renderer = s7_figure.canvas.get_renderer()
    s7_text_checks = []
    for artist in s7_figure.findobj(Text):
        if not artist.get_visible() or not artist.get_text():
            continue
        bounds = artist.get_window_extent(s7_renderer).extents / s7_figure.dpi
        assert artist.get_fontsize() in [6, 7] and artist.get_fontfamily() == ["Open Sans"]
        assert bounds[0] >= s7_bbox.x0 and bounds[1] >= s7_bbox.y0, artist.get_text()
        assert bounds[2] <= s7_bbox.x1 and bounds[3] <= s7_bbox.y1, artist.get_text()
        s7_text_checks.append({"text": artist.get_text(), "bounds_inches": bounds.tolist(), "size_pt": artist.get_fontsize()})
    for axis in s7_figure.axes:
        assert not any(spine.get_visible() for spine in axis.spines.values())
        for line in axis.lines:
            assert np.asarray(line.get_ydata()).min() >= axis.get_ylim()[0]
            assert np.asarray(line.get_ydata()).max() <= axis.get_ylim()[1]
        for collection in axis.collections:
            assert collection.get_offsets()[:, 1].min() >= axis.get_ylim()[0]
            assert collection.get_offsets()[:, 1].max() <= axis.get_ylim()[1]
    s7b_lines = [line for axis in s7_figure.axes for line in axis.lines
                  if (line.get_gid() or "").startswith("S7b-")]
    assert len(s7b_lines) == 4
    assert all(line.get_linestyle() == "-" and line.get_marker() == "o" for line in s7b_lines)
    s7_checks = {"traces": s7_artist_checks, "text": s7_text_checks,
        "canvas_mm": [180, 160], "export_mm": [s7_width_mm, s7_height_mm],
        "crop_original_viewbox": s7_crop, "all_1202_plotted_values_exact": True,
        "panel_b_lines_and_markers": True,
        "clustering_rerun": False, "axes": [{"xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
            "position": list(ax.get_position().bounds)} for ax in s7_figure.axes]}


# %%
# Export into the chosen destination; existing artwork is never replaced.
if BUILD_FIGURE_S7:
    s7_stem = "FigureS7"
    s7_outputs = []
    s7_work.mkdir(parents=True, exist_ok=True)
    for folder, extension in [('svg', 'svg'), ('pdf', 'pdf'), ('previews', 'png')]:
        (s7_final_dir / folder).mkdir(parents=True, exist_ok=True)
        assert not (s7_final_dir / folder / f'{s7_stem}.{extension}').exists(), 'Choose a new output directory'
    for folder, extension in [("svg", "svg"), ("pdf", "pdf"), ("previews", "png")]:
        path = s7_final_dir / folder / f"{s7_stem}.{extension}"
        s7_figure.savefig(path, format=extension, dpi=600,
            bbox_inches=s7_bbox, pad_inches=0)
        s7_outputs.append({"Figure": s7_stem, "Path": str(path), "Format": extension,
            "SHA256": hashlib.sha256(path.read_bytes()).hexdigest(), "Bytes": path.stat().st_size,
            "WidthMM": s7_width_mm, "HeightMM": s7_height_mm,
            "PNG_DPI": 600 if extension == "png" else "vector"})
    (s7_work / "FigureS7.manifest.json").write_text(
        json.dumps({"outputs": s7_outputs, "checks": s7_checks}, indent=2) + "\n")
    print("Written:\n" + "\n".join(row["Path"] for row in s7_outputs))


# %%
# Per-figure metadata stays with this export; no historical collection is refreshed.
if BUILD_FIGURE_S7:
    from PIL import Image as PILImage
    with PILImage.open(s7_final_dir / 'previews/FigureS7.png') as image:
        assert all(abs(dpi - 600) < .01 for dpi in image.info['dpi'])
        s7_dimensions = {'Figure': 'FigureS7', 'WidthMM': s7_width_mm, 'HeightMM': s7_height_mm,
                         'PixelWidth': image.width, 'PixelHeight': image.height, 'DPI': 600}
    pd.DataFrame([s7_dimensions]).to_csv(s7_work / 'dimensions.tsv', sep='\t', index=False)

