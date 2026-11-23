"""Shared physical placement, exact table loading and vector export checks."""
import hashlib
import json
import os
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.text import Text

NOTEBOOK_DIR = Path(__file__).resolve().parents[1]
FINAL = Path(os.environ.get('MANUSCRIPT_FINAL_FIGURES', NOTEBOOK_DIR / 'final_figures/new_run')).resolve()
FINALIZATION = FINAL / 'render_work'
INK = '#2a3f5f'
SIZES = {'panel': 7, 'axis': 7, 'heading': 7, 'legend_heading': 7,
         'legend': 6, 'tick': 6, 'annotation': 6, 'subtitle': 6}


def load_table(name, settings=None, input_dir=None):
    from manuscript_paths import figure_source
    path = Path(input_dir) / f'{name}.tsv' if input_dir is not None else figure_source(name)
    if settings is not None:
        return pd.read_table(path, float_precision='round_trip'), settings
    metadata = json.loads(path.with_suffix('.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == metadata['sha256']
    table = pd.read_table(path, float_precision='round_trip')
    assert len(table) == metadata['rows'] and list(table) == metadata['columns']
    return table, metadata


def axes_mm(figure, x, top, width, height, grid='y', hide_y=False):
    canvas_width, canvas_height = figure.get_size_inches() * 25.4
    axis = figure.add_axes([x / canvas_width, 1 - (top + height) / canvas_height,
                           width / canvas_width, height / canvas_height])
    axis.spines[:].set_visible(False)
    axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False,
                     right=False, labelleft=not hide_y, pad=1.5, colors=INK)
    axis.set_axisbelow(True)
    if grid:
        axis.grid(True, axis=grid, color='#dbe4ef', linewidth=0.35)
    for text in [axis.title, axis.xaxis.label, axis.yaxis.label]:
        text.set_color(INK)
    return axis


def text_mm(figure, x, top, label, role='heading', **kwargs):
    width, height = figure.get_size_inches() * 25.4
    if role == 'panel':
        kwargs.setdefault('weight', 'bold')
    artist = figure.text(x / width, 1 - top / height, label,
                        color=INK, fontsize=SIZES[role], **kwargs)
    artist.set_gid(f'role-{role}-{len(figure.texts)}')
    return artist


def check_xy(artist, x, y, kind='line'):
    expected = np.column_stack([x, y])
    actual = artist.get_offsets() if kind == 'scatter' else artist.get_xydata()
    np.testing.assert_array_equal(np.asarray(actual), expected)


def inline_svg_markers(path):
    """Preserve native paths; share their geometry through internal XML entities.

    Unlike <use>, these expand to ordinary paths in the available librsvg.
    Each marker keeps its own transform/style, including alpha compositing.
    No external DTD or external entity is used.
    """
    svg_ns = 'http://www.w3.org/2000/svg'
    href = '{http://www.w3.org/1999/xlink}href'
    ET.register_namespace('', svg_ns)
    ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
    tree = ET.parse(path)
    definitions = {e.get('id'): e for e in tree.iter(f'{{{svg_ns}}}path') if e.get('id')}
    geometries = {}
    count = 0
    for element in tree.iter(f'{{{svg_ns}}}use'):
        original = dict(element.attrib)
        marker = definitions[original[href][1:]]
        geometry = marker.attrib['d']
        key = geometries.setdefault(geometry, f'marker{len(geometries)}')
        element.tag = marker.tag
        element.attrib.clear()
        element.set('d', f'__{key}__')
        translation = f'translate({original.get("x", "0")} {original.get("y", "0")})'
        element.set('transform', (original.get('transform', '') + ' ' + translation).strip())
        element.set('style', marker.get('style', '') + ';' + original.get('style', ''))
        count += 1
    temporary = path.with_suffix('.serialization.tmp')
    tree.write(temporary, encoding='utf-8', xml_declaration=False)
    with path.open('w') as output, temporary.open() as source:
        output.write('<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE svg [\n')
        for geometry, key in geometries.items():
            assert '"' not in geometry and '&' not in geometry and '<' not in geometry
            output.write(f'<!ENTITY {key} "{geometry}">\n')
        output.write(']>\n')
        for line in source:
            for key in geometries.values():
                line = line.replace(f'__{key}__', f'&{key};')
            output.write(line)
    temporary.unlink()
    return count


def export_figure(figure, name, data_checks):
    """Keep the full physical canvas; audit original artists before vector export."""
    for folder in ['svg', 'pdf', 'previews']:
        (FINAL / folder).mkdir(parents=True, exist_ok=True)
    (FINALIZATION / 'evidence').mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    width, height = figure.get_size_inches() * 25.4
    texts, axes = [], []
    for artist in figure.findobj(Text):
        if not artist.get_visible() or not artist.get_text():
            continue
        compact_s18 = name == 'FigureS18' and artist.get_gid() == 's18-author-approved-compact-text'
        assert (0 < artist.get_fontsize() <= 7) if compact_s18 else artist.get_fontsize() in [6, 7], (artist.get_text(), artist.get_fontsize())
        assert artist.get_fontfamily() == ['Open Sans']
        bounds = artist.get_window_extent(renderer).extents / figure.dpi * 25.4
        texts.append({'text': artist.get_text(), 'size_pt': artist.get_fontsize(),
                      'bounds_mm': bounds.tolist(), 'rotation': artist.get_rotation(),
                      'gid': artist.get_gid(), 'weight': artist.get_fontweight()})
    for axis in figure.axes:
        assert not any(s.get_visible() for s in axis.spines.values())
        ticks = (axis.xaxis.get_major_ticks() + axis.xaxis.get_minor_ticks()
                 + axis.yaxis.get_major_ticks() + axis.yaxis.get_minor_ticks())
        assert not any(t.tick1line.get_visible() or t.tick2line.get_visible() for t in ticks)
        assert not any(a.get_rasterized() for a in axis.get_children())
        axes.append({'position_fraction': list(axis.get_position().bounds),
                     'xlim': list(axis.get_xlim()), 'ylim': list(axis.get_ylim()),
                     'xscale': axis.get_xscale(), 'yscale': axis.get_yscale(),
                     'visible_frames': False, 'visible_tick_marks': False})
    clipped = [t for t in texts if t['bounds_mm'][0] < 0 or t['bounds_mm'][1] < 0
               or t['bounds_mm'][2] > width or t['bounds_mm'][3] > height]
    assert not clipped, [(t['text'], t['bounds_mm']) for t in clipped]
    for extension in ['svg', 'pdf']:
        figure.savefig(FINAL / extension / f'{name}.{extension}', metadata={'Creator': 'Combinatorics manuscript figures'})
    if data_checks.get('points', 0) > 500000:
        # Serialization only: identical marker geometry, position, color and opacity.
        data_checks['native_svg_marker_paths_inlined'] = inline_svg_markers(FINAL / 'svg' / f'{name}.svg')
    report = {'canvas_mm': [width, height], 'texts': texts, 'axes': axes,
              'outside_canvas': clipped, 'data_checks': data_checks,
              'render_export_seconds': time.perf_counter() - started}
    (FINALIZATION / 'evidence' / f'{name}_artists.json').write_text(json.dumps(report, indent=2))
    print(f'{name}: full-canvas vector SVG/PDF at {width:g} x {height:g} mm', flush=True)
    figure.savefig(FINAL / 'previews' / f'{name}.png', dpi=300)
    figure.savefig(FINALIZATION / 'evidence' / f'{name}_print_scale.png', dpi=96)
    plt.close(figure)


def continuous_colorbar(axis, scale, vmin, vmax, ticks):
    """Standard continuous-normalization colorbar, with seam-free vector solids."""
    cmap = mpl.colors.LinearSegmentedColormap.from_list('reference', scale, N=1024)
    scalar = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin, vmax), cmap=cmap)
    colorbar = axis.figure.colorbar(scalar, cax=axis, ticks=ticks)
    colorbar.formatter = mpl.ticker.FuncFormatter(lambda value, _: f'{value:g}'.replace('-', '−'))
    colorbar.update_ticks()
    colorbar.solids.set_rasterized(False)
    colorbar.solids.set_edgecolor('face')
    colorbar.solids.set_linewidth(.1)
    colorbar.solids.set_antialiased(False)
    axis.spines[:].set_visible(False)
    axis.tick_params(which='both', left=False, right=False, top=False, bottom=False)


def check_layout_spacing(figure):
    """Measure external heading/tick clearance and adjacent facet row clearance."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    mm = 25.4 / figure.dpi
    headings = [t for t in figure.texts + [a.title for a in figure.axes]
                if t.get_visible() and t.get_text() and t.get_rotation() == 0]
    ticks = [t for a in figure.axes for t in a.get_xticklabels()
             if t.get_visible() and t.get_text()]
    heading_boxes = [(t, t.get_window_extent(renderer)) for t in headings]
    tick_boxes = [(t, t.get_window_extent(renderer)) for t in ticks]
    for title, box in heading_boxes:
        for tick, tick_box in tick_boxes:
            assert not box.overlaps(tick_box), ('heading overlaps X tick', title.get_text(), tick.get_text())
    row_gaps = []
    for lower in figure.axes:
        lower_box = lower.get_window_extent(renderer)
        titles = [(t, b) for t, b in heading_boxes
                  if lower_box.x0 < (b.x0+b.x1)/2 < lower_box.x1
                  and -.1 <= (b.y0-lower_box.y1)*mm <= 6]
        if not titles:
            continue
        candidates = [upper for upper in figure.axes
                      if upper.get_window_extent(renderer).y0 > lower_box.y1
                      and upper.get_window_extent(renderer).x0 < lower_box.x1
                      and upper.get_window_extent(renderer).x1 > lower_box.x0]
        if not candidates:
            continue
        upper = min(candidates, key=lambda a: a.get_window_extent(renderer).y0)
        upper_ticks = [t.get_window_extent(renderer) for t in upper.get_xticklabels()
                       if t.get_visible() and t.get_text()]
        boundary = min([b.y0 for b in upper_ticks] or [upper.get_window_extent(renderer).y0])
        for title, box in titles:
            gap = (boundary-box.y1)*mm
            assert gap >= 1.5, ('facet row clearance < 1.5 mm', title.get_text(), gap)
            row_gaps.append({'lower_title': title.get_text(), 'clearance_mm': gap})
    return {'heading_tick_pairs_checked': len(heading_boxes)*len(tick_boxes),
            'heading_tick_overlaps': 0, 'adjacent_row_gaps': row_gaps}


def close_facet_title(axis, label, grid_top=None):
    """Put the text bottom 0.6 mm above the visible plotting area's top."""
    height_mm = axis.get_position().height * axis.figure.get_figheight() * 25.4
    plot_top = 1 if grid_top is None else axis.transAxes.inverted().transform(
        axis.transData.transform((axis.get_xlim()[0], grid_top)))[1]
    title = axis.set_title(label, fontsize=7, y=plot_top + .6/height_mm,
                          pad=0, va='bottom')
    return axis, title, plot_top


def check_facet_title_attachment(facets):
    """Require titles to sit closer to their own plot than the preceding row."""
    figure = facets[0][0].figure
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    mm = 25.4/figure.dpi
    report = []
    for axis, title, plot_top in facets:
        box = title.get_window_extent(renderer)
        top_px = axis.transAxes.transform((.5, plot_top))[1]
        own_gap = (box.y0-top_px)*mm
        assert .55 <= own_gap <= .65, (title.get_text(), 'own plot gap', own_gap)
        bounds = axis.get_window_extent(renderer)
        preceding = [a for a, _, _ in facets if a.get_window_extent(renderer).y0 > bounds.y1
                     and a.get_window_extent(renderer).x0 < bounds.x1
                     and a.get_window_extent(renderer).x1 > bounds.x0]
        gap = None
        if preceding:
            above = min(preceding, key=lambda a: a.get_window_extent(renderer).y0)
            labels = [t.get_window_extent(renderer).y0 for t in above.get_xticklabels()
                      if t.get_visible() and t.get_text()]
            gap = (min(labels)-box.y1)*mm
            assert gap >= 4 and gap >= 4*own_gap, (title.get_text(), 'preceding row gap', gap)
        report.append({'title':title.get_text(), 'own_plot_gap_mm':own_gap,
                       'preceding_row_clearance_mm':gap})
    return report


def annotate_rank_markers(axis, details, line, platform, gene):
    """Jointly place nearby labels, retaining marker order and clear association."""
    from matplotlib.path import Path as ArtistPath
    from itertools import product
    renderer = axis.figure.canvas.get_renderer()
    curve = ArtistPath(line.get_transform().transform(line.get_xydata()))
    bounds = axis.get_window_extent(renderer)
    candidates = [(dx, dy, 'left', 'center') for dx in [4, 6, 8]
                  for dy in [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5,
                             -6, 6, -7, 7, -8, 8, -9, 9, -10, 10]
                  if np.hypot(dx, dy)*25.4/72 <= 4]
    candidates += [(dx, dy, 'center', 'top') for dy in [-4, -6, -8]
                   for dx in [0, 4, 8] if np.hypot(dx, dy)*25.4/72 <= 4]
    options = []
    points = list(details.itertuples())
    marker_positions = axis.transData.transform([(p.RankLabel,p.CumulativeExpression) for p in points])
    mm = 25.4/axis.figure.dpi
    for point_index, point in enumerate(points):
        valid = []
        for dx, dy, horizontal, vertical in candidates:
            label = axis.annotate(point.Label, (point.RankLabel, point.CumulativeExpression),
                xytext=(dx, dy), textcoords='offset points', fontsize=6,
                ha=horizontal, va=vertical, color=INK)
            box = label.get_window_extent(renderer).padded(.7 * axis.figure.dpi / 72)
            inside = bounds.contains(box.x0, box.y0) and bounds.contains(box.x1, box.y1)
            distances = [np.hypot(max(box.x0-x, 0, x-box.x1),
                                   max(box.y0-y, 0, y-box.y1))*mm for x,y in marker_positions]
            association_clear = distances[point_index] <= min(distances)+.05
            if inside and association_clear and not curve.intersects_bbox(box, filled=False):
                valid.append((dx, dy, horizontal, vertical, box, distances[point_index]))
            label.remove()
        assert valid, (platform, gene, point.RankLabel, 'no curve-free label candidate')
        options.append(valid)
    # Small shifts of several labels are preferable to one distant label.
    ordered_trials = (trial for trial in product(*options)
        if not any(trial[i][4].overlaps(trial[j][4]) for i,j in [(0,1),(0,2),(1,2)])
        and all((trial[i][4].y0+trial[i][4].y1) < (trial[i+1][4].y0+trial[i+1][4].y1)
                for i in [0,1]))
    placements = min(ordered_trials, key=lambda trial:
        (max(np.hypot(t[0],t[1]) for t in trial), sum(np.hypot(t[0],t[1]) for t in trial)), default=None)
    assert placements is not None, (platform, gene, 'rank labels overlap')
    report = []
    for point, (dx,dy,horizontal,vertical,box,distance) in zip(points, placements):
        axis.annotate(point.Label, (point.RankLabel, point.CumulativeExpression),
            xytext=(dx,dy), textcoords='offset points', fontsize=6,
            ha=horizontal, va=vertical, color=INK)
        report.append({'platform':platform,'gene':gene,'rank':int(point.RankLabel),
            'value':float(point.CumulativeExpression),'offset_pt':[dx,dy],
            'text_bounds_px':list(box.extents),'curve_or_text_collision':False})
        report[-1].update({'marker_to_label_edge_mm':distance,
            'offset_distance_mm':float(np.hypot(dx,dy)*25.4/72),
            'own_marker_is_nearest':True,'vertical_rank_order_preserved':True})
    return report
