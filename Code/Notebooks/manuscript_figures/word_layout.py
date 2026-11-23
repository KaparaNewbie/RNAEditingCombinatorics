"""Presentation-only export and measured legend clearance for 14 October."""
from pathlib import Path
import hashlib
import numpy as np
from . import batch_helpers

ROOT = batch_helpers.FINAL / 'render_work'


def place_clear_legends(figure, genes):
    """Choose the first empty rendered legend rectangle in a fixed candidate order."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    results = []
    for axis, gene in zip(figure.axes, genes, strict=True):
        legend = axis.get_legend()
        original = legend.get_window_extent(renderer).frozen()
        points, radii = [], []
        for collection in axis.collections:
            xy = collection.get_offset_transform().transform(collection.get_offsets())
            points.extend(xy)
            radii.extend(np.broadcast_to(np.sqrt(collection.get_sizes())/2*figure.dpi/72,
                                        (len(xy),)))
        points, radii = np.asarray(points), np.asarray(radii)
        clearance = radii + .4*figure.dpi/25.4
        bounds = axis.get_window_extent(renderer)
        inset = .5*figure.dpi/25.4
        width, height = original.width, original.height
        # Prefer common corners, then a reproducible top-to-bottom/right-to-left grid.
        fractions = [(1,1),(0,1),(1,0),(0,0)]
        fractions += [(x,y) for y in np.linspace(1,0,21) for x in np.linspace(1,0,31)]
        legend.set_loc('lower left')
        legend.borderaxespad = 0
        chosen = None
        for fx,fy in fractions:
            x = bounds.x0+inset+fx*(bounds.width-width-2*inset)
            y = bounds.y0+inset+fy*(bounds.height-height-2*inset)
            anchor = axis.transAxes.inverted().transform((x,y))
            legend.set_bbox_to_anchor(anchor, transform=axis.transAxes)
            box = legend.get_window_extent(renderer).padded(.3*figure.dpi/72)
            distances = np.hypot(np.maximum.reduce([box.x0-points[:,0],points[:,0]-box.x1,np.zeros(len(points))]),
                                 np.maximum.reduce([box.y0-points[:,1],points[:,1]-box.y1,np.zeros(len(points))]))
            if np.all(distances > clearance):
                chosen = {'gene':gene,'anchor_axes':anchor.tolist(),
                          'legend_bounds_mm':(box.extents*25.4/figure.dpi).tolist(),
                          'minimum_marker_clearance_mm':float(np.min(distances-radii)*25.4/figure.dpi),
                          'points_checked':len(points),'occluded_markers':0}
                break
        assert chosen is not None, (gene,'No empty in-panel legend position; allocate external space.')
        old_distance = np.hypot(np.maximum.reduce([original.x0-points[:,0],points[:,0]-original.x1,np.zeros(len(points))]),
                                np.maximum.reduce([original.y0-points[:,1],points[:,1]-original.y1,np.zeros(len(points))]))
        chosen['original_occluded_markers'] = int((old_distance <= radii).sum())
        results.append(chosen)
    figure.canvas.draw()
    return results


def export_layout_figure(figure, name, data_checks):
    """Reuse established validation/export, retaining raw native-size vector pages."""
    previous_final, previous_evidence = batch_helpers.FINAL, batch_helpers.FINALIZATION
    for folder in ['pdf','svg','previews']:
        (ROOT/'raw'/folder).mkdir(parents=True,exist_ok=True)
    # Each scientific trace is fingerprinted before any exterior cropping.
    arrays = []
    for index, axis in enumerate(figure.axes):
        for kind, artists in [('line',axis.lines),('collection',axis.collections)]:
            for artist in artists:
                values = artist.get_xydata() if kind=='line' else artist.get_offsets()
                array = np.asarray(values)
                arrays.append({'axis':index,'kind':kind,'shape':list(array.shape),
                               'sha256':hashlib.sha256(array.tobytes()).hexdigest()})
    data_checks['artist_array_fingerprints'] = arrays
    batch_helpers.FINAL, batch_helpers.FINALIZATION = ROOT/'raw', ROOT
    try:
        batch_helpers.export_figure(figure,name,data_checks)
    finally:
        batch_helpers.FINAL, batch_helpers.FINALIZATION = previous_final, previous_evidence


def balance_middle_panel(figure, upper, middle, lower):
    """Translate the complete middle group until both visible inter-panel gaps match."""
    from matplotlib.axes import Axes
    from matplotlib.patches import Rectangle
    from matplotlib.text import Text
    from matplotlib.transforms import Bbox

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    def bounds(group):
        boxes = [artist.get_tightbbox(renderer) for artist in group]
        assert boxes and all(box is not None and np.isfinite(box.extents).all() for box in boxes)
        return Bbox.union(boxes)

    a, b, c = [bounds(group) for group in [upper, middle, lower]]
    before = np.array([a.y0-b.y1, b.y0-c.y1]) * 25.4 / figure.dpi
    shift_pixels = ((a.y0-b.y1) - (b.y0-c.y1)) / 2
    shift_fraction = shift_pixels / figure.bbox.height
    for artist in middle:
        if isinstance(artist, Axes):
            x, y, width, height = artist.get_position().bounds
            artist.set_position([x, y+shift_fraction, width, height])
        elif isinstance(artist, Text):
            x, y = artist.get_position()
            artist.set_position((x, y+shift_fraction))
        elif isinstance(artist, Rectangle):
            artist.set_y(artist.get_y()+shift_fraction)
        else:
            raise TypeError(f'Unsupported middle-panel artist: {type(artist).__name__}')
    figure.canvas.draw()
    a, b, c = [bounds(group) for group in [upper, middle, lower]]
    after = np.array([a.y0-b.y1, b.y0-c.y1]) * 25.4 / figure.dpi
    np.testing.assert_allclose(after[0], after[1], atol=1e-8)
    assert after.min() >= 1.5
    return {'before_gaps_mm':before.tolist(), 'after_gaps_mm':after.tolist(),
            'middle_shift_down_mm':float(-shift_pixels*25.4/figure.dpi),
            'includes_axes_titles_ticks_and_legends':True}
