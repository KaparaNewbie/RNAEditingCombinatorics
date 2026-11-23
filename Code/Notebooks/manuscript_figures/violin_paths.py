"""Read the original absolute SVG Bézier contours without re-estimating a KDE."""
import re
from matplotlib.path import Path


def original_violin_path(data):
    tokens = re.findall(r'[A-Za-z]|[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?', data)
    vertices, codes = [], []
    index = 0
    commands = {'M': (1, Path.MOVETO), 'L': (1, Path.LINETO),
                'Q': (2, Path.CURVE3), 'C': (3, Path.CURVE4)}
    while index < len(tokens):
        command = tokens[index]
        index += 1
        if command == 'Z':
            vertices.append(vertices[0])
            codes.append(Path.CLOSEPOLY)
            continue
        assert command in commands, command
        count, code = commands[command]
        values = list(map(float, tokens[index:index+count*2]))
        assert len(values) == count*2
        vertices.extend(zip(values[::2], values[1::2]))
        codes.extend([code]*count)
        index += count*2
    return Path(vertices, codes)
