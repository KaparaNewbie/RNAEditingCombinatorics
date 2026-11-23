"""Export a Plotly SVG at declared physical size, preserving vector text in PDF."""
from pathlib import Path
import os
import subprocess
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape
import plotly.io as pio
from .vector_pdf import svg_to_pdf

NOTEBOOK_DIR = Path(__file__).resolve().parents[1]
FINAL = Path(os.environ.get('MANUSCRIPT_FINAL_FIGURES', NOTEBOOK_DIR / 'final_figures/new_run')).resolve()


def export_plotly(figure, name, width_mm, height_mm, *, output_dir=None, font_dir=None, browser=None):
    # Four SVG coordinates per mm; notebook font sizes use the same conversion.
    assert figure.layout.width == width_mm * 4
    assert figure.layout.height == height_mm * 4
    final = Path(output_dir) if output_dir is not None else FINAL
    evidence = final / 'render_work/plotly'
    evidence.mkdir(parents=True, exist_ok=True)
    for folder in ['svg', 'pdf', 'previews']:
        (final / folder).mkdir(parents=True, exist_ok=True)
    chrome = Path(browser or os.environ.get('BROWSER_PATH', Path.home() / '.local/share/choreographer/deps/chrome-linux64/chrome'))
    assert chrome.is_file(), 'Existing Chrome executable is required by installed Kaleido'
    os.environ['BROWSER_PATH'] = str(chrome)
    fontconfig = evidence / 'fontconfig.xml'
    fonts = Path(font_dir) if font_dir is not None else Path(os.environ['MANUSCRIPT_FONTS'])
    fontconfig.write_text('<?xml version="1.0"?><fontconfig><dir>'+escape(str(fonts))+
                         '</dir><cachedir>'+escape(str(evidence/'fontconfig-cache'))+'</cachedir></fontconfig>')
    environment = dict(os.environ, FONTCONFIG_FILE=str(fontconfig))
    os.environ['FONTCONFIG_FILE'] = str(fontconfig)
    pio.defaults.mathjax = None
    svg = pio.to_image(figure, format='svg')
    root = ET.fromstring(svg)
    assert root.get('viewBox') == f'0 0 {int(width_mm*4)} {int(height_mm*4)}'
    # Set physical units without altering any path, text, or plot coordinates.
    # This shared export step is repeated by the notebook, not a manual graphic edit.
    root.set('width', f'{width_mm}mm')
    root.set('height', f'{height_mm}mm')
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    svg_path = final / 'svg' / f'{name}.svg'
    ET.ElementTree(root).write(svg_path, encoding='utf-8', xml_declaration=True)
    svg_to_pdf(svg_path, final/'pdf'/f'{name}.pdf', width_mm, height_mm)
    subprocess.run(['rsvg-convert','-d','300','-p','300','-o',str(final/'previews'/f'{name}.png'),str(svg_path)],
                   env=environment, check=True, timeout=60)
    figure.write_json(evidence/f'{name}_plotly.json')
    print(f'{name}: Plotly vector SVG/PDF and 300 dpi PNG at {width_mm} × {height_mm} mm.')
