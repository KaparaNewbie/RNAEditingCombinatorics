"""Crop one vector page at scale one; run in the documented PyMuPDF image.

Input PDF/SVG and fonts are read-only. The destination must be a new directory.
The rectangle is the approved physical viewport in PDF points, top-left origin.
"""
from pathlib import Path
import argparse
import base64
import re
import shutil
import subprocess
import sys
import pymupdf as fitz
from PIL import Image

MM = 25.4 / 72

def crop_svg(source, destination, crop, original_page, font_dir):
    """Change the viewport only; retain every vector/text element byte for byte."""
    svg = source.read_text()
    match = re.search(r'<svg\b[^>]*>',svg)
    header = match.group()
    view = [float(v) for v in re.search(r'viewBox="([^"]+)"',header)[1].split()]
    scale_x,scale_y = view[2]/original_page.width,view[3]/original_page.height
    assert abs(scale_x/scale_y-1)<1e-6
    new_view = [view[0]+crop.x0*scale_x,view[1]+crop.y0*scale_y,
                crop.width*scale_x,crop.height*scale_y]
    new = re.sub(r'\bwidth="[^"]+"',f'width="{crop.width*MM:.9f}mm"',header,count=1)
    new = re.sub(r'\bheight="[^"]+"',f'height="{crop.height*MM:.9f}mm"',new,count=1)
    new = re.sub(r'\bviewBox="[^"]+"','viewBox="'+' '.join(f'{v:.10f}' for v in new_view)+'"',new,count=1)
    result = svg[:match.start()]+new+svg[match.end():]
    # Embed the existing local fonts so review does not depend on installed Open Sans.
    if ('Open Sans' in result or source.stem=='Figure2b') and '@font-face' not in result:
        styles=[]
        family='Arial' if source.stem=='Figure2b' else 'Open Sans'
        font_files=[('Arial.TTF',400)] if family=='Arial' else [('OpenSans-Regular.ttf',400),('OpenSans-Bold.ttf',700)]
        for file,weight in font_files:
            encoded=base64.b64encode((font_dir/file).read_bytes()).decode()
            styles.append('@font-face{font-family:"'+family+'";font-weight:'+str(weight)+';src:url(data:font/ttf;base64,'+encoded+') format("truetype");}')
        position=result.index('>',result.index('<svg'))+1
        result=result[:position]+'<defs><style>'+''.join(styles)+'</style></defs>'+result[position:]
    destination.write_text(result)
    counts={tag:len(re.findall(r'<(?:\w+:)?'+tag+r'\b',result)) for tag in ['text','path','image']}
    for tag,count in counts.items():
        assert count==len(re.findall(r'<(?:\w+:)?'+tag+r'\b',svg)),(tag,count)
    assert counts['text']>0
    return counts


def crop_page(source_dir, output_dir, name, rectangle, font_dir):
    output_dir.mkdir(parents=True, exist_ok=False)
    for folder in ['pdf', 'svg', 'previews']:
        (output_dir / folder).mkdir()
    with fitz.open(source_dir / 'pdf' / f'{name}.pdf') as original:
        assert len(original) == 1 and original[0].rotation == 0
        page = original[0]
        crop = fitz.Rect(rectangle)
        assert page.rect.contains(crop)
        destination = output_dir / 'pdf' / f'{name}.pdf'
        document = fitz.open()
        new = document.new_page(width=crop.width, height=crop.height)
        xref = new.show_pdf_page(new.rect, original, 0, clip=crop, keep_proportion=True)
        assert document.xref_stream(xref) == page.read_contents()
        document.set_metadata(original.metadata)
        document.save(destination, garbage=3, deflate=True, no_new_id=True)
        document.close()
        crop_svg(source_dir / 'svg' / f'{name}.svg', output_dir / 'svg' / f'{name}.svg', crop, page.rect, font_dir)
    # A fresh process prevents MuPDF font-state leakage between source and final PDF.
    subprocess.run([sys.executable, __file__, 'render', str(output_dir), name], check=True, timeout=180)


def render_page(output_dir, name):
    with fitz.open(output_dir / 'pdf' / f'{name}.pdf') as document:
        page = document[0]
        pixmap = page.get_pixmap(dpi=600, alpha=False, colorspace=fitz.csRGB)
        image = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
        expected = (round(page.rect.width * 600 / 72), round(page.rect.height * 600 / 72))
        assert abs(image.width - expected[0]) <= 1 and abs(image.height - expected[1]) <= 1
        canvas = Image.new('RGB', expected, 'white')
        canvas.paste(image, (0, 0))  # At most one blank outer pixel; no resampling.
        canvas.save(output_dir / 'previews' / f'{name}.png', dpi=(600, 600), compress_level=6)


if __name__ == '__main__':
    if sys.argv[1] == 'render':
        render_page(Path(sys.argv[2]), sys.argv[3])
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--source', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument('--name', required=True)
        parser.add_argument('--crop', nargs=4, type=float, required=True)
        parser.add_argument('--fonts', type=Path, required=True)
        args = parser.parse_args()
        crop_page(args.source, args.output, args.name, args.crop, args.fonts)
