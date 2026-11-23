"""Convert PDF graphics with MuPDF; retain Unicode text and the embedded glyphs.

BioRender's PDF subsets omit cmap/name tables. MuPDF's SVG writer conflates the
regular and bold subsets, so text must be paired with the PDF text trace. Restore
only font metadata/Unicode maps; never change glyph outlines, advances or positions.
"""
from pathlib import Path
from io import BytesIO
import base64
import json
import re
import xml.etree.ElementTree as ET
import numpy as np
from fontTools.ttLib import TTFont
from fontTools.fontBuilder import FontBuilder

SVG = 'http://www.w3.org/2000/svg'


def export_live_svg(document, destination, fonts_dir, evidence_dir):
    page = document[0]
    tree = ET.fromstring(page.get_svg_image(text_as_path=False))
    traces = page.get_texttrace()
    texts = tree.findall(f'.//{{{SVG}}}text')
    assert len(texts) == len(traces)
    font_maps = {}
    for span in traces:
        assert span['type'] == 0 and span['opacity'] == 1 and span['wmode'] == 0
        mapping = font_maps.setdefault(span['font'], {})
        for unicode, glyph, origin, box in span['chars']:
            assert unicode >= 32 and (unicode not in mapping or mapping[unicode] == glyph)
            mapping[unicode] = glyph
    fonts = {}
    for record in page.get_fonts(full=True):
        basename, extension, kind, content = document.extract_font(record[0])
        name = basename.split('+')[-1]
        if name in fonts:
            assert fonts[name] == content
        fonts[name] = content
    assert set(font_maps) == set(fonts)
    stylesheet, aliases, records = [], {}, []
    for index, (name, mapping) in enumerate(font_maps.items()):
        font = TTFont(BytesIO(fonts[name]))
        font.recalcTimestamp = False
        glyph_order = font.getGlyphOrder()
        # FontTools may compress repeated hmtx advances; compare metric values,
        # not that table's byte encoding or its resulting numberOfHMetrics.
        outlines_before = font.getTableData('glyf')
        advances_before = dict(font['hmtx'].metrics)
        units_before = font['head'].unitsPerEm
        vertical_before = (font['hhea'].ascent,font['hhea'].descent,font['hhea'].lineGap)
        builder = FontBuilder(font=font)
        builder.setupCharacterMap({code:glyph_order[glyph] for code,glyph in mapping.items()})
        alias = f'Figure2Font{index}'
        aliases[name] = alias
        builder.setupNameTable({'familyName':alias, 'styleName':'Regular',
            'uniqueFontIdentifier':alias, 'fullName':name, 'psName':alias,
            'version':'Version 1.0; embedded PDF glyphs'})
        if 'OS/2' not in font:
            builder.setupOS2(sTypoAscender=font['hhea'].ascent,
                sTypoDescender=font['hhea'].descent, sTypoLineGap=font['hhea'].lineGap,
                usWinAscent=max(0,font['head'].yMax), usWinDescent=max(0,-font['head'].yMin))
        if 'post' not in font:
            builder.setupPost()
        assert font.getTableData('glyf') == outlines_before
        assert font['hmtx'].metrics == advances_before
        assert font['head'].unitsPerEm == units_before
        assert (font['hhea'].ascent,font['hhea'].descent,font['hhea'].lineGap) == vertical_before
        font.recalcTimestamp = False
        buffer = BytesIO(); font.save(buffer)
        content = buffer.getvalue()
        (fonts_dir/f'{alias}.ttf').write_bytes(content)
        encoded = base64.b64encode(content).decode('ascii')
        stylesheet.append(f'@font-face{{font-family:"{alias}";src:url(data:font/ttf;base64,{encoded}) format("truetype");font-weight:400;font-style:normal;}}')
        records.append({'PDF_font':name, 'SVG_family':alias, 'unicode_glyph_map':mapping,
                        'original_outlines_and_advances_retained':True})

    changed_strings = 0
    max_origin_error = 0
    for text, span in zip(texts, traces):
        tspans = list(text)
        assert len(tspans) == 1 and tspans[0].tag == f'{{{SVG}}}tspan'
        tspan = tspans[0]
        matrix = np.array([float(v) for v in re.findall(r'[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?',text.get('transform'))])
        assert matrix.shape == (6,)
        x = np.fromstring(tspan.get('x'), sep=' ')
        y = np.fromstring(tspan.get('y'), sep=' ')
        y = np.repeat(y,len(x)) if len(y) == 1 else y
        assert len(x) == len(y) == len(span['chars'])
        origins = np.column_stack([matrix[0]*x+matrix[2]*y+matrix[4],
                                   matrix[1]*x+matrix[3]*y+matrix[5]])
        expected = np.array([char[2] for char in span['chars']])
        error = float(np.max(np.abs(origins-expected)))
        max_origin_error = max(max_origin_error,error)
        assert error < .0002, (span['font'],error)
        expected_size = float(text.get('font-size'))*np.hypot(matrix[0],matrix[1])
        assert abs(expected_size-span['size']) < .0001
        string = ''.join(chr(char[0]) for char in span['chars'])
        changed_strings += string != tspan.text
        tspan.text = string
        text.set('font-family',aliases[span['font']])
        text.set('font-weight','400')
        text.set('style','font-synthesis:none;font-kerning:none;font-variant-ligatures:none')
    style = ET.Element(f'{{{SVG}}}style')
    style.text = '\n'.join(stylesheet)
    tree.insert(0,style)
    width,height = [float(v) for v in tree.get('viewBox').split()][2:]
    # PDF soft masks use page coordinates. SVG's default objectBoundingBox clips
    # horizontal/vertical paths to an empty box and hides the original RNA lines.
    masks = tree.findall(f'.//{{{SVG}}}mask')
    for mask in masks:
        mask.set('maskUnits','userSpaceOnUse')
        mask.set('maskContentUnits','userSpaceOnUse')
        mask.set('x','0'); mask.set('y','0')
        mask.set('width',str(width)); mask.set('height',str(height))
    tree.set('width',f'{width*25.4/72:.12g}mm')
    tree.set('height',f'{height*25.4/72:.12g}mm')
    ET.register_namespace('',SVG)
    ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
    ET.ElementTree(tree).write(destination,encoding='utf-8',xml_declaration=True)
    record = {'live_text_spans':len(texts),'live_characters':sum(len(s['chars']) for s in traces),
              'source_unicode_restored_in_spans':changed_strings,
              'max_text_origin_error_pt':max_origin_error,'fonts':records,
              'page_coordinate_soft_masks':len(masks),
              'text_converted_to_outlines':False}
    (evidence_dir/'svg_font_validation.json').write_text(json.dumps(record,indent=2))
    return record
