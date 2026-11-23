"""Container-only PDF composition, vector/text SVG export and exact PDF checks."""
from pathlib import Path
import csv, hashlib, json, math, re, sys
import xml.etree.ElementTree as ET
import numpy as np
import pymupdf as fitz
from PIL import Image
from live_svg import export_live_svg

ROOT = None
EXPORTS = None
CONFIG = {'sources': {'reference': {'file': 'inputs/BioRender_reference.pdf',
                           'sha256': '488f0c583be6c2c43a9aab1fa757284d20784929d783a7a4eca364805e9b1af0'},
             'base': {'file': 'inputs/BioRender_base.pdf',
                      'sha256': '14b9bbe48488f80cc4aa4612e129472e340699ccccd5cadfc956b048a2625e44'},
             'panel': {'file': 'inputs/Figure2b_generated.pdf'}},
 'page': {'MediaBox_pt': [0, 0, 510, 629.04],
          'CropBox_pt': [0, 0, 510, 629.04],
          'intended_mm': [180, 222],
          'preserve_actual_boxes': True},
 'placement': {'coordinate_system': 'PDF page coordinates in points, origin '
                                    'TOP LEFT, x right, y down',
               'top_left_pt': [367.08, 10.065],
               'size_mm': [49, 45.5],
               'whole_source_page': True,
               'rotation_degrees': 0,
               'preserve_aspect_ratio': True},
 'reference_panel_bbox_pt': [367.08, 10.065, 506.16, 139.185],
 'container_image': 'combinatorics-figure2:2026-10-13',
 'container_image_id': 'sha256:e1fcbeb34fbe87a26e76d553e8987277b29e3792eba52b7555f96fc1d5a47c97'}
PT_PER_MM = 72/25.4


def read_sources():
    documents = {}
    for role, entry in CONFIG['sources'].items():
        path = ROOT/entry['file']
        expected = (json.loads((ROOT/'evidence/panel_source.json').read_text())['sha256']
                    if role == 'panel' else entry['sha256'])
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected,role
        doc = fitz.open(path)
        assert len(doc)==1 and doc[0].rotation==0
        documents[role]=doc
    return documents


def placement_rect():
    x,y=CONFIG['placement']['top_left_pt']
    w,h=np.array(CONFIG['placement']['size_mm'])*PT_PER_MM
    return fitz.Rect(x,y,x+w,y+h)


def assemble():
    documents=read_sources()
    base,panel,reference=[documents[k] for k in ['base','panel','reference']]
    np.testing.assert_allclose(base[0].mediabox, CONFIG['page']['MediaBox_pt'],atol=.0001,rtol=0)
    np.testing.assert_allclose(base[0].cropbox, CONFIG['page']['CropBox_pt'],atol=.0001,rtol=0)
    np.testing.assert_allclose([panel[0].rect.width,panel[0].rect.height],
        np.array([49,45.5])*PT_PER_MM,atol=.0001,rtol=0)
    box=placement_rect()
    assert not panel[0].get_image_info()
    assert len(base[0].get_image_info())==6 and len(reference[0].get_image_info())==7
    assert all(not box.intersects(fitz.Rect(item['bbox'])) for item in base[0].get_image_info())
    assert sum(word[4]=='b' for word in base[0].get_text('words'))==1
    base[0].show_pdf_page(box,panel,0,keep_proportion=True,overlay=True)
    base.save(EXPORTS/'Figure2.pdf',garbage=3,deflate=True,no_new_id=True)
    (ROOT/'evidence/placement.json').write_text(json.dumps({
        'rect_pt_top_left':list(box),'rect_mm_top_left':(np.array(box)/PT_PER_MM).tolist(),
        'size_mm':[49,45.5],'uniform_scale':1,
        'whole_asset_page_placed':True,'b_letter_added':False},indent=2))
    print('Assembled vector PDF:',list(box),flush=True)


def export():
    document=fitz.open(EXPORTS/'Figure2.pdf')
    export_live_svg(document,EXPORTS/'Figure2.svg',ROOT/'fonts',ROOT/'evidence')
    pixmap=document[0].get_pixmap(dpi=600,alpha=False,colorspace=fitz.csRGB)
    assert (pixmap.width,pixmap.height)==(4250,5242)
    pixmap.save(EXPORTS/'Figure2.png')
    document[0].get_pixmap(dpi=150,alpha=False).save(ROOT/'previews/Figure2_review.png')
    document[0].get_pixmap(dpi=300,alpha=False).save(ROOT/'previews/Figure2_PDF_300dpi.png')
    print('Exported live-text SVG and full-page white-background PNG:',pixmap.width,pixmap.height,flush=True)


def pixels(page,dpi=300):
    pix=page.get_pixmap(dpi=dpi,alpha=False,colorspace=fitz.csRGB)
    return np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,3)


def validate():
    docs=read_sources()
    output=fitz.open(EXPORTS/'Figure2.pdf')
    page,base,panel,reference=output[0],docs['base'][0],docs['panel'][0],docs['reference'][0]
    assert tuple(page.mediabox)==tuple(base.mediabox) and tuple(page.cropbox)==tuple(base.cropbox)
    for key in ['MediaBox','CropBox']:
        assert output.xref_get_key(page.xref,key)==docs['base'].xref_get_key(base.xref,key)
    bbox=placement_rect()
    # The adopted clean base is the comparison authority outside the replaced panel.
    original_image=pixels(base);final_image=pixels(page);reference_image=pixels(reference)
    outside=np.ones(original_image.shape[:2],dtype=bool)
    original_box=CONFIG['reference_panel_bbox_pt']
    union=fitz.Rect(original_box)|bbox
    x0,y0=np.floor(np.array(union[:2])*300/72).astype(int)-1
    x1,y1=np.ceil(np.array(union[2:])*300/72).astype(int)+1
    outside[y0:y1,x0:x1]=False
    changed=np.any(original_image!=final_image,axis=2)&outside
    assert not changed.any(),int(changed.sum())
    reference_diff=np.any(reference_image!=original_image,axis=2)&outside
    # Persist inherited base/reference differences, never attribute them to assembly.
    Image.fromarray(np.where(reference_diff[:,:,None],original_image,255).astype('uint8')).save(
        ROOT/'previews/inherited_base_reference_difference.png')
    assert page.get_text().split()==base.get_text().split()+panel.get_text().split()
    assert sum(w[4]=='b' for w in page.get_text('words'))==1
    panel_spans=panel.get_texttrace();final_spans=page.get_texttrace()[-len(panel_spans):]
    for old,new in zip(panel_spans,final_spans):
        assert old['font']==new['font'] and abs(old['size']-new['size'])<.0001
        assert [v[:2] for v in old['chars']]==[v[:2] for v in new['chars']]
        np.testing.assert_allclose(np.array([v[2] for v in new['chars']]),
            np.array([v[2] for v in old['chars']])+bbox[:2],atol=.0001,rtol=0)
        for char in new['chars']:
            assert bbox.contains(fitz.Rect(char[3]))
    before=base.get_drawings();after=page.get_drawings()
    assert len(after)==len(before)+len(panel.get_drawings())
    assert after[:len(before)]==before
    # show_pdf_page carries the original page content stream as a Form XObject.
    original_stream=b'\n'.join(docs['panel'].xref_stream(x) for x in panel.get_contents())
    forms=page.get_xobjects()
    assert any(output.xref_stream(row[0])==original_stream for row in forms)
    green_black_gray=[]
    for drawing in panel.get_drawings():
        if len(drawing['items'])>100:
            green_black_gray.append(sum(item[0]=='re' for item in drawing['items']))
    assert sum(green_black_gray)==14400,green_black_gray

    image_records=[]
    base_images=base.get_image_info(hashes=True)
    final_images=page.get_image_info(hashes=True)
    assert len(base_images)==len(final_images)==6
    labels=[('a–d','Curved connection arrow'),('c','Pairwise relation mark 1'),
            ('c','Pairwise relation mark 2'),('c','Pairwise relation mark 3'),
            ('e–f','Curved assignment arrow'),('d–e','Horizontal arrow')]
    for old,new,(panel_name,label) in zip(base_images,final_images,labels):
        assert old['digest']==new['digest'] and old['bbox']==new['bbox'] and old['transform']==new['transform']
        assert not bbox.intersects(fitz.Rect(new['bbox']))
        a,b,c,d,e,f=new['transform']
        ppi_x=new['width']*72/math.hypot(a,b);ppi_y=new['height']*72/math.hypot(c,d)
        np.testing.assert_allclose([ppi_x,ppi_y],[600,600],atol=.002,rtol=0)
        image_records.append({'panel':panel_name,'component':label,'width_px':new['width'],
            'height_px':new['height'],'effective_ppi_x':ppi_x,'effective_ppi_y':ppi_y,
            'bbox_pt':json.dumps(new['bbox']),'pixels_preserved':True})
    with (ROOT/'reports/remaining_raster.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(image_records[0]),delimiter='\t')
        writer.writeheader();writer.writerows(image_records)
    fonts=[]
    for record in page.get_fonts(full=True):
        name,extension,kind,data=output.extract_font(record[0])
        assert extension=='ttf' and len(data)>0
        fonts.append(name)
    tree=ET.parse(EXPORTS/'Figure2.svg').getroot()
    ns={'s':'http://www.w3.org/2000/svg'}
    svg_text=[''.join(node.itertext()) for node in tree.findall('.//s:text',ns)]
    trace_text=[''.join(chr(c[0]) for c in span['chars']) for span in page.get_texttrace()]
    assert svg_text==trace_text
    assert tree.find('s:style',ns).text.count('@font-face')==3
    matrix_paths=[]
    for path in tree.findall('.//s:path',ns):
        rectangles=len(re.findall(r'[Mm]',path.get('d','')))
        if rectangles>100:matrix_paths.append((path.get('fill','#000000'),rectangles))
    assert sorted(matrix_paths)==[('#000000',10294),('#4af626',1851),('#c0c0c0',2255)]
    assert len(tree.findall('.//s:image',ns))==12  # six original images plus alpha masks
    with Image.open(EXPORTS/'Figure2.png') as image:
        assert image.size==(4250,5242) and image.mode=='RGB'
        assert abs(image.info['dpi'][0]-600)<.02
    record={'MediaBox_pt':list(page.mediabox),'CropBox_pt':list(page.cropbox),
        'page_mm':[510*25.4/72,629.04*25.4/72], 'panel_mm':[49,45.5],
        'placement_pt':list(bbox),'outside_2b_changed_pixels_300dpi':int(changed.sum()),
        'inherited_base_vs_reference_changed_pixels_300dpi':int(reference_diff.sum()),
        'all_original_panel_content_stream_bytes_preserved':True,
        'matrix_rectangles':sum(green_black_gray),'matrix_raster_images':0,
        'SVG_matrix_vector_subpaths':dict(matrix_paths),'SVG_existing_images_and_masks':12,
        'live_SVG_text_spans':len(svg_text),'live_characters':sum(map(len,svg_text)),
        'embedded_PDF_fonts':sorted(set(fonts)),'other_raster_components':image_records,
        'PNG_pixels':[4250,5242],'PNG_DPI':600,'panel_letters_b':1,
        'scientific_recomputation':False}
    (ROOT/'evidence/validation.json').write_text(json.dumps(record,indent=2))
    print('Passed PDF, page, placement, all vector cells, fonts, live text and unchanged-base checks.',flush=True)


if __name__=='__main__':
    ROOT = Path(sys.argv[2])
    EXPORTS = ROOT / 'exports'
    assert ROOT.is_dir()
    {'assemble':assemble,'export':export,'validate':validate}[sys.argv[1]]()
