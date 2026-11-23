"""Write scoped S1–S4 numeric cells while preserving the author workbook template."""
from pathlib import Path
import math
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET
import openpyxl
import numpy as np
import pandas as pd
import copy


def export_sample_sd_workbook(source, destination, summary):
    source, destination = Path(source), Path(destination)
    assert source.resolve() != destination.resolve()
    assert not destination.exists(), destination
    workbook = openpyxl.load_workbook(source)
    sheet = workbook['Table S2']
    assert sheet['Q2'].value in {
        'STD Of Recoding Substitutions Per Read (Minimal Estimate)',
        'SD of number of recoded AAs per read (minimum estimate)'}
    assert summary.Gene.is_unique and len(summary) == 19
    summary = summary.set_index('Gene')
    updates = {}
    for row in range(3, 22):
        gene = sheet.cell(row, 3).value
        value = float(summary.loc[gene, 'SD_ddof1_sample'])
        assert math.isfinite(value) and value >= 0
        assert math.isclose(sheet.cell(row, 16).value, summary.loc[gene, 'Mean'], abs_tol=5e-6)
        if not math.isclose(sheet.cell(row, 17).value, value, abs_tol=1e-12, rel_tol=0):
            updates[f'Q{row}'] = value
    if not updates:
        shutil.copy2(source, destination)
        assert source.read_bytes() == destination.read_bytes()
        return []
    with zipfile.ZipFile(source) as archive:
        ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        book = ET.fromstring(archive.read('xl/workbook.xml'))
        sheet_id = next(s for s in book.find('m:sheets', ns) if s.get('name') == 'Table S2').get(
            '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        rels = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        target = next(r.get('Target') for r in rels if r.get('Id') == sheet_id)
        member = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        xml = archive.read(member).decode()
        for address, value in updates.items():
            pattern = rf'(<c\b[^>]*\br="{address}"[^>]*>\s*<v>)[^<]+(</v>\s*</c>)'
            xml, count = re.subn(pattern, lambda m: m[1] + repr(value) + m[2], xml)
            assert count == 1, (address, 'Expected one existing numeric SD cell')
        with zipfile.ZipFile(destination, 'w') as output:
            for item in archive.infolist():
                output.writestr(item, xml.encode() if item.filename == member else archive.read(item.filename))
    with zipfile.ZipFile(source) as old, zipfile.ZipFile(destination) as new:
        assert old.namelist() == new.namelist()
        assert all(old.read(name) == new.read(name) for name in old.namelist() if name != member)
    return list(updates)


def export_fisher_workbook(source, destination, summary):
    """Update S4 numeric cells/caches; retain formulas, styles and other sheets."""
    workbook, destination = Path(source), Path(destination)
    assert workbook.resolve() != destination.resolve() and not destination.exists()
    assert len(summary) == 21 and not summary.duplicated(['Gene', 'ReadType']).any()
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(workbook) as archive:
        strings = [''.join(t.text or '' for t in node.findall('.//m:t', ns))
                   for node in ET.fromstring(archive.read('xl/sharedStrings.xml'))]
        xml = archive.read('xl/worksheets/sheet4.xml').decode('utf-8')
        sheet = ET.fromstring(xml)
        cell_values = {}
        for cell in sheet.findall('.//m:c', ns):
            value = cell.find('m:v', ns)
            text = value.text if value is not None else ''
            cell_values[cell.get('r')] = strings[int(text)] if cell.get('t') == 's' else text
        headers = [cell_values[f'{letter}1'] for letter in 'ABCDEFGHIJKLM']
        expected_headers = ['Gene', 'ReadType', 'NumOfEditingSites', 'NumOfUntestableSites', 'NumOfTestableSites',
            'NumOfSignificantSites', 'NumOfSignificantSitesHigherInRare', 'NumOfSignificantSitesLowerInRare',
            'NumOfInsignificantSites', '%OfSignificantSites', '%OfSignificantSitesHigherInRare',
            '%OfSignificantSitesLowerInRare', '%OfInsignificantSites']
        assert headers == expected_headers
        replacements = {}
        for row in range(2, 23):
            match = summary.loc[summary.Gene.eq(cell_values[f'A{row}']) & summary.ReadType.eq(cell_values[f'B{row}'])]
            assert len(match) == 1
            for letter, label in zip('CDEFGHIJKLM', headers[2:]):
                replacements[f'{letter}{row}'] = float(match.iloc[0][label])
        for letter, label in zip('JKLM', headers[9:]):
            replacements[f'{letter}24'] = float(summary[label].mean())
            replacements[f'{letter}25'] = float(summary[label].std(ddof=0))
        updated = xml
        for address, value in replacements.items():
            assert np.isfinite(value)
            pattern = rf'(<c\b[^>]*\br="{address}"[^>]*>.*?<v>)([^<]*)(</v>)'
            rendered = format(value, '.17g')
            updated, count = re.subn(pattern, lambda m: m[1] + rendered + m[3], updated, flags=re.S)
            assert count == 1, address
            cell_values[address] = rendered
        with zipfile.ZipFile(destination, 'w') as output:
            for member in archive.infolist():
                output.writestr(copy.copy(member), updated.encode('utf-8') if member.filename == 'xl/worksheets/sheet4.xml'
                                else archive.read(member.filename))
        with zipfile.ZipFile(destination) as output:
            assert output.namelist() == archive.namelist()
            unchanged = [n for n in archive.namelist() if n != 'xl/worksheets/sheet4.xml']
            assert all(archive.read(n) == output.read(n) for n in unchanged)
            # Removing only permitted values must leave S4 byte-identical, including formulas/styles.
            def neutralize(text):
                for address in replacements:
                    text = re.sub(rf'(<c\b[^>]*\br="{address}"[^>]*>.*?<v>)[^<]*(</v>)',
                                  lambda m: m[1] + 'VALUE' + m[2], text, flags=re.S)
                return text
            assert neutralize(xml) == neutralize(updated)
    return replacements


def export_pacbio_metrics_workbook(source, destination, metrics, assignments_dir):
    """Write native S1/S3 metrics into the current author template, preserving XML.

S1's two historical filtered-read columns both refer to the ORF population.
The accompanying native TSV also reports all mapped reads separately. S1 recoding
mean/SD retain their documented six-decimal display; other new values are unrounded.
Primers, sample metadata, merged cells, notes, formulas and other sheets are retained.
    """
    from manuscript_inputs import selections, read_assignment
    source, destination = Path(source), Path(destination)
    assert source.resolve() != destination.resolve() and not destination.exists()
    assert not metrics.duplicated(['Dataset', 'Gene', 'Replicate']).any() and len(metrics) == 14
    workbook = openpyxl.load_workbook(source, data_only=True)
    updates = {'Table S1': {}, 'Table S3': {}}
    s1 = workbook['Table S1']
    assert [s1.cell(2, c).value for c in [3, 4, 18, 22, 23, 24]] == [
        'Gene', 'Transcript ID', 'Total CCS Reads', 'Reads retained for pileup analysis',
        'Unique Reads (unique editing status)', 'Distinct encoded protein isoforms ']
    mapping = {'E': 'Start', 'F': 'End', 'G': 'ORFAminoAcids', 'K': 'KnownSites', 'L': 'DetectedSites',
        'M': 'MeanEditedSitesPerRead', 'N': 'SDEditedSitesPerRead', 'O': 'EditableAAs',
        'P': 'MeanMinimalRecodingSubstitutionsPerRead', 'Q': 'SDMinimalRecodingSubstitutionsPerRead',
        'R': 'TotalCCSReads', 'S': 'AlignedORFReads', 'T': 'FilteredORFReads', 'U': 'FilteredORFReads',
        'V': 'PileupReads', 'W': 'UniqueReads', 'X': 'DistinctProteins'}
    for row, gene in [(3, 'GRIA2'), (4, 'PCLO')]:
        data = metrics.loc[metrics.Dataset.eq('PB1') & metrics.Gene.eq(gene)].iloc[0]
        assert s1.cell(row, 3).value == gene and s1.cell(row, 4).value == data.Chrom
        for letter, column in mapping.items():
            value = float(data[column])
            updates['Table S1'][f'{letter}{row}'] = round(value, 6) if letter in ['P', 'Q'] else value
        for letter, metric in [('Y', 'Grantham100'), ('Z', 'Miyata')]:
            selected = selections().loc[lambda t: t.Platform.eq('PacBio') & t.Gene.eq(gene) & t.Metric.eq(metric)].iloc[0]
            updates['Table S1'][f'{letter}{row}'] = len(read_assignment(selected, assignments_dir))
    s3 = workbook['Table S3']
    assert [s3.cell(1, c).value for c in [3, 4, 5, 24, 25]] == [
        'Gene', 'Replicate', 'Transcript ID', 'Reads retained for pileup analysis', 'Unique Reads (unique editing status)']
    mapping = {'P': 'MeanEditedSitesPerRead', 'Q': 'SDEditedSitesPerRead',
        'S': 'MeanMinimalRecodingSubstitutionsPerRead', 'T': 'SDMinimalRecodingSubstitutionsPerRead',
        'U': 'AlignedORFReads', 'V': 'FilteredMappedReads', 'W': 'FilteredORFReads',
        'X': 'PileupReads', 'Y': 'UniqueReads'}
    for start, gene in [(2, 'GRIA2'), (6, 'ADAR1'), (10, 'IQEC1')]:
        assert s3.cell(start, 3).value == gene
        for offset, replicate in enumerate(['1', '2', '3', 'Pooled']):
            row = start + offset
            assert str(s3.cell(row, 4).value) == replicate
            data = metrics.loc[metrics.Dataset.eq('PB3') & metrics.Gene.eq(gene) & metrics.Replicate.astype(str).eq(replicate)].iloc[0]
            for letter, column in mapping.items(): updates['Table S3'][f'{letter}{row}'] = float(data[column])
            if replicate == 'Pooled':
                for letter, column in [('N', 'KnownSites'), ('O', 'DetectedSites'), ('R', 'EditableAAs'), ('Z', 'DistinctProteins')]:
                    if pd.notna(data[column]): updates['Table S3'][f'{letter}{row}'] = float(data[column])
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(source) as archive:
        book = ET.fromstring(archive.read('xl/workbook.xml'))
        rels = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        rewritten = {}
        for sheet_name, values in updates.items():
            sheet_id = next(s for s in book.find('m:sheets', ns) if s.get('name') == sheet_name).get(
                '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
            target = next(r.get('Target') for r in rels if r.get('Id') == sheet_id)
            member = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            xml = archive.read(member).decode()
            for address, value in values.items():
                assert math.isfinite(value)
                pattern = rf'(<c\b[^>]*\br="{address}"[^>]*>.*?<v>)[^<]*(</v>)'
                xml, count = re.subn(pattern, lambda m: m[1] + format(value, '.17g') + m[2], xml, flags=re.S)
                assert count == 1, (sheet_name, address)
            rewritten[member] = xml.encode()
        with zipfile.ZipFile(destination, 'w') as output:
            for member in archive.infolist(): output.writestr(copy.copy(member), rewritten.get(member.filename, archive.read(member.filename)))
        with zipfile.ZipFile(destination) as output:
            assert output.namelist() == archive.namelist()
            assert all(output.read(n) == archive.read(n) for n in archive.namelist() if n not in rewritten)
    return updates
