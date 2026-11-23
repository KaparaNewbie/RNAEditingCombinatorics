"""Bounded per-figure export using the existing, versioned rendering image."""
from pathlib import Path
import os
import subprocess

IMAGE = 'combinatorics-figure2:2026-10-13'
IMAGE_ID = 'sha256:e1fcbeb34fbe87a26e76d553e8987277b29e3792eba52b7555f96fc1d5a47c97'


def export_cropped_page(source_dir, output_dir, name, rectangle, font_dir):
    """Write one new figure directory. No approved collection is replaced."""
    source_dir, output_dir, font_dir = map(lambda p: Path(p).resolve(), [source_dir, output_dir, font_dir])
    assert not output_dir.exists(), output_dir
    assert all((source_dir / ext / f'{name}.{ext}').is_file() for ext in ['pdf', 'svg'])
    installed = subprocess.check_output(['docker', 'image', 'inspect', '--format', '{{.Id}}', IMAGE], text=True).strip()
    assert installed == IMAGE_ID
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    code = Path(__file__).resolve().parent
    command = ['docker', 'run', '--rm', '--network', 'none', '--cpus', '2', '--memory', '4g',
        '--pids-limit', '128', '--read-only', '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges',
        '--tmpfs', '/tmp:rw,size=256m', '--user', f'{os.getuid()}:{os.getgid()}',
        '-v', f'{code}:/source:ro', '-v', f'{source_dir}:/inputs:ro',
        '-v', f'{font_dir}:/fonts:ro', '-v', f'{output_dir.parent}:/outputs',
        IMAGE, '/source/page_export.py', '--source', '/inputs', '--output', f'/outputs/{output_dir.name}',
        '--name', name, '--crop', *map(str, rectangle), '--fonts', '/fonts']
    subprocess.run(command, check=True, timeout=360)
    return output_dir


def prepare_figure2(source_dir, work_dir, panel_pdf):
    """Bind external artwork and this run's native-data panel in a new workspace."""
    import hashlib
    import json
    import shutil
    source_dir, work_dir, panel_pdf = map(Path, [source_dir, work_dir, panel_pdf])
    hashes = {'BioRender_base.pdf': '14b9bbe48488f80cc4aa4612e129472e340699ccccd5cadfc956b048a2625e44',
        'BioRender_reference.pdf': '488f0c583be6c2c43a9aab1fa757284d20784929d783a7a4eca364805e9b1af0'}
    for name, digest in hashes.items():
        assert hashlib.sha256((source_dir / name).read_bytes()).hexdigest() == digest, name
    assert panel_pdf.is_file(), 'Render Figure2b from this run\'s PCLO table first'
    work_dir.mkdir(parents=True, exist_ok=False)
    for folder in ['inputs', 'exports', 'evidence', 'previews', 'reports', 'fonts']:
        (work_dir / folder).mkdir()
    for name in hashes:
        shutil.copy2(source_dir / name, work_dir / 'inputs' / name)
    shutil.copy2(panel_pdf, work_dir / 'inputs/Figure2b_generated.pdf')
    (work_dir / 'evidence/panel_source.json').write_text(json.dumps({
        'source': str(panel_pdf.resolve()),
        'sha256': hashlib.sha256(panel_pdf.read_bytes()).hexdigest()}, indent=2))


def run_figure2_stage(work_dir, stage):
    assert stage in ['assemble', 'export', 'validate']
    assert subprocess.check_output(['docker', 'image', 'inspect', '--format', '{{.Id}}', IMAGE], text=True).strip() == IMAGE_ID
    code = Path(__file__).resolve().parent
    work_dir = Path(work_dir).resolve()
    command = ['docker', 'run', '--rm', '--network', 'none', '--cpus', '2', '--memory', '4g',
        '--pids-limit', '128', '--read-only', '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges',
        '--tmpfs', '/tmp:rw,size=256m', '--user', f'{os.getuid()}:{os.getgid()}',
        '-v', f'{code}:/source:ro', '-v', f'{work_dir}:/work', IMAGE,
        '/source/figure2_assembly.py', stage, '/work']
    result = subprocess.run(command, capture_output=True, text=True, timeout=180)
    (work_dir / 'evidence' / f'{stage}.log').write_text(result.stdout + result.stderr)
    result.check_returncode()


def copy_new_figure(source_dir, destination, name):
    """Copy completed figure artifacts only when all destination files are absent."""
    import shutil
    source_dir, destination = Path(source_dir), Path(destination)
    pairs = [(source_dir / f'{name}.{ext}', destination / folder / f'{name}.{ext}')
        for ext, folder in [('pdf', 'pdf'), ('svg', 'svg'), ('png', 'previews')]]
    assert all(src.is_file() and not dst.exists() for src, dst in pairs)
    for src, dst in pairs:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
