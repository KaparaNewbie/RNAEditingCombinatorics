from pathlib import Path
import subprocess


def multiqc(multiqc_path: Path, data_dir: Path):
    # `--force` is needed to overwrite the existing report
    subprocess.run(f"{multiqc_path} --force {data_dir}", shell=True, cwd=str(data_dir))
