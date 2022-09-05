# helper functions to set args' types

from pathlib import Path


def abs_path_from_str(rel_str_path: str) -> Path:
    return Path(rel_str_path).absolute()


def expanded_path_from_str(unexpanded_str_path: str) -> Path:
    return Path(unexpanded_str_path).expanduser()
