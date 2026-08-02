from pathlib import Path
import shutil

PROJECT = Path.cwd()

SAVE_DIR = PROJECT / "results"
FIG_DIR = PROJECT / "figures"
CALC_DIR = PROJECT / "calculations"
VIEW_DIR = FIG_DIR / "molecules"
XYZ_DIR = PROJECT / "molecules"

for d in [SAVE_DIR, FIG_DIR, CALC_DIR, VIEW_DIR, XYZ_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def calc_folder(name, clean=False):
    folder = CALC_DIR / name

    if clean:
        shutil.rmtree(folder, ignore_errors=True)

    folder.mkdir(parents=True, exist_ok=True)
    return folder
