# %%
from pathlib import Path
import os
from opi.input.structures.structure import Structure
PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)

SAVE_DIR = PROJECT / "results"
FIG_DIR = PROJECT / "figures"
CALC_DIR = PROJECT / "calculations"

SAVE_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)
CALC_DIR.mkdir(exist_ok=True)

print("Project:", PROJECT)
print("Current dir:", Path.cwd())

structure = Structure.from_xyz(SAVE_DIR / "water.xyz")
