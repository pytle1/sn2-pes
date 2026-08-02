# > import of required data structures
from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Method, Sqm
from opi.input.structures import Structure

import html
import py3Dmol
import matplotlib.pyplot as plt

import shutil
import sys
from pathlib import Path
import os

from paths import calc_folder, SAVE_DIR, FIG_DIR, VIEW_DIR, XYZ_DIR, CALC_DIR

print(Path.cwd())

PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)
workdir = calc_folder("initial_complex_hf", clean=True)
workdir

# > initializing calculator
calc = Calculator(basename="job", working_dir=workdir)

# > adding geometry of water from smiles
calc.structure = Structure.from_xyz(XYZ_DIR / "initial_complex.xyz", charge=-1, multiplicity=1)

# > setting method Hartree-Fock Hartree-Fock with basis set def2-TZVP
calc.input.add_simple_keywords(
        Method.HF,
        BasisSet.DEF2_TZVP
)

# > doing calculations
calc.input.ncores = 4
# calc.input.maxcore = 2000

# > writing ORCA input file
calc.write_input()

# > executing orca job
calc.run()

# > get output object
output = calc.get_output()

# > check if orca terminated successfully
if not output.terminated_normally():
    print("exit")
    sys.exit()

# > parse results from orca's JSON files
output.parse()

# > verify that the scf converged
if not output.scf_converged():
    sys.exit()

# > print the final energy in Hartree
output.get_final_energy()
print("Final energy:", output.get_final_energy())

