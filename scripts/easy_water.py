# > import of required data structures
from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Method
from opi.input.structures import Structure
import os
import sys
from pathlib import Path

print(Path.cwd())

PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)
print(Path.cwd())

from paths import CALC_DIR, FIG_DIR, XYZ_DIR, calc_folder

workdir = calc_folder("easy_water", clean=True)
workdir

# > initializing calculator
calc = Calculator(basename="job", working_dir=workdir)

# > adding geometry of water from smiles
water_smiles = Structure.from_smiles("O", charge=0, multiplicity=1)
calc.structure = water_smiles

# > setting method Hartree-Fock Hartree-Fock with basis set def2-TZVP
calc.input.add_simple_keywords(
        Method.HF,
        BasisSet.DEF2_TZVP
)

# > doing calculations
calc.input.ncores = 4
calc.input.maxcore = 2000

# > writing ORCA input file
calc.write_input()

# > executing orca job
calc.run()

# > get output object
output = calc.get_output()

# > check if orca terminated successfully
if not output.terminated_normally():
    print("terminated")
    sys.exit()

# > parse results from orca's JSON files
output.parse()

# > verify that the scf converged
if not output.scf_converged():
    print("not converged")
    sys.exit()

# > print the final energy in Hartree
output.get_final_energy()
print("Final energy:", output.get_final_energy())

