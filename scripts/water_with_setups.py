from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import Dft, Task
from opi.input.structures.structure import Structure
import os
import sys
from pathlib import Path

def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 1) -> Calculator:
    # > Set up a Calculator object, the basename for the ORCA calculation is also set
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > Assign structure to calculator
    calc.structure = structure

    # Define a level of theory: we use r2SCAN-3c for geometry optimization
    sk_list = [
        Dft.R2SCAN_3C, # > r2SCAN-3c method (Comes with a predefined basis set)
        Task.SP # > Perform the singlepoint
    ]

    # > Use simple keywords in calculator
    calc.input.add_simple_keywords(*sk_list)

    # > Define number of CPUs for the calcualtion
    calc.input.ncores = ncores # > CPUs for this ORCA run (default: 1)

    return calc


def run_calc(calc: Calculator) -> Output:
    # > Write the ORCA input file
    calc.write_input()
    # > Run the ORCA calculation
    print("Running ORCA calculation ...", end="")
    calc.run()
    print("   Done")

    # > Get the output object
    output = calc.get_output()
    
    return output

def check_and_parse_output(output: Output):
    # > Check for proper termination of ORCA
    status = output.terminated_normally()
    if not status:
        # > ORCA did not terminate normally
        raise RuntimeError("ORCA did not terminate normally. Please check RUN/single_point.out")
    else:
        # > ORCA did terminate normally so we can parse the output
        output.parse()

    # Now check for convergence of the SCF
    if output.results_properties.geometries[0].single_point_data.converged:
        print("SCF CONVERGED")
    else:
        raise RuntimeError("SCF DID NOT CONVERGE")


PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)
print(Path.cwd())

from paths import CALC_DIR, FIG_DIR, XYZ_DIR, calc_folder

workdir = calc_folder("easy_water_with_setups", clean=True)
workdir

# > adding geometry of water from smiles
structure = Structure.from_smiles("O", charge=0, multiplicity=1)

calc = setup_calc("single_point", working_dir=workdir, structure=structure)
output = run_calc(calc)
check_and_parse_output(output)

# > Total energy in Hartee
Etot = output.results_properties.geometries[0].single_point_data.finalenergy

# > Print total energy
print("Final energy in Hartee: ",Etot)

