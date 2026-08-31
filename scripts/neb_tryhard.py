import os
from pathlib import Path

from opi.core import Calculator
from opi.output.core import Output
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import Dft, Task, Scf, BasisSet, SimpleKeyword, Neb, DispersionCorrection, Wft
from opi.input.blocks import BlockNeb, BlockScf

from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
import matplotlib.pyplot as plt



SIMPLE_KEYWORDS = [
        Wft.HF,
        BasisSet.STO-3G,
        Task.FREQ,
]

NCORES = 4


# > Define setup functions for all three calculator types
def setup_opt_calc(basename : str, xyz_file: str, base_dir: Path, ncores: int, sk_list: list[SimpleKeyword]) -> Calculator:
    # > Setup subdirectory for geometry optimization calculations
    working_dir = base_dir / "geometry_optimization"
    working_dir.mkdir(parents=True, exist_ok=True)  # Ensure the subdirectory exists

    # > Specify geometry optimization as task
    sk_list = sk_list.copy()
    sk_list.append(Task.OPT)

    # > Set up a Calculator object
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > Assign structure to calculator
    calc.structure = Structure.from_xyz(base_dir / xyz_file)
    # > Add simple keywords to calculator
    calc.input.add_simple_keywords(*sk_list)
    # > Set the number of CPUs for the calculation
    calc.input.ncores = ncores

    return calc

def setup_neb_calc(basename: str, r_xyz_file: str, p_xyz_file: str, base_dir: Path, ncores: int, sk_list: list[SimpleKeyword]) -> Calculator:
    working_dir = base_dir / "nudged_elastic_band"
    working_dir.mkdir(parents=True, exist_ok=True)
    # > Copy the product XYZ file to the calculation directory so ORCA can find it
    p_xyz_path = base_dir / p_xyz_file
    dest_path = working_dir / p_xyz_path.name
    shutil.copy2(p_xyz_path, dest_path)

    sk_list = sk_list.copy()
    sk_list.append(Neb.NEB_TS)

    calc = Calculator(basename=basename, working_dir=working_dir)
    calc.structure = Structure.from_xyz(base_dir / r_xyz_file)
    # > Define input block for product structure
    calc.input.add_blocks(BlockNeb(neb_end_xyzfile=p_xyz_path.name))
    # > Add options to the calculator
    calc.input.add_simple_keywords(*sk_list)
    calc.input.ncores = ncores

    return calc

def setup_irc_calc(basename: str, xyz_file: str, base_dir: Path, ncores: int, sk_list: list[SimpleKeyword]) -> Calculator:
    working_dir = base_dir / "intrinsic_reaction_coordinate"
    working_dir.mkdir(parents=True, exist_ok=True)  # Ensure the subdirectory exists

    sk_list = sk_list.copy()
    sk_list.append(Neb.IRC)

    calc = Calculator(basename=basename, working_dir=working_dir)
    calc.structure = Structure.from_xyz(base_dir / xyz_file)
    calc.input.add_blocks(BlockScf(maxiter=200))
    calc.input.add_simple_keywords(*sk_list)
    calc.input.ncores = ncores

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
        raise RuntimeError(f"ORCA did not terminate normally. Please check {output.basename}.out file")
    else:
        # > ORCA did terminate normally so we can parse the output
        output.parse()
    # > Now check for convergence of the SCF
    if output.results_properties.geometries[-1].single_point_data.converged:
        print(f"{output.basename}: SCF CONVERGED")
    else:
        raise RuntimeError("SCF DID NOT CONVERGE")

r_opt_calc = setup_opt_calc("reactant", "reactant.xyz", base_dir, ncores=NCORES, sk_list=SIMPLE_KEYWORDS)
r_opt_out = run_calc(r_opt_calc)
check_and_parse_output(r_opt_out)

p_opt_calc = setup_opt_calc("product", "product.xyz", base_dir, ncores=NCORES, sk_list=SIMPLE_KEYWORDS)
p_opt_out = run_calc(p_opt_calc)
check_and_parse_output(p_opt_out)

neb_ts_calc = setup_neb_calc("ts", "geometry_optimization/reactant.xyz", "geometry_optimization/product.xyz", base_dir, ncores=NCORES, sk_list=SIMPLE_KEYWORDS)
neb_ts_out = run_calc(neb_ts_calc)
check_and_parse_output(neb_ts_out)

irc_calc = setup_irc_calc("irc", "nudged_elastic_band/ts.xyz", base_dir, ncores=NCORES, sk_list=SIMPLE_KEYWORDS)
irc_out = run_calc(irc_calc)
check_and_parse_output(irc_out)









