from pathlib import Path
import shutil
import os

from opi.core import Calculator
from opi.output.core import Output
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import (
    Wft,
    Task,
    BasisSet,
    SimpleKeyword,
    Neb,
)
from opi.input.blocks import BlockNeb, BlockScf

from paths import (
    calc_folder,
    SAVE_DIR,
    FIG_DIR,
    VIEW_DIR,
    XYZ_DIR,
    CALC_DIR,
)

PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)

NCORES = 4

SIMPLE_KEYWORDS = [
    Wft.HF,
    BasisSet.STO_3G,
]

REACTION_DIR = calc_folder("neb_reaction_ordered", clean=True)

GEOM_DIR = REACTION_DIR / "geometry_optimization"
NEB_DIR = REACTION_DIR / "nudged_elastic_band"
IRC_DIR = REACTION_DIR / "intrinsic_reaction_coordinate"

for d in [GEOM_DIR, NEB_DIR, IRC_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def setup_opt_calc(
    name: str,
    xyz_path: Path,
    ncores: int,
    sk_list: list[SimpleKeyword],
) -> Calculator:
    working_dir = GEOM_DIR / name
    working_dir.mkdir(parents=True, exist_ok=True)

    keywords = sk_list.copy()
    keywords.append(Task.OPT)

    calc = Calculator(
        basename=name,
        working_dir=working_dir,
    )

    calc.structure = Structure.from_xyz(xyz_path, charge=-1, multiplicity=1)
    calc.input.add_simple_keywords(*keywords)
    calc.input.ncores = ncores

    return calc

def setup_neb_calc(
    basename: str,
    reactant_xyz: Path,
    product_xyz: Path,
    ncores: int,
    sk_list: list[SimpleKeyword],
) -> Calculator:
    NEB_DIR.mkdir(parents=True, exist_ok=True)

    product_dest = NEB_DIR / "final_complex.xyz"
    shutil.copy2(product_xyz, product_dest)

    keywords = sk_list.copy()
    keywords.append(Neb.NEB_TS)

    calc = Calculator(
        basename=basename,
        working_dir=NEB_DIR,
    )

    calc.structure = Structure.from_xyz(reactant_xyz, charge=-1, multiplicity=1)

    calc.input.add_blocks(
        BlockNeb(neb_end_xyzfile="final_complex.xyz")
    )

    calc.input.add_simple_keywords(*keywords)
    calc.input.ncores = ncores

    return calc

def setup_irc_calc(
    basename: str,
    ts_xyz: Path,
    ncores: int,
    sk_list: list[SimpleKeyword],
) -> Calculator:
    IRC_DIR.mkdir(parents=True, exist_ok=True)

    keywords = sk_list.copy()
    keywords.append(Neb.IRC)

    calc = Calculator(
        basename=basename,
        working_dir=IRC_DIR,
    )

    calc.structure = Structure.from_xyz(ts_xyz, charge=-1, multiplicity=1)
    calc.input.add_blocks(BlockScf(maxiter=200))
    calc.input.add_simple_keywords(*keywords)
    calc.input.ncores = ncores

    return calc

def run_calc(calc: Calculator) -> Output:
    calc.write_input()

    print(f"Running {calc.basename} in {calc.working_dir} ...")
    calc.run()
    print("Done")

    return calc.get_output()


def check_and_parse_output(output: Output):
    if not output.terminated_normally():
        raise RuntimeError(
            f"ORCA did not terminate normally. Check {output.basename}.out"
        )

    output.parse()

    if output.results_properties.geometries[-1].single_point_data.converged:
        print(f"{output.basename}: SCF CONVERGED")
    else:
        raise RuntimeError(f"{output.basename}: SCF DID NOT CONVERGE")

reactant_xyz = XYZ_DIR / "initial_complex_ordered.xyz"
product_xyz = XYZ_DIR / "final_complex_ordered.xyz"

r_opt_calc = setup_opt_calc(
    name="initial_complex",
    xyz_path=reactant_xyz,
    ncores=NCORES,
    sk_list=SIMPLE_KEYWORDS,
)

r_opt_out = run_calc(r_opt_calc)
check_and_parse_output(r_opt_out)

p_opt_calc = setup_opt_calc(
    name="final_complex",
    xyz_path=product_xyz,
    ncores=NCORES,
    sk_list=SIMPLE_KEYWORDS,
)

p_opt_out = run_calc(p_opt_calc)
check_and_parse_output(p_opt_out)

optimized_reactant_xyz = GEOM_DIR / "initial_complex" / "initial_complex.xyz"
optimized_product_xyz = GEOM_DIR / "final_complex" / "final_complex.xyz"

neb_ts_calc = setup_neb_calc(
    basename="ts",
    reactant_xyz=optimized_reactant_xyz,
    product_xyz=optimized_product_xyz,
    ncores=NCORES,
    sk_list=SIMPLE_KEYWORDS,
)

neb_ts_out = run_calc(neb_ts_calc)
check_and_parse_output(neb_ts_out)

ts_xyz = NEB_DIR / "ts.xyz"

irc_calc = setup_irc_calc(
    basename="irc",
    ts_xyz=ts_xyz,
    ncores=NCORES,
    sk_list=SIMPLE_KEYWORDS,
)

irc_out = run_calc(irc_calc)
check_and_parse_output(irc_out)
