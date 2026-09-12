from pathlib import Path

from opi.input.simple_keywords import Task

from setups import read_structure, make_cn_ccl_constraints, setup_calc, run_calc, check_and_parse_output


PROJECT = Path(__file__).resolve().parents[1]

xyz_file = PROJECT / "structures" / "radical" / "D_seed.xyz"

structure = read_structure(
        xyz_file,
        charge=0,
        multiplicity=2,
)


geom = make_cn_ccl_constraints(
        cn_target=1.48,
        ccl_target=2.63,
)


workdir = (
    PROJECT
    / "calculations"
    / "radical"
    / "D_reproduce"
    / "constrained"
)

constrained_calc = setup_calc(
        basename="D_constrained",
        working_dir=workdir,
        structure=structure,
        ncores=4,
        calculation_keyword=Task.OPT,
        geom=geom,
)


constrained_output = run_calc(constrained_calc)
check_and_parse_output(constrained_output)
constrained_structure = constrained_output.get_structure()

if constrained_structure is None:
    raise RuntimeError("Could not obtain constrained structure")


opt_workdir = (
    PROJECT
    / "calculations"
    / "radical"
    / "D_reproduce"
    / "unconstrained"
)

opt_calc = setup_calc(
        basename="D_opt",
        working_dir=opt_workdir,
        structure=constrained_structure,
        calculation_keyword=Task.OPT,
        ncores=4,
)


opt_output = run_calc(opt_calc)
check_and_parse_output(opt_output)
opt_structure = opt_output.get_structure()

if opt_structure is None:
    raise RuntimeError("Could not obtain optimized D structure")


freq_workdir = (
    PROJECT
    / "calculations"
    / "radical"
    / "D_reproduce"
    / "freq"
)

freq_calc = setup_calc(
        basename="D_freq",
        working_dir=freq_workdir,
        structure=opt_structure,
        calculation_keyword=Task.FREQ,
        ncores=4,
)


freq_output = run_calc(freq_calc)
check_and_parse_output(freq_output)
freq_structure = freq_output.get_structure()

if opt_structure is None:
    raise RuntimeError("Could not obtain frequency D structure")
