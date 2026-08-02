from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, Scf, AuxBasisSet, Approximation, Dft, Task, OutputControl, Grid
from opi.input.structures.structure import Structure

import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np

def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 4) -> Calculator:
    # > Set up a Calculator object, the basename for the ORCA calculation is set to LED
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > Assign structure to calculator
    calc.structure = structure

    # define a level of theory: wB97M-V/def2-TZVPPD with simple
    sk_list = [
        Dft.WB97M_V, # > wB97M-V method
        BasisSet.DEF2_TZVPD, # > def2-TZVPD basis set
        AuxBasisSet.DEF2_J, # > auxiliary basis set
        Scf.TIGHTSCF, # > tight SCF
        Approximation.RIJCOSX, # > Use RIJ and COSX
        Grid.DEFGRID3, # > Large grid for exchange-correlation and COSX
        Task.ENGRAD, # > Calculate energy and gradient
        OutputControl.PRINTGAP, # > Print HOMO-LUMO gap
    ]

    # use simple keywords in calculator
    calc.input.add_simple_keywords(*sk_list)

    # Define number of CPUs for the calcualtion
    calc.input.ncores = ncores # > 4 CPUs for this ORCA run

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

def process_output(output: Output) -> dict:
    # > Conversion factor for Hartree -> eV
    h_to_eV = 27.2114
    # > Conversion factor for Bohr -> Angstrom
    bohr_to_A = 0.529177

    # > Energy and gradient information
    # > Total energy in eV
    Etot = output.results_properties.geometries[0].single_point_data.finalenergy*h_to_eV
    # > VV10 energy part in eV
    # > Note that VV10 is self-consistent in gradient computations and thus ecnl is used
    E_VV10 = output.results_properties.geometries[0].dft_energy.ecnl*h_to_eV
    # > Gradient in eV/Bohr (list of size n_atoms containing one list per atom containing three entries)
    Grad = output.results_properties.geometries[0].nuclear_gradient[0].grad
    Grad = [[x * (h_to_eV/bohr_to_A) for x in nuc_grad] for nuc_grad in Grad]

    # > General information
    # > Charge of molecule (float value)
    tot_chrg = output.results_properties.calculation_info.charge
    # > Spin of molecule (float value)
    tot_spin = output.results_properties.calculation_info.mult
    # > Number of atoms (float value)
    n_at = output.results_properties.calculation_info.numofatoms
    # > Number of electrons (float value)
    n_el = output.results_properties.calculation_info.numofelectrons
    # > Number of basis functions (float value)
    n_bf = output.results_properties.calculation_info.numofbasisfuncts
    # > Number of alpha electrons (float value)
    n_el_alp = output.results_properties.geometries[0].dft_energy.nalphael
    # > Number of beta electrons (float value)
    n_el_beta = output.results_properties.geometries[0].dft_energy.nbetael
    
    # > Population analysis
    # > Loewdin charges (list of size n_atoms containing lists with 3 entries)
    chrgs_lowedin = output.results_properties.geometries[0].loewdin_population_analysis[0].atomiccharges
    # > Mulliken charges (list of size n_atoms containing lists with 3 entries)
    chrgs_mulliken = output.results_properties.geometries[0].mulliken_population_analysis[0].atomiccharges
    # > Mayer bond order (list of size of the number of bonds)
    bond_order = output.results_properties.geometries[0].mayer_population_analysis[0].bondorders
    # > Mayer's total valence (list of size n_atoms)
    va_mayer = output.results_properties.geometries[0].mayer_population_analysis[0].va

    # > Dipole moments
    # > Total dipole moments (list with 3 entries: x, y, z)
    dip_tot = output.results_properties.geometries[0].dipole_moment[0].dipoletotal
    # > Magnitude of dipole (float value)
    dip_mag = output.results_properties.geometries[0].dipole_moment[0].dipolemagnitude
    # > Electronic dipole contributions (list with 3 entries: x, y, z)
    dip_el = output.results_properties.geometries[0].dipole_moment[0].dipoleeleccontrib
    # > Nuclear dipole contributions (list with 3 entries: x, y, z)
    dip_nu = output.results_properties.geometries[0].dipole_moment[0].dipolenuccontrib

    # > MO handling
    mos = output.results_gbw.molecule.molecularorbitals.mos # > saves a list of MOs with size n_atoms
    mo_energies = []
    mo_occupancy = []
    for mo in mos:
        mo_energies.append(mo.orbitalenergy) # > Extract the MO energies
        mo_occupancy.append(mo.occupancy) # > Extract the occupation numbers
    mo_energies = np.array(mo_energies) # > Make float list
    mo_occupancy = np.array(mo_occupancy) # > Make float list
    if np.all(np.isin(mo_occupancy, [1.0, 0.0])):
        print("Unrestricted calculation detected.")
    elif np.all(np.isin(mo_occupancy, [2.0, 0.0])):
        print("Restricted calculation detected.")

    HOMO_energy = np.max(mo_energies[(mo_occupancy == 1.0) | (mo_occupancy == 2.0)]) # > Detect HOMO energy
    HOMO_energy = HOMO_energy * h_to_eV
    LUMO_energy = np.min(mo_energies[mo_occupancy == 0.0]) # > Detect LUMO energy
    LUMO_energy = LUMO_energy * h_to_eV

    summary_data = {
    "Total Energy [eV]": [Etot],
    "Gradient [ev/A]": [str(Grad[0][0])[:10]+', '+str(Grad[1][0])[:10]+', '+str(Grad[2][0])[:10]],
    "VV10 Energy [eV]": [E_VV10],
    "Dipole Magnitude": [dip_mag],
    "Total Charge": [tot_chrg],
    "Spin Multiplicity": [tot_spin],
    "Number of Atoms": [n_at],
    "Number of Electrons": [n_el],
    "Alpha Electrons": [n_el_alp],
    "Beta Electrons": [n_el_beta],
    "Basis Functions": [n_bf],
    "HOMO Energy [eV]": [HOMO_energy],
    "LUMO Energy [eV]": [LUMO_energy],
    }
    
    return summary_data


PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)
print(Path.cwd())

from paths import CALC_DIR, FIG_DIR, XYZ_DIR, calc_folder

workdir = calc_folder("ml_properties", clean=True)
workdir

structure = Structure.from_xyz(XYZ_DIR / "4nitrobiphenyl.xyz", charge=0, multiplicity=1)

calc = setup_calc("single_point", working_dir=workdir, structure=structure)
output = run_calc(calc)
check_and_parse_output(output)
summary_data = process_output(output)
df_summary = pd.DataFrame(summary_data, index=[""])
df_summary.T
