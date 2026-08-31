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

print(Path.cwd())

PROJECT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT)

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

def plot_mo_diagram(
    output,
    out_path=None,
    title="MO Diagram",
    window_ev=0.5,
):
    mo_list = output.get_mos()["mo"]
    energies = [mo.orbitalenergy for mo in mo_list]
    occupations = [mo.occupancy for mo in mo_list]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title(title)
    ax.set_ylabel("Energy / eV")

    lumo_id = occupations.index(0)
    homo_energy = energies[lumo_id-1]
    lumo_energy = energies[lumo_id]
    
    prev_e = None
    max_x = 0
    for i, (e, occ) in enumerate(zip(energies, occupations)):
        if prev_e is not None and abs(e-prev_e) < 0.01:
            x += 2
        else:
            x = 0
        if x > max_x: max_x = x
        prev_e = e
        ax.hlines(e, x - 0.4, x + 0.4, color='k', linewidth=2)

        if occ == 2:
            ax.annotate("↑", (x, e), textcoords="offset points", xytext=(-10, -5), ha='center', fontsize=12)
            ax.annotate("↓", (x, e), textcoords="offset points", xytext=(10, -5), ha='center', fontsize=12)
        elif occ == 1:
            raise ValueError("This function does not support plotting of UHF type wavefunctions.")

        # MO index on the right, start counting at 0
        if e > homo_energy - 0.5 and e < lumo_energy + 0.5: 
            label = f"MO {i}"
            ax.text(x + 0.6, e, label, va='center', fontsize=10)
        
    ax.set_xlim(-1, max_x + 1.25)
    ax.set_ylim(homo_energy - window_ev, lumo_energy + window_ev)
    ax.set_xticks([])
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5)
    plt.tight_layout()

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"saved {out_path}")

    return fig, ax

def save_mo_views_html(
    output,
    xyz_data,
    plot_list,
    out_path,
    resolution=40,
    gbw_type="gbw",
    isoval=0.005,
    width=300,
    height=300,
    title="Molecular orbitals",
):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    html_blocks = []

    for mo_index in plot_list:
        cube_output = output.plot_mo(
            mo_index,
            resolution=resolution,
            gbw_type=gbw_type,
        )
        cube_data = cube_output.cube

        view = py3Dmol.view(width=width, height=height)
        view.addModel(xyz_data, "xyz")
        view.setStyle({
            "stick": {"radius": 0.1},
            "sphere": {"scale": 0.2},
        })
        view.addVolumetricData(
            cube_data,
            "cube",
            {"isoval": isoval, "color": "blue", "opacity": 0.8},
        )
        view.addVolumetricData(
            cube_data,
            "cube",
            {"isoval": -isoval, "color": "red", "opacity": 0.8},
        )
        view.zoomTo()

        viewer_html = view._make_html()

        html_blocks.append(f"""
        <div class="viewer">
          <h3>MO {mo_index}</h3>
          {viewer_html}
        </div>
        """)

    page = f"""
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, sans-serif;
      margin: 2rem;
      background: #111;
      color: #eee;
    }}
    .grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 1rem;
      align-items: flex-start;
    }}
    .viewer {{
      background: #1b1b1b;
      border: 1px solid #333;
      border-radius: 12px;
      padding: 1rem;
    }}
    h1, h3 {{
      margin-top: 0;
    }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <p>gbw type: <code>{html.escape(gbw_type)}</code>, isoval: <code>{isoval}</code></p>
  <div class="grid">
    {''.join(html_blocks)}
  </div>
</body>
</html>
"""

    out_path.write_text(page, encoding="utf-8")
    print(f"saved {out_path}")
    return out_path

workdir = calc_folder("4nitrochlorobiphenyl_fast", clean=True)
workdir

# > initializing calculator
calc = Calculator(basename="job", working_dir=workdir)

# > adding geometry of water from smiles
calc.structure = Structure.from_xyz(XYZ_DIR / "4nitrochlorobiphenyl.xyz", charge=-1, multiplicity=1)

# > setting method Hartree-Fock Hartree-Fock with basis set def2-TZVP
calc.input.add_simple_keywords(
        Sqm.NATIVE_GFN2_XTB,
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

mo_list = output.get_mos()["mo"]
energies = [mo.orbitalenergy for mo in mo_list]
occupations = [mo.occupancy for mo in mo_list]   

# MO energy diagram
plot_mo_diagram(
    output,
    out_path=FIG_DIR / "4_nitrochlorobiphenyl_mo_diagram.png",
    title="Nitrochlorobiphenyl MO diagram",
)

# > Plot orbitals 
nel, _ = output.get_nelectrons()
homo_index = nel // 2 - 1
plot_list = [homo_index-2,homo_index-1,homo_index,homo_index+1]
xyz_data = (XYZ_DIR / "4nitrochlorobiphenyl.xyz").read_text()

save_mo_views_html(
    output=output,
    xyz_data=xyz_data,
    plot_list=plot_list,
    out_path=VIEW_DIR / "4nitrochlorobiphenyl_mos_gbw.html",
    resolution=40,
    gbw_type="gbw",
)

save_mo_views_html(
    output=output,
    xyz_data=xyz_data,
    plot_list=plot_list,
    out_path=VIEW_DIR / "4nitrochlorobiphenyl_mos_loc.html",
    resolution=40,
    gbw_type="loc",
)
