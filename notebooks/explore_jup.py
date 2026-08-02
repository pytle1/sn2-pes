# %%
from pathlib import Path
import os
import opi
from opi.input.structures.structure import Structure
import py3Dmol
import sys, pkgutil, importlib.metadata as md

PROJECT = Path.cwd()

SAVE_DIR = PROJECT / "results"
FIG_DIR = PROJECT / "figures"
CALC_DIR = PROJECT / "calculations"
VIEW_DIR = PROJECT / "figures" / "molecules"
XYZ_DIR = PROJECT / "molecules"

SAVE_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)
CALC_DIR.mkdir(exist_ok=True)
VIEW_DIR.mkdir(parents=True, exist_ok=True)
XYZ_DIR.mkdir(exist_ok=True)

print("Project:", PROJECT)
print("Current dir:", Path.cwd())

# > Define cartesian coordinates in angstrom as python string 
water_data = """\
3

O      0.00000   -0.00000    0.00000
H      0.00000    0.96899    0.00000
H      0.93966   -0.23409    0.03434\n
"""

# > Visualize the input structure

name = "4nitrobiphenyl"
name2 = name + ".xyz"

xyz = Path(XYZ_DIR / name2).read_text()

xyz
print(xyz)

def show_xyz(name):
    name2 = name + ".xyz"
    name3 = name + ".html"
    xyz = Path(XYZ_DIR / name2).read_text()

    view = py3Dmol.view(width=1000, height=700)
    view.addModel(xyz, 'xyz')
    view.setStyle({}, {'stick': {}})
    view.addPropertyLabels(
        "index",
        "",
        {"fontSize": 14}
    )
    view.zoomTo()

    html = view._make_html()

    out = VIEW_DIR / name3
    out.write_text(html, encoding="utf-8")

show_xyz("4nitrodiphenyl")

with open(SAVE_DIR / "4nitrobiphenyl_view.html", "w") as f:
    f.write(xyz)

# read structure into object

structure = Structure.from_xyz(SAVE_DIR / "water.xyz")
